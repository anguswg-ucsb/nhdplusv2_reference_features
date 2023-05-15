find_file_path = function(VPU, files, new_dir){
  f = grep(VPU, basename(files))
  tmp01 = files[f]
  tmp02 = paste0(new_dir,  basename(files[f]))
  
  if(!file.exists(tmp02)){
    log_info('\tCopying ', VPU, " to reference.")
    file.copy(grep(tmp01, files, value = TRUE), tmp02)
  }
  
  tmp02
  
}

drop_zm_dims <- function(shp) {
  
  # Drop the geometry column and rename column names to upper case
  tmp <- st_drop_geometry(shp)
  names(tmp) <- toupper(names(tmp))
  
  # drop Z/M dimensions and bring back original geometry
  shp <- st_zm(
    st_as_sf(
      cbind(tmp, st_geometry(shp))
    )
  )
  
  return(shp)
  
}

# quickly check and validate invalid geometries only
fast_validity_check <- function(x){
  
  bool    = st_is_valid(x)
  valid   = filter(x, bool)
  invalid = st_make_valid(filter(x, !bool)) %>% 
    st_cast("POLYGON")
  
  return(bind_rows(valid, invalid))
  
}

clean_geometry2 <- function(catchments,
                            ID = "ID",
                            keep = NULL,
                            crs = 5070,
                            grid = .0009,
                            sys = NULL) {
  
  # keep an original count of # rows in catchments
  MASTER_COUNT = nrow(catchments)
  # use system mapshaper or not
  if(Sys.getenv("TURN_OFF_SYS_MAPSHAPER") == "YUP")  { sys <- FALSE }
  
  if(is.null(sys)) {
    sys <- FALSE
    try(sys <- is.character(check_sys_mapshaper(verbose = FALSE)))
  }
  
  # set crs variable to crs of catchments
  if(!is.null(crs)){  crs = st_crs(catchments)  }
  
  catchments = lwgeom::st_snap_to_grid(
    st_transform(catchments, 5070), 
    size = grid
  )
  
  # cast MPs to POLYGONS and add featureid count column
  polygons = suppressWarnings({
    catchments %>%
      st_cast("POLYGON") %>%
      fast_validity_check() %>%
      add_count(!!sym(ID)) %>%
      mutate(areasqkm = add_areasqkm(.), tmpID = 1:n()) %>%
      rename_geometry("geometry")
  })
  
  if(any(st_geometry_type(polygons) != "POLYGON") | nrow(polygons) != MASTER_COUNT){
    
    # separate polygons with more than 1 feature counts
    extra_parts = filter(polygons, n != 1)
    
    # dissolve, and explode if necessary
    try(
      extra_parts <- ms_explode(
        ms_dissolve(
          extra_parts, ID, copy_fields = names(extra_parts)
        )
      ),
      silent = TRUE
    )
    
    # recalculate area
    extra_parts <- 
      extra_parts %>% 
      group_by(.data[[ID]]) %>%
      mutate(
        areasqkm = add_areasqkm(extra_parts), 
        newID    = row_number(desc(areasqkm))
      ) %>% 
      ungroup()
    
    # get the biggest parts by area in each catchment and bind with rest of good_to_go catchments
    main_parts <-
      extra_parts %>%
      group_by(.data[[ID]]) %>%
      slice_max(areasqkm, with_ties = FALSE) %>%
      ungroup()
    
    small_parts <-
      extra_parts %>%
      filter(newID != main_parts$newID)
    
    if(!sum(nrow(main_parts)) + nrow(filter(polygons, n == 1)) == MASTER_COUNT){
      
      stop() 
      
    }
    
    main_parts =  bind_rows(main_parts, filter(polygons, n == 1))
    
    if(nrow(small_parts) > 0){
      # dissolve, and explode if necessary
      small_parts <- tryCatch(
        ms_explode(
          ms_dissolve(small_parts, ID, copy_fields = names(small_parts))
        ), error = function(e){ 
          NULL
        }, 
        warning = function(w){
          NULL
        }
      )
      # add area
      small_parts =  mutate(
        small_parts,
        areasqkm = add_areasqkm(small_parts),
        newID    = 1:n()
      ) %>%
        select(newID)
      # get the intersection between big parts and small parts and pull out the LINESTRINGs
      out = tryCatch({
        suppressWarnings({
          st_intersection(small_parts, st_make_valid(main_parts)) %>%
            st_collection_extract("LINESTRING")
        })
      }, error = function(e) {
        NULL
      })
      
      ints =  
        out %>%
        mutate(l = st_length(.)) %>%
        group_by(newID) %>%
        slice_max(l, with_ties = FALSE) %>%
        ungroup()
      
      tj <- 
        right_join(
          small_parts,
          select(st_drop_geometry(ints), featureid, newID),
          by = "newID"
        ) %>%
        bind_rows(main_parts) %>%
        select(-areasqkm, -tmpID, -newID) %>%
        group_by(featureid) %>%
        mutate(n = n()) %>%
        ungroup() %>%
        rename_geometry("geometry")
      
      in_cat <-
        union_polygons(
          filter(tj, .data$n > 1),
          "featureid"
        ) %>%
        bind_rows(
          select(
            filter(tj, .data$n == 1),
            "featureid")
        ) %>%
        mutate(tmpID = 1:n()) %>%
        fast_validity_check()
      
    } else {
      
      in_cat = fast_validity_check(main_parts)
      
    }
    
    if(all(st_is_valid(in_cat)) & all(st_geometry_type(in_cat) == "POLYGON")){
      
      if(!is.null(keep)){ in_cat = simplify_process(in_cat, keep, sys) }
      
    } else {
      
      warning ("Invalid geometries found.", call. = FALSE)
      
    }
    
    return(
      
      mutate(in_cat, areasqkm = add_areasqkm(in_cat)) %>%
        st_transform(crs) %>%
        select("{ID}" := ID,areasqkm)  %>%
        left_join(st_drop_geometry(select(catchments, -areasqkm)), by = ID)
      
    )
  } else {
    return(
      polygons
    )
  }
}
# 
clean_geometry2 <- function(catchments,
                            ID = "ID",
                            keep = NULL,
                            crs = 5070,
                            grid = .0009,
                            sys = NULL
) {
  
  # keep an original count of # rows in catchments
  MASTER_COUNT = nrow(catchments)
  
  # use system mapshaper or not
  if(Sys.getenv("TURN_OFF_SYS_MAPSHAPER") == "YUP")  {
    
    sys <- FALSE
    
  }
  
  if(is.null(sys)) {
    sys <- FALSE
    try(sys <- is.character(check_sys_mapshaper(verbose = FALSE)))
  }
  
  # set crs variable to crs of catchments
  if(!is.null(crs)){
    
    crs = st_crs(catchments)
    
  }
  
  catchments = lwgeom::st_snap_to_grid(
    st_transform(catchments, 5070),
    grid
  )
  
  # quickly check and validate invalid geometries only
  fast_validity_check <- function(x){
    
    bool    = st_is_valid(x)
    valid   = filter(x, bool)
    invalid = st_make_valid(filter(x, !bool)) %>%
      st_cast("POLYGON")
    
    return(bind_rows(valid, invalid))
    
  }
  
  
  message("\t--- Casting to Polygon")
  
  # cast MPs to POLYGONS and add featureid count column
  pgs = suppressWarnings({
    catchments %>%
      st_cast("POLYGON") %>%
      fast_validity_check() %>%
      add_count(!!sym(ID)) %>%
      mutate(areasqkm = add_areasqkm(.), tmpID = 1:n()) %>%
      rename_geometry("geometry")
  })
  
  
  if(any(st_geometry_type(pgs) != "POLYGON") | nrow(pgs) != MASTER_COUNT){
    
    message("\t--- Refining Parts")
    
    # separate polygons with more than 1 feature counts
    parts = filter(pgs, n != 1)
    
    # dissolve, and explode if necessary
    try(
      parts <- ms_explode(
        ms_dissolve(parts, ID, copy_fields = names(parts))
      ),
      silent = TRUE
    )
    
    # recalculate area
    parts <- mutate(parts, areasqkm = add_areasqkm(parts))
    
    parts <-
      parts %>%
      dplyr::group_by(.data[[ID]]) %>%
      dplyr::mutate(
        newID = dplyr::row_number(dplyr::desc(areasqkm))
      ) %>%
      dplyr::ungroup()
    
    # get the biggest parts by area in each catchment and bind with rest of good_to_go catchments
    big_parts <-
      parts %>%
      group_by(.data[[ID]]) %>%
      slice_max(areasqkm, with_ties = FALSE) %>%
      # slice_min(newID, with_ties = FALSE) %>%
      ungroup() %>%
      bind_rows(filter(pgs, n == 1))
    
    small_parts <-
      parts %>%
      group_by(.data[[ID]]) %>%
      dplyr::filter(newID != 1) %>%
      ungroup()
    
    # dissolve, and explode if necessary
    try(
      # small_parts <-  ms_explode(small_parts),
      small_parts <-  ms_explode(
        ms_dissolve(small_parts, ID, copy_fields = names(small_parts))
      ),
      silent = F
    )
    
    # add area
    small_parts =  mutate(
      small_parts,
      areasqkm = add_areasqkm(small_parts),
      newID    = 1:n()
    ) %>%
      dplyr::select(newID)
    
    message("\t--- Computing intersections")
    
    # get the intersection between big parts and small parts and pull out the LINESTRINGs
    out = tryCatch({
      suppressWarnings({
        
        st_intersection(small_parts, st_make_valid(big_parts)) %>%
          st_collection_extract("LINESTRING")
      })
    }, error = function(e) {
      NULL
    })
    
    ints =
      out %>%
      mutate(l = st_length(.)) %>%
      group_by(.data$newID) %>%
      slice_max(.data$l, with_ties = FALSE) %>%
      ungroup()
    
    tj = right_join(
      small_parts,
      select(st_drop_geometry(ints), featureid, newID),
      by = "newID"
    ) %>%
      bind_rows(big_parts) %>%
      select(-areasqkm, -tmpID, -newID) %>%
      group_by(featureid) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      rename_geometry('geometry')
    
    message("\t--- Dissolving catchment fragments")
    
    in_cat <-
      union_polygons(
        filter(tj, .data$n > 1),
        'featureid'
      ) %>%
      bind_rows(
        select(
          filter(tj, .data$n == 1),
          'featureid')
      ) %>%
      mutate(tmpID = 1:n()) %>%
      fast_validity_check()
    
    
    if(all(st_is_valid(in_cat)) & all(st_geometry_type(in_cat) == "POLYGON")){
      
      message("\t--- Simplifying")
      
      if(!is.null(keep)){ in_cat = simplify_process(in_cat, keep, sys) }
      
    } else {
      warning ("Invalid geometries found.", call. = FALSE)
    }
    
    return(
      mutate(in_cat, areasqkm = add_areasqkm(in_cat)) %>%
        st_transform(crs) %>%
        select("{ID}" := ID,areasqkm)  %>%
        left_join(
          st_drop_geometry(select(catchments, -areasqkm)),
          by = ID
        )
    )
    
  }
}

simplify_process = function(catchments, keep, sys){
  
  # simplify catchments
  catchments =  ms_simplify(catchments, keep = keep, keep_shapes = TRUE, sys = sys)
  
  # mark valid/invalid geoms
  bool = st_is_valid(catchments)
  
  # make invalid geoms valid
  invalids = st_make_valid(filter(catchments, !bool))
  
  # if not all polygons get returned, try different simplification keep value
  if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){
    
    warning("Invalid geometries found. Trying new keep of:",  keep + ((1-keep) / 2) , call. = FALSE)
    
    # try simplification again
    catchments = ms_simplify(catchments, keep =  keep + ((1-keep) / 2), keep_shapes = TRUE, sys = sys)
    
    # mark valid/invalid geoms
    bool = st_is_valid(catchments)
    
    # make invalid geoms valid
    invalids = st_make_valid(filter(catchments, !bool))
    
    # if catchments still containing non POLYGON geometries, return original data
    if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){ 
      
      warning("Invalid geometries found. Original catchments returned." , call. = FALSE) 
      
      return(catchments)
      
    } else {
      
      # combine corrected invalids with valids and recalc area
      return(
        mutate(
          bind_rows(invalids, filter(catchments, bool)), 
          areasqkm = add_areasqkm(.)
        ) 
      )
    }
    
  } else {
    
    # combine corrected invalids with valids and recalc area
    return(
      mutate(
        bind_rows(invalids, filter(catchments, bool)), 
        areasqkm = add_areasqkm(.)
      ) 
    )
  }
}
# 
# simplify_process = function(catch, keep, sys){
#   
#   catch =  ms_simplify(catch, keep = keep, keep_shapes = TRUE, sys = sys)
#   
#   bool = st_is_valid(catch)
#   invalids = st_make_valid(filter(catch, !bool))
#   
#   if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){ 
#     warning("Invalid geometries found. Trying new keep of:",  keep + ((1-keep) / 2) , call. = FALSE) 
#     catch =  ms_simplify(catch, keep =  keep + ((1-keep) / 2), keep_shapes = TRUE, sys = sys)
#     bool = st_is_valid(catch)
#     invalids = st_make_valid(filter(catch, !bool))
#     
#     if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){ 
#       warning("Invalid geometries found. Original catchments returned." , call. = FALSE) 
#       return(catch)
#     } else {
#       bind_rows(invalids, filter(catch, bool)) %>% 
#         mutate(areasqkm = add_areasqkm(.)) 
#     }
#   } else {
#     bind_rows(invalids, filter(catch, bool)) %>% 
#       mutate(areasqkm = add_areasqkm(.)) 
#   }
# }
# 
# simplify_process = function(catch, keep, sys){
#   
#   catch =  ms_simplify(catch, keep = keep, keep_shapes = TRUE, sys = sys)
#   
#   bool = st_is_valid(catch)
#   invalids = st_make_valid(filter(catch, !bool))
#   
#   if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){ 
#     
#     warning("Invalid geometries found. Trying new keep of:",  keep + ((1-keep) / 2) , call. = FALSE) 
#     
#     catch =  ms_simplify(catch, keep =  keep + ((1-keep) / 2), keep_shapes = TRUE, sys = sys)
#     
#     bool = st_is_valid(catch)
#     invalids = st_make_valid(filter(catch, !bool))
#     
#     if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){ 
#       warning("Invalid geometries found. Original catchments returned." , call. = FALSE) 
#       return(catch)
#     } else {
#       bind_rows(invalids, filter(catch, bool)) %>% 
#         mutate(areasqkm = add_areasqkm(.)) 
#     }
#   } else {
#     bind_rows(invalids, filter(catch, bool)) %>% 
#       mutate(areasqkm = add_areasqkm(.)) 
#   }
# }
# simplify_process = function(catchments, keep, sys){
#   
#   catchments =  ms_simplify(catchments, keep = keep, keep_shapes = TRUE, sys = sys)
#   #suppressWarnings({ st_set_crs( read_sf(out_tmp[1]), 5070) })
#   
#   bool = st_is_valid(catchments)
#   invalids = st_make_valid(filter(catchments, !bool))
#   
#   if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){
#     return(NULL)
#   } else {
#     bind_rows(invalids, filter(catchments, bool)) %>% 
#       mutate(areasqkm = add_areasqkm(.)) %>% 
#       select(-tmpID)
#   }
# }
# simplify_process <- function(catchments = out, out_geojson, out_tmp, num){
#   
#   unlink(out_geojson[1])
#   unlink(out_tmp[1])
#   
#   write_sf(catchments, out_geojson[1])
# 
#   if (Sys.info()['sysname'] == "Windows") {
#     
#     system( paste0(
#       "mapshaper ", out_geojson ,' -simplify ', num, 
#       '% keep-shapes -o ',  "format=geojson ",  out_tmp)
#     )
#   
#   } else{ 
#     
#     system(paste0('node  --max-old-space-size=16000 `which mapshaper` ',
#                   out_geojson, ' -simplify ',
#                   num, '% keep-shapes -o ',
#                   out_tmp))
#     }
#   
# 
#   catchments = suppressWarnings({ st_set_crs( read_sf(out_tmp[1]), 5070) })
#   
#   bool = st_is_valid(catchments)
#   invalids = st_make_valid(filter(catchments, !bool))
#   
#   
#   if(nrow(filter(invalids, st_geometry_type(invalids) != "POLYGON")) > 0){
#     return(NULL)
#   } else {
#     bind_rows(invalids, filter(catchments, bool)) %>% 
#       mutate(areasqkm = add_areasqkm(.)) %>% 
#       select(-tmpID)
#   }
# }

explode_collections <- function(catch) {
  
  # separate out GEOMETRYCOLLECTIONs, extract polygons and bind together
  cat_gc <-
    catch %>%
    dplyr::filter(st_geometry_type(catch) == "GEOMETRYCOLLECTION") %>% 
    sf::st_collection_extract() %>% 
    sf::st_cast("POLYGON") %>% 
    dplyr::group_by(featureid, gridcode, sourcefc) %>% 
    dplyr::summarise() %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(areasqkm = add_areasqkm(.)) %>% 
    dplyr::relocate(featureid, areasqkm, gridcode, sourcefc, geometry)
  
  # filter out GEOMETRYCOLLECTIONs and bind with cat_gc
  catch <- 
    catch %>% 
    dplyr::filter(st_geometry_type(catch) != "GEOMETRYCOLLECTION") %>% 
    dplyr::bind_rows(cat_gc)
  
  # return the updated catchments object
  return(catch)
}
# test

# get_nhdplus(comid = 3645, realization = "all") %>%
#   mapview::mapview()

# get_nhdplus(comid = 3981, realization = "all") %>%
#   mapview::mapview()
# i = 2
# flowlines = st_transform(read_sf(flowpath_files[i]), 5070)
# ID = "featureid"

# Fixes issues with MULTIPOLYGON geoms in a sf catchment/flowline datasets.
# First identifies the problematic "MULTIPOLYGON" geoms and converts them to "POLYGON" geometries. 
# Then intersects each "POLYGON" geometry with the flowlines and IDs the longest "LINESTRING" geometry after intersection.
# If a "POLYGON" geometry fails to intersect with any flowlines or the intersection result is a "POINT", "MULTIPOINT", or "POLYGON", the function logs a warning message. 
# Otherwise, it dissolves the "POLYGON" geometry with the longest "LINESTRING" geometry and updates the dataset. 
# Last step in function checks if the number and type of geoms in the updated dataset matches those in the original dataset 
# and if the check passes, the updated dataset is returned, otherwise function stops
fix_mp_issues = function(
    catch, 
    flowlines,
    ID = "ID"
){
  
  # catch = catchments_valid
  # ID = "featureid"
  
  mps = filter(catch, st_geometry_type(catch) == "MULTIPOLYGON")
  mps =  suppressWarnings({ st_cast(mps, "POLYGON") }) %>%
    filter(!duplicated(.))
  
  mps <- 
    mps %>%
    rmapshaper::ms_dissolve("featureid", copy_fields = names(mps)) 
  
  
  good_to_go = 
    catch %>% 
    filter(st_geometry_type(catch) == "POLYGON") %>%
    nhdplusTools::rename_geometry("geometry") %>% 
    dplyr::bind_rows(filter(mps, st_geometry_type(catch) == "POLYGON"))
  
  new = list()
  
  if(nrow(mps) > 1){
    
    tmp  = mutate(mps, areasqkm = add_areasqkm(mps))
    
    uids = unique(tmp[[ID]])
    message("Exploding ", length(uids), " unique IDs.")
    
    if (length(uids) > 0) {
      for (i in 1:length(uids)) {
        
        i = 15
        
        here = 
          filter(tmp, get(ID) == uids[i])  %>%
          mutate(newID = 1:n())
        mapview::mapview(here)
        plot(here$geometry)
        # plot(here$geometry[2])
        
        here %>%
          rmapshaper::ms_dissolve("featureid", copy_fields = names(here)) %>%
          .$geometry%>%
          plot()
        
        imap =
          here %>% 
          st_intersects(
            filter(flowlines, COMID == here[[ID]][1])
          ) %>%
          lengths()
        
        if(sum(imap == 1) > 1){
          
          ll = vector()
          
          suppressWarnings({
            for(z in 1:nrow(here)){
              tmp_l = st_intersection(here[z,], filter(flowlines, COMID == here[[ID]][1])) %>%
                st_length()
              
              if(length(tmp_l) == 0){
                ll[z] =  0
              } else {
                ll[z] =  tmp_l
              }
            }
          })
          
          imap = rep(0, length(ll))
          imap[which.max(ll)] = 1
        } else if(sum(imap == 1) < 1){
          
          imap = rep(0, nrow(here))
          imap[1] = 1
        }
        
        new[[i]] = filter(here, imap == 1)
        
        to_dissolve = filter(here, imap != 1)
        
        
        for(j in 1:nrow(to_dissolve)){
          # j = 1
          # filter good to go polygons if they touch any of the polygons that are to be dissolved
          tmap = st_filter(
            good_to_go,
            to_dissolve[j,],
            .predicate = st_touches
          )
          # mapview::mapview(tmap) + here + to_dissolve + opt
          # # plot(tmap$geometry)
          # tmap = st_filter(
          #   good_to_go,
          #   to_dissolve[j,],
          #   .predicate = st_intersects
          # )
          # 
          
          suppressWarnings({
            xx  <- sf::st_snap(tmap, to_dissolve[j,], tolerance = 0.009)
            
            opt <- st_intersection(to_dissolve[j,], xx)
            
            gt  <- sf::st_geometry_type(opt)
            
            message(paste0("length == 1 check: ", length(gt) == 1 && grepl("POINT", gt)))
            message(paste0("length > 0 check: ", length(gt) > 0 && grepl("POINT", gt)))
            
            if (length(gt) == 1 && grepl("POINT", gt)) {
              
              message(paste0("--> length == ", length(gt)))
              opt <- sf::st_cast(opt, "LINESTRING")
              
            }  else if (length(gt) > 0 && all(grepl("POINT", gt))) {
              # } else if (length(gt) > 1 && grepl("POINT", gt)) {
              
              message(paste0("--> length == ", length(gt)))
              
              opt <- sf::st_cast(sf::st_cast(opt, "POINT"), "LINESTRING")
              
            }
            
            
          })
          
          opt <-
            opt %>%
            st_collection_extract("LINESTRING") %>%
            mutate(l = st_length(.)) %>%
            slice_max(.data$l, with_ties = FALSE)
          
          ind = which(
            good_to_go[[ID]] == opt[[paste0(ID, ".1")]]
          )
          
          # mapview::mapview(here, col.region = "red") + tmap + good_to_go[ind,] + opt + to_dissolve
          
          good_to_go$geometry[ind] = st_union(
            st_union(
              st_geometry(to_dissolve[j,])
            ),
            st_geometry(good_to_go[ind,])
          )
          
        }
        
        message(i)
      }
    }
    
    new_adj = 
      bind_rows(new) %>%
      mutate(newID = NULL)
    
    xx = bind_rows(good_to_go, new_adj)
    
  } else {
    xx = good_to_go
  }
  
  if(nrow(xx) == length(unique(catch[[ID]])) &
     sum(st_geometry_type(xx) == "POLYGON") == length(unique(catch[[ID]]))
  ) {
    
    return(xx)
    
  } else {
    stop()
  }
}


fix_mp_issues3 = function(
    catch, 
    flowlines,
    ID = "ID"
){
  
  mps2 = filter(catch, st_geometry_type(catch) == "MULTIPOLYGON")
  mps2 =  suppressWarnings({ st_cast(mps2, "POLYGON") }) %>%
    filter(!duplicated(.))
  
  good_to_go = 
    catch %>% 
    filter(st_geometry_type(catch) == "POLYGON") %>%
    nhdplusTools::rename_geometry("geometry")
  
  new = list()
  
  if(nrow(mps2) > 1){
    
    tmp  = mutate(mps2, areasqkm = add_areasqkm(mps2))
    
    uids = unique(tmp[[ID]])
    message("Exploding ", length(uids), " unique IDs.")
    # i = 4
    
    if (length(uids) > 0) {
      for (i in 1:length(uids)) {
        
        message(paste0(i, "/", length(uids)))
        
        here = 
          filter(tmp, get(ID) == uids[i])  %>%
          mutate(newID = 1:n())
        
        imap =
          here %>% 
          st_intersects(
            filter(flowlines, COMID == here[[ID]][1])
          ) %>%
          lengths()
        
        if(sum(imap == 1) > 1){
          
          ll = vector()
          
          suppressWarnings({
            for(z in 1:nrow(here)){
              tmp_l = st_intersection(here[z,], filter(flowlines, COMID == here[[ID]][1])) %>%
                st_length()
              
              if(length(tmp_l) == 0){
                ll[z] =  0
              } else {
                ll[z] =  tmp_l
              }
            }
          })
          
          imap = rep(0, length(ll))
          imap[which.max(ll)] = 1
          
        } else if(sum(imap == 1) < 1){
          
          imap = rep(0, nrow(here))
          imap[1] = 1
          
        }
        
        new[[i]]    = filter(here, imap == 1)
        
        to_dissolve = filter(here, imap != 1)
        
        for(j in 1:nrow(to_dissolve)){
          
          
          tmap = sf::st_filter(
            good_to_go,
            to_dissolve[j,],
            .predicate = st_intersects
          )
          
          message(paste0("# of neighbors: ", nrow(tmap)))
          
          suppressWarnings({
            
            opt <- st_intersection(to_dissolve[j,], tmap)
            gt  <- sf::st_geometry_type(opt)
            
            
            if (length(gt) == 1 && grepl("POINT", gt)) {
              
              message(paste0("--> length == ", length(gt)))
              
              opt <- sf::st_cast(opt, "LINESTRING")
              
            } else if (length(gt) > 1 & all(grepl("POINT", gt))) {
              
              opt <- sf::st_cast(sf::st_cast(opt, "POINT"), "LINESTRING")
              
            } 
            
          })
          
          opt <- st_cast(opt)
          
          if (!all(sf::st_geometry_type(opt) %in% c("POINT"))) {
            
            message(paste0("Discarding POINT geometry"))
            
            # opt <- 
            #   opt %>% 
            #   dplyr::filter(sf::st_geometry_type(opt) != "POINT")
            
            if (any(sf::st_geometry_type(opt) %in% c("GEOMETRYCOLLECTION"))) {
              
              message(paste0("GEOMETRYCOLLECTION geometry detected"))
              # opt %>% 
              # st_collection_extract("LINESTRING")
              
              opt <-
                opt %>% 
                dplyr::filter(sf::st_geometry_type(opt) != "POINT") %>% 
                sf::st_cast() %>% 
                dplyr::filter(sf::st_geometry_type(.) != "POINT") %>% 
                sf::st_cast("MULTILINESTRING") %>% 
                sf::st_cast("LINESTRING") %>% 
                dplyr::group_by(featureid.1) %>%
                dplyr::mutate(
                  l    = st_length(geometry),
                  suml = sum(l)
                ) %>% 
                dplyr::ungroup() %>%
                slice_max(.data$suml, with_ties = FALSE)
              
            } else if (any(sf::st_geometry_type(opt) %in% c("MULTIPOLYGON"))) {
              
              message(paste0("MULTIPOLYGON geometry detected"))
              
              opt <- 
                opt %>% 
                dplyr::filter(sf::st_geometry_type(opt) != "POINT") %>% 
                sf::st_cast("MULTILINESTRING") %>% 
                dplyr::group_by(featureid.1) %>%
                dplyr::mutate(
                  l    = st_length(geometry),
                  suml = sum(l)
                ) %>% 
                dplyr::ungroup() %>%
                slice_max(.data$suml, with_ties = FALSE)
              
            } else {
              
              message(paste0("POLYGON/LINESTRING geometry detected"))
              
              opt <- 
                opt %>% 
                dplyr::filter(sf::st_geometry_type(opt) != "POINT") %>%
                sf::st_cast("LINESTRING") %>% 
                # st_collection_extract("LINESTRING")
                mutate(l = st_length(.)) %>%
                slice_max(.data$l, with_ties = FALSE)
              
              
            }
            
          } else {
            
            message(paste0("ONLY POINT geometry detected"))
            
            opt <- 
              opt %>% 
              sf::st_cast("LINESTRING") %>% 
              mutate(l = st_length(.)) %>%
              slice_max(.data$l, with_ties = FALSE)
            
          }
          
          ind = which(
            good_to_go[[ID]] == opt[[paste0(ID, ".1")]]
          )
          
          good_to_go$geometry[ind] = st_union(
            st_union(
              st_geometry(to_dissolve[j,])
            ),
            st_geometry(good_to_go[ind,])
          )
          
        }
        
      }
    }
    
    new_adj = 
      bind_rows(new) %>%
      mutate(newID = NULL)
    
    xx = bind_rows(good_to_go, new_adj)
    
  } else {
    xx = good_to_go
  }
  
  if(nrow(xx) == length(unique(catch[[ID]])) &
     sum(st_geometry_type(xx) == "POLYGON") == length(unique(catch[[ID]]))
  ) {
    
    return(xx)
    
  } else {
    stop()
  }
}

fix_mp_issues2 = function(
    catch, 
    flowlines,
    ID = "ID"
){
  
  # catch = catchments
  # flowlines = st_transform(read_sf(flowpath_files[i]), 5070)
  # ID = "featureid"
  mps = filter(catch, st_geometry_type(catch) != "POLYGON")
  mps =  suppressWarnings({ st_cast(mps, "POLYGON") }) %>%
    filter(!duplicated(.))
  
  good_to_go = 
    catch %>% 
    filter(st_geometry_type(catch) == "POLYGON") %>%
    nhdplusTools::rename_geometry("geometry") 
  # dplyr::bind_rows(
  #   nhdplusTools::rename_geometry(
  #     st_cast(
  #       filter(catch, st_geometry_type(catch) == "GEOMETRYCOLLECTION"), 
  #       "POLYGON",
  #       do_split = F
  #     ), 
  #     "geometry"
  #   )
  # )
  
  # new = list()
  new2 = list()
  
  if(nrow(mps) > 1){
    
    tmp  = mutate(mps, areasqkm = add_areasqkm(mps))
    
    uids = unique(tmp[[ID]])
    message("Exploding ", length(uids), " unique IDs.")
    
    if (length(uids) > 0) {
      for (i in 1:length(uids)) {
        message(paste0("========================"))
        message(paste0("Index: ", i))
        message(paste0("UID: ", uids[i]))
        message(paste0("========================"))
        
        here = 
          filter(tmp, get(ID) == uids[i])  %>%
          mutate(newID = 1:n())
        
        imap =
          here %>% 
          st_intersects(
            filter(flowlines, COMID == here[[ID]][1])
          ) %>%
          lengths()
        
        if(sum(imap == 1) > 1){
          
          ll = vector()
          
          suppressWarnings({
            for(z in 1:nrow(here)){
              tmp_l = st_intersection(here[z,], filter(flowlines, COMID == here[[ID]][1])) %>%
                st_length()
              
              if(length(tmp_l) == 0){
                ll[z] =  0
              } else {
                ll[z] =  tmp_l
              }
            }
          })
          
          imap = rep(0, length(ll))
          imap[which.max(ll)] = 1
        } else if(sum(imap == 1) < 1){
          
          imap = rep(0, nrow(here))
          imap[1] = 1
        }
        
        # new[[i]] =  filter(here, imap == 1)
        
        to_dissolve = filter(here, imap != 1)
        
        # resolve internal boundaries and return a single POLYGON 
        mls_lst <- resolve_internal_bounds(poly = here)
        
        if (sf::st_geometry_type(mls_lst) == "MULTIPOLYGON") {
          
          message(paste0(
            "CATCH STILL MP, UNIONING\n---> INSERTING INTO 'good_to_go' geometry via 'st_touches'"
          ))
          
          for(j in 1:nrow(to_dissolve)){
            
            # filter good to go polygons if they touch any of the polygons that are to be dissolved
            tmap = st_filter(
              good_to_go, 
              to_dissolve[j,], 
              .predicate = st_touches
            )
            
            suppressWarnings({
              
              opt <- st_intersection(to_dissolve[j,], tmap)
              gt  <- sf::st_geometry_type(opt)
              
              message(paste0("length == 1 check: ", length(gt) == 1 && grepl("POINT", gt)))
              message(paste0("length > 0 check: ", length(gt) > 0 && grepl("POINT", gt)))
              
              if (length(gt) == 1 && grepl("POINT", gt)) {
                
                message(paste0("--> length == ", length(gt)))
                opt <- sf::st_cast(opt, "LINESTRING")
                
              } else if (length(gt) > 0 && all(grepl("POINT", gt))) {
                message(paste0("--> length == ", length(gt)))
                
                opt <- sf::st_cast(sf::st_cast(opt, "POINT"), "LINESTRING")
                
              }
              
              
            })
            
            opt <-
              opt %>%
              st_collection_extract("LINESTRING") %>%
              mutate(l = st_length(.)) %>%
              slice_max(.data$l, with_ties = FALSE)
            
            ind = which(
              good_to_go[[ID]] == opt[[paste0(ID, ".1")]]
            )
            
            good_to_go$geometry[ind] = st_union(
              st_union(
                st_geometry(to_dissolve[j,])
              ),
              st_geometry(good_to_go[ind,])
            )
            
            here <- resolve_internal_bounds(
              sf::st_difference(here,   
                                st_union(
                                  st_geometry(to_dissolve[j,])
                                )
              )
            )
            
            
          }
          
          new2[[i]] <- here
          
        } else {
          
          new2[[i]] <- mls_lst
          
        }
        
        
      }
      
    }
    
    new_adj = 
      bind_rows(new2) %>%
      mutate(newID = NULL)
    
    xx = bind_rows(good_to_go, new_adj)
    
  } else {
    
    message(paste0("ALL GOOD TO GO!"))
    
    xx = good_to_go
  }
  
  if(nrow(xx) == length(unique(catch[[ID]])) &
     sum(st_geometry_type(xx) == "POLYGON") == length(unique(catch[[ID]]))
  ) {
    
    return(xx)
    
  } else {
    stop()
  }
  
}

resolve_internal_bounds <- function(poly) {
  
  tryCatch({
    
    mls <- 
      poly %>%
      # st_transform(5070) %>% 
      st_cast("MULTILINESTRING") %>% 
      st_union() %>% 
      st_polygonize() %>% 
      st_sfc() %>% 
      st_sf() 
    
  }, error = function(e) 
    
  {  
    
    mls <-
      st_cast(
        # filter(
        st_transform(
          sf::st_as_sf(poly),
          5070
        ),
        # st_geometry_type(poly) == "GEOMETRYCOLLECTION"),
        "POLYGON",
        do_split = F
      ) %>%
      st_cast("MULTILINESTRING") %>% 
      st_union() %>% 
      st_polygonize() %>% 
      st_sfc() %>% 
      st_sf() 
    
    # mls
    
  }
  
  )
  # loop over each geometry in GEOMETRY COLLECTION and make a an SF geom,
  # rebind, and union into single POLYGON
  mls_lst <- lapply(1:length(mls$geometry[[1]]), function(i) {
    
    mls$geometry[[1]][i] %>% 
      st_sfc() %>% 
      st_geometry() %>% 
      st_set_crs(st_crs(mls$geometry)) %>% 
      st_as_sf() %>% 
      mutate(tmp_id = i)
    
  }) %>% 
    bind_rows()
  
  mls_lst <-
    mls_lst %>%
    sf::st_snap(mls_lst, tolerance = 1) %>%
    sf::st_union() %>% 
    sf::st_as_sf() %>% 
    dplyr::rename(geometry = x)
  
  mls_lst <-
    mls_lst %>% 
    dplyr::mutate(
      areasqkm = add_areasqkm(mls_lst)
    ) %>% 
    dplyr::bind_cols(
      dplyr::select(
        sf::st_drop_geometry(
          dplyr::slice_max(poly, areasqkm)
        ),
        -areasqkm
      )
    ) %>% 
    dplyr::relocate(names(poly))
  
  if (sf::st_geometry_type(mls_lst) == "MULTIPOLYGON") {
    
    warning(paste0("Failed to resolve internal boundaries, return is still a MULTIPOLYGON."))
    
    return(mls_lst)
    
  } else {
    
    message(paste0("Internal boundaries resolved, returning POLYGON."))
    
    return(mls_lst)
    
  }
  
  
}

try_resolve <- function(polys) {
  
  mls <- 
    polys %>%
    st_as_sf() %>% 
    st_transform(5070) %>% 
    st_collection_extract("POLYGON") %>%
    st_cast("MULTILINESTRING") %>% 
    st_union() %>% 
    st_cast("POLYGON") %>%
    st_sfc() %>%
    st_sf()
  
  return(mls)
}


resolve_geoc <- function(polys) {
  
  message(paste0("Resolving ", nrow(polys)," GEOMETRYCOLLECTIONS"))
  
  
  res <- lapply(1:nrow(polys), function(z) {
    
    if (sf::st_geometry_type(polys[z, ]) != "POLYGON") {
      
      # message(paste0(z, "/", nrow(polys)))
      
      mls <- try_resolve(polys[z, ])
      
      mls <- 
        mls %>%
        sf::st_snap(mls, tolerance = 5) %>%
        sf::st_union() %>%
        sf::st_as_sf() %>%
        dplyr::rename(geometry = x)
      
      mls <- 
        mls %>% 
        dplyr::mutate(
          areasqkm = add_areasqkm(mls)
        ) %>% 
        dplyr::bind_cols(
          dplyr::select(
            sf::st_drop_geometry(polys[z, ]),
            -areasqkm
          )
        ) %>% 
        dplyr::relocate(names(polys[z, ])) 
      
      mls
      
    } else {
      
      # message(paste0("SKIPPING ", z, " already a POLYGON"))
      
      sf::st_transform(polys[z, ], 5070)
      
    }
    
  }) %>% 
    dplyr::bind_rows()
  
  return(res)
  
}

resolve_mp2 <- function(polys) {
  
  res <- lapply(1:nrow(polys), function(z) {
    
    message(paste0(z, "/", nrow(polys)))
    
    # cast MULTIPOLYGONS as MULTILINESTRING 
    mls <- 
      polys[z, ] %>%
      sf::st_as_sf() %>% 
      sf::st_transform(5070) %>% 
      sf::st_collection_extract("POLYGON") %>%
      sf::st_cast("POLYGON", do_split = T) %>%
      sf::st_cast("MULTILINESTRING") %>% 
      sf::st_union() %>% 
      # sf::st_cast("POLYGON") %>%
      sf::st_polygonize() %>%
      sf::st_sfc() %>%
      sf::st_sf()
    
    mls_lst <- lapply(1:length(mls$geometry[[1]]), function(k) {
      
      mls$geometry[[1]][k] %>%
        st_sfc() %>%
        st_geometry() %>%
        st_set_crs(st_crs(mls$geometry)) %>%
        st_as_sf() %>%
        mutate(tmp_id = k)
      
    }) %>%
      dplyr::bind_rows() 
    # sf::st_union() %>%
    # sf::st_as_sf() %>%
    # dplyr::rename(geometry = x)
    # 
    mls_lst <-
      mls_lst %>%
      sf::st_snap(mls_lst, tolerance = 5) %>%
      sf::st_union() %>%
      sf::st_as_sf() %>%
      dplyr::rename(geometry = x)
    
    mls_lst <-
      mls_lst %>% 
      dplyr::mutate(
        areasqkm = add_areasqkm(mls_lst)
      ) %>% 
      dplyr::bind_cols(
        dplyr::select(
          sf::st_drop_geometry(polys[z, ]),
          -areasqkm
        )
      ) %>% 
      dplyr::relocate(names(polys[z, ]))
    
    mls_lst
    
    
  }) %>% 
    dplyr::bind_rows()
  
  return(res)
}

resolve_mp <- function(polys) {
  
  message(paste0("Resolving ", nrow(polys)," MULTIPOLYGONS"))
  # polys <- geoc
  
  res <- lapply(1:nrow(polys), function(z) {
    if (sf::st_geometry_type(polys[z, ]) != "POLYGON") {
      
      # cast MULTIPOLYGONS as MULTILINESTRING 
      mls <- 
        polys[z, ] %>%
        sf::st_as_sf() %>% 
        sf::st_transform(5070) %>% 
        sf::st_collection_extract("POLYGON") %>%
        sf::st_cast("POLYGON", do_split = T) %>%
        sf::st_cast("MULTILINESTRING") %>% 
        sf::st_union() %>% 
        # sf::st_cast("POLYGON") %>%
        sf::st_polygonize() %>%
        sf::st_sfc() %>%
        sf::st_sf()
      
      # loop over each geometry in GEOMETRY COLLECTION and make a an SF geom,
      # rebind, and union into single POLYGON
      mls_lst <- lapply(1:length(mls$geometry[[1]]), function(k) {
        
        mls$geometry[[1]][k] %>%
          st_sfc() %>%
          st_geometry() %>%
          st_set_crs(st_crs(mls$geometry)) %>%
          st_as_sf() %>%
          mutate(tmp_id = k)
        
      }) %>%
        dplyr::bind_rows() 
      
      # 
      mls_lst <-
        mls_lst %>%
        sf::st_snap(mls_lst, tolerance = 1) %>%
        sf::st_union() %>%
        sf::st_as_sf() %>%
        dplyr::rename(geometry = x)
      
      mls_lst <-
        mls_lst %>% 
        dplyr::mutate(
          areasqkm = add_areasqkm(mls_lst)
        ) %>% 
        dplyr::bind_cols(
          dplyr::select(
            sf::st_drop_geometry(polys[z, ]),
            -areasqkm
          )
        ) %>% 
        dplyr::relocate(names(polys[z, ]))
      
      if (sf::st_geometry_type(mls_lst) == "MULTIPOLYGON") {
        warning(paste0("Failed to resolve internal boundaries, return is still a MULTIPOLYGON."))
      } else {
        message(paste0("Internal boundaries resolved, returning POLYGON."))
      }
      
      mls_lst
      
    } else  { 
      # message(paste0("SKIPPING ", z, " already a POLYGON"))
      
      sf::st_transform(polys[z, ], 5070)
      
    }
    
  }) %>% 
    dplyr::bind_rows()
  
  return(res)
  
}
