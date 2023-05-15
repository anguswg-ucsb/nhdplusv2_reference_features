# This R code seems to be part of a data processing pipeline. Here's what it does:

# 1. Loads a configuration file called "config.R" from the "workflow/nhdplusv2" directory. This file sets up the file structure for the rest of the pipeline.

# 2. Creates a list of file paths for files in the "fl_dir" and "ble_dir" directories.

# 3. Reads in a Parquet file containing new attributes.

# 4. Tries to get the VAA using the "nhdplusTools::get_vaa()" function. If this fails, it catches the error, unlinks the file specified by "nhdplusTools::get_vaa_path()", and retries getting VAA data

# 5. Sets the "par" variable to 3.

# 6. Loops over each file path in the "fl_paths" list:
# - For each file path, it extracts the VPU name from the file path,
# - finds the corresponding file path in the "ble_paths" list, 
# - creates an output file path for the flowlines.

# 7. If the output file does NOT already exist:
# - reads in the flowlines and boundary line files
# - aligns their names
# - joins them with the "vaa" object
# - adds the new attributes 
# - drops the geometry from the BLE file
# - left joins it with the flowlines on the "COMID" column
# - replaces the geometry of certain flowlines based on the "StartFlag" and "Divergence" columns.

# 8. creates a custom network by: 
#  - selecting the "COMID", "FromNode", "ToNode", and "Divergence" columns
#  - getting the tocomid for each comid using the "nhdplusTools::get_tocomid()" function
#  - selecting the "comid" and "override_tocomid" columns.

# 9. left joins the flowlines with the custom network on the "COMID" column and replaces any "override_tocomid" values of 0 with the corresponding "tocomid" values.

# 10. Check if each flowline is a headwater and for sure flows to something. If not, it stops with an error message.

# 11. Match each "override_tocomid" value with the corresponding "COMID" value in the flowlines and creates a list of flowlines and network data frames.

# 12. Creates a cluster of parallel workers using the "makeCluster()"with the number of workers specified by the "par" variable.

# 13. Applies the "nhdplusTools::fix_flowdir()" function to each pair of flowlines and networks using the "pblapply()" function, to run the function in parallel

# 14. Stop the cluster of parallel workers.

# 15. IDs any flowlines that caused errors in the "nhdplusTools::fix_flowdir()" function and sets their corresponding values in the "check" logical vector to FALSE.

# 16. Check if any errors other than empty geometry occurred during the "nhdplusTools::fix_flowdir()" function. If so, it stops with an error message.

# 17. Replace the geometry of all flowlines that were NOT IDs as errors with the fixed geometry produced by the "nhdplusTools::fix_flowdir()" function.

# 18. Select the flowlines with "COMID" values in the "new_atts$comid" vector, removes the "override_tocomid" column, and writes the resulting sf dataframe to the output file path.

# load config.R file to set up file structure
source("workflow/nhdplusv2/config.R")

fl_paths  = list.files(fl_dir,  full.names = TRUE)
ble_paths = list.files(ble_dir, full.names = TRUE)

new_atts = read_parquet(glue("{base_dir}/enhd_nhdplusatts.parquet"))

tryCatch({
  
  vaa <- nhdplusTools::get_vaa()
  
}, error = function(e) {
  
  message("Error caught: ", conditionMessage(e))
  
  # Attempt to unlink the file
  unlink(nhdplusTools::get_vaa_path())
  
  # Retry getting VAA
  vaa <- nhdplusTools::get_vaa()
})

par = 3

# Loop through all the flowline paths in fl_paths
for(i in 1:length(fl_paths)){
  
  # Get the current flowline path
  fl_path   = fl_paths[i]
  
  # Get the VPU identifier for the current flowline file
  which_VPU = gsub(".gpkg", "", gsub("NHDPlus", "", basename(fl_path)))
  
  # Get the BLE path for the current VPU
  ble_path  = ble_paths[grep(which_VPU, basename(ble_paths))]
  
  # Define the output file path for the processed flowlines
  outfile   = glue("{reference_dir}flowlines_{which_VPU}.gpkg")
  
  # Check if the output file already exists, and if not, process the flowlines
  if(!file.exists(outfile)){
    
    # Read the flowline file and perform some preprocessing steps
    nhd <- 
      fl_path %>% 
      read_sf() %>% 
      st_zm() %>%
      align_nhdplus_names() %>%
      left_join(
        vaa, 
        by = c("COMID" = "comid")
      ) %>%
      select(COMID, fromnode, tonode, startflag, streamcalc, divergence, dnminorhyd) %>%
      left_join(
        new_atts,
        by = c("COMID" = "comid")
      ) %>%
      align_nhdplus_names() %>%
      mutate(LENGTHKM  = add_lengthkm(.))
    
    # Read the BLE file and perform some preprocessing steps
    ble <- 
      ble_path %>% 
      read_sf() %>%
      rename(COMID = LineID)
    
    # Join the NHDPlus flowlines with the BLE lines where applicable
    ble <- 
      left_join(
        select(st_drop_geometry(nhd), COMID), 
        ble, 
        by = "COMID"
      ) %>%
      st_as_sf() %>%
      st_zm()
    
    # Define a flag to indicate where the BLE lines will be used
    flag = !st_is_empty(st_geometry(ble)) & (nhd$StartFlag == 1 | nhd$Divergence == 2)
    
    # Replace the flowline geometry with the BLE line geometry where applicable
    st_geometry(nhd)[flag] <- st_geometry(ble)[flag]
    
    # Generate a custom network for the flowlines based on the override_tocomid values
    custom_net <- 
      nhd %>% 
      st_drop_geometry() %>%
      select(COMID, FromNode, ToNode, Divergence) %>%
      nhdplusTools::get_tocomid(remove_coastal = FALSE) %>%
      select(comid, override_tocomid = tocomid)
    
    # Update the override_tocomid values to account for nodes with no downstream flow
    nhd <- 
      left_join(
        nhd, 
        custom_net, 
        by = c("COMID" = "comid")
      ) %>%
      mutate(
        override_tocomid = ifelse(toCOMID == 0, override_tocomid, toCOMID)
      )
    
    # Check for any invalid override_tocomid values, that each flowline is a headwater and leads somewhere
    check <- !nhd$COMID %in% nhd$override_tocomid &
      !(nhd$override_tocomid == 0 | is.na(nhd$override_tocomid) |
          !nhd$override_tocomid %in% nhd$COMID)
    
    # Filter the flowlines to only those with
    check_direction <- filter(nhd, check)
    
    if(!all(check_direction$override_tocomid[check_direction$override_tocomid != 0] %in% nhd$COMID)){
      stop("this won't work")
    }
    
    matcher <- match(check_direction$override_tocomid, nhd$COMID)
    
    matcher <- nhd[matcher, ]
    
    fn_list <- Map(list,
                   split(check_direction, seq_len(nrow(check_direction))),
                   split(matcher, seq_len(nrow(matcher))))
    # creates a cluster of parallel workers 
    cl <- makeCluster(par)
    
    # check and fix each pair of flowlines and networks using nhdplusTools::fix_flowdir()
    new_geom <- pblapply(
      cl  = cl,
      X   = fn_list,
      FUN = function(x) {
        nhdplusTools::fix_flowdir(
          comid   = x[[1]]$COMID,
          fn_list = list(
            flowline  = x[[1]],
            network   = x[[2]],
            check_end = "end"
          )
        )
      }
    )
    
    # stop workers
    stopCluster(cl)
    
    # ID any problems
    error_index <- sapply(new_geom, inherits, what = "try-error")
    errors      <- filter(nhd, COMID %in% check_direction$COMID[error_index])
    check[which(nhd$COMID %in% errors$COMID)] <- FALSE
    
    if(!all(sapply(st_geometry(errors), st_is_empty))) {
      stop("Errors other than empty geometry from fix flowdir")
    }
    
    st_geometry(nhd)[check] <- do.call(c, new_geom[!error_index])
    
    # select relavent COMIDs and save out
    nhd %>% 
      filter(COMID %in% new_atts$comid) %>%
      select(-override_tocomid) %>%
      write_sf(outfile, "flowlines")
  }
  
  logger::log_info("Finished VPU ", which_VPU, " flowlines")
}
