# Script Summary:
# downloading and processing spatial data related to the National Hydrography Dataset Plus (NHDPlus) from the Environmental Protection Agency (EPA).
# 
# script load a config.R file that sets up file struture and obtains a data frame listing objects from the EPA bucket. Then, it defines and downloads various files related to catchments, flowlines, and water bodies, from which shapefiles are extracted, and data is written to geopackages using the ogr2ogr tool. Finally, the code downloads and extracts burn component files, and the resultant data is written to geopackages

# load config.R file to set up file structure
source("workflow/nhdplusv2/config.R")

# get EPA bucket from S#
epa = get_bucket_df(epa_bucket, max = Inf)

# New Attributes  ---------------------------------------------------------

# parquet file name
n = 'enhd_nhdplusatts.parquet'

# download parquet file and save to base_dir
item_file_download(
  sb_id        = '63cb311ed34e06fef14f40a3',
  names        = n,
  destinations = glue("{base_dir}/{n}")
)

# Catchments --------------------------------------------------------------

# Define
cats = grep("_NHDPlusCatchment_", epa$Key, value = TRUE)

# make dataframe with paths to S3 Reousrces for catchments
all_cats = data.frame(
  key    = cats,
  VPU    = sapply(strsplit(cats, "_"), `[`, 3),
  region = sapply(strsplit(cats, "_"), `[`, 2),
  link   = glue('https://{epa_bucket}.s3.amazonaws.com/{cats}')) %>%
  mutate(
    outfile = glue('{epa_download}/NHDPlus{VPU}.7z'),
    shp     = glue("{epa_download}NHDPlus{region}/NHDPlus{VPU}/NHDPlusCatchment/Catchment.shp")
  ) %>%
  filter(!shp %in% list.files(epa_download, recursive  = TRUE, pattern = ".shp$"))

# Download & Extract
# Loop through catchments dataframe and get objects from S3
for (i in 1:nrow(all_cats)) {
  save_object(object = all_cats$key[i], bucket = epa_bucket, file = all_cats$outfile[i])
  
  archive_extract(all_cats$outfile[i], dir = epa_download)
  
  unlink(all_cats$outfile[i])
  log_info("Downloaded and extracted ", all_cats$VPU[i])
}

#### Save data to geopackages

# add paths to geopackage data
all_cats =
  all_cats %>%
  mutate(gpkg = glue('{catchments_dir}NHDPlus{all_cats$VPU}.gpkg')) %>%
  filter(!gpkg %in% list.files(catchments_dir, full.names = TRUE))

# if the number of rows in all_cats are >0, then make new ogr2ogr geopackage calls and execute
if(nrow(all_cats) > 0){
  
  # adding "-nlt MULTIPOLYGON" option as the warning suggests, this did make the warning go away...
  calls = paste('ogr2ogr -f GPKG -nlt MULTIPOLYGON', all_cats$gpkg, all_cats$shp)
  
  # original ogr2ogr command that produces warning (below)
  # calls = paste('ogr2ogr -f GPKG', all_cats$gpkg, all_cats$shp)
  
  for(i in 1:length(calls)){
    system(calls[i])
    log_info("Produced ", basename(all_cats$gpkg[i]))
  }
}

### Warning I get when I don't have the "-nlt MULTIPOLYGON" option set as the warning suggests...
# Warning 1: A geometry of type MULTIPOLYGON is inserted into layer Catchment of geometry type POLYGON, which is not normally allowed by the GeoPackage specification, but the driver will however do it. To create a conformant GeoPackage, if using ogr2ogr, the -nlt option can be used to override the layer geometry type. This warning will no longer be emitted for this combination of layer and feature geometry type.

# Flowlines + Waterbodies --------------------------------------------------------------

snap = grep("_NHDSnapshot_", epa$Key, value =TRUE)

# NHDSnapshot dataframe of file paths to get and save data from/to
all_snap = data.frame(
  key    = snap,
  VPU    = sapply(strsplit(snap, "_"), `[`, 3),
  region = sapply(strsplit(snap, "_"), `[`, 2),
  link   = glue('https://{epa_bucket}.s3.amazonaws.com/{snap}')
) %>%
  mutate(
    outfile = glue("{epa_download}NHDPlusSnapshot{VPU}.7z"),
    fl_shp  = glue("{epa_download}NHDPlus{region}/NHDPlus{VPU}/NHDSnapshot/Hydrography/NHDFlowline.shp"),
    wb_shp  = glue("{epa_download}NHDPlus{region}/NHDPlus{VPU}/NHDSnapshot/Hydrography/NHDWaterbody.shp")
  ) %>%
  filter(!fl_shp %in% list.files(epa_download, recursive  = TRUE, pattern = ".shp$"))

# Download Flowlines ----------------
# loop through Snap dataframe
#  --> get/save objects from S3 EPA Bucket
#  --> archive data to outfile directory
for(i in 1:nrow(all_snap)){
  
  save_object(
    object = all_snap$key[i],
    bucket = epa_bucket, 
    file   = all_snap$outfile[i]
  )
  
  archive_extract(all_snap$outfile[i], dir = epa_download)
  
  unlink(all_snap$outfile[i])
  log_info("Downloaded and extracted ", all_snap$VPU[i])
}

# Make geopackages -----------------

# add paths flowline and water body geopackage 
all_snap = 
  all_snap %>%
  mutate(
    fl_gpkg = glue("{fl_dir}NHDPlus{VPU}.gpkg"),
    wb_gpkg = glue("{wb_dir}NHDPlus{VPU}.gpkg")
  ) %>%
  filter(!wb_gpkg %in% list.files(wb_dir, full.names = TRUE))

# if the number of rows in all_snap are >0, then make new ogr2ogr geopackage calls and execute
if(nrow(all_snap) > 0){
  calls = c(
    paste('ogr2ogr -f GPKG -nlt MULTIPOLYGON', all_snap$wb_gpkg, all_snap$wb_shp),
    paste('ogr2ogr -f GPKG -nlt MULTILINESTRING', all_snap$fl_gpkg, all_snap$fl_shp)
  )
  for(i in 1:length(calls)){
    system(calls[i])
    log_info("Downloaded and extracted ", i, "/", length(calls))
  }
}

# BLE --------------------------------------------------------------

# get paths to NHD Burn components
burn = grep("_NHDPlusBurnComponents_", epa$Key, value =TRUE)

# make dataframe for burn paths 
all_burn = data.frame(
  key    = burn,
  VPU    = sapply(strsplit(burn, "_"), `[`, 3),
  region = sapply(strsplit(burn, "_"), `[`, 2),
  link   = glue('https://{epa_bucket}.s3.amazonaws.com/{burn}')
) |>
  mutate(
    outfile  = glue("{epa_download}NHDPlusBurnComponents{VPU}.7z"),
    ble_shp  = glue("{epa_download}NHDPlus{region}/NHDPlus{VPU}/NHDPlusBurnComponents/BurnAddLine.shp")
  ) %>%
  filter(!ble_shp %in% list.files(epa_download, recursive  = TRUE, pattern = ".shp$"))

# Download BLE ----------------
# loop through burn dataframe
#  --> get/save objects from S3 EPA Bucket
#  --> archive data to outfile directory
for(i in 1:nrow(all_burn)){
  
  save_object(
    object = all_burn$key[i], 
    bucket = epa_bucket, 
    file   = all_burn$outfile[i]
  )
  
  archive_extract(all_burn$outfile[i], dir = epa_download)
  
  unlink(all_burn$outfile[i])
  log_info("Downloaded and extracted ", all_burn$VPU[i])
}

# Make BLE geopackages -----------------

# add paths to burn geopackages
all_burn = 
  all_burn %>%
  mutate(
    ble_gpkg = paste0(ble_dir, "NHDPlus", VPU, ".gpkg")
  ) %>%
  filter(!ble_gpkg %in% list.files(ble_dir, full.names = TRUE))

# if the number of rows in all_burn are >0, then make new ogr2ogr geopackage calls and execute
if(nrow(all_burn) > 0){
  calls = paste('ogr2ogr -f GPKG -nlt MULTILINESTRING', all_burn$ble_gpkg, all_burn$ble_shp)
  
  for(i in 1:length(calls)){
    system(calls[i])
    log_info("Downloaded and extracted ", i, "/", length(calls))
  }
}

