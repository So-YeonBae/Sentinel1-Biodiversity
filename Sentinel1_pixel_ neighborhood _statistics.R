# Written by Soyeon Bae and Benjamin Leutner 
# Bae et al. (2019) Radar vision in the mapping of forest biodiversity from space (submitted to Nature Communications)
# Supplementary Information
# III. Supplementary Methods 
# S3.3 Derivation of pixel- and neighborhood-based summary statistics 

# Call the packages (if you didn't install them yet, install them first by using the following line)
# install.packages(c("rgeos","raster","rgdal","glcm")) # If you already installed them, please skip this line.
library(rgeos)
library(raster)
library(rgdal)
library(glcm)

## Define your input directory containing example input files
inputDirectory <- "D:/Dropbox/TREESCAPE/Writing/Sentinel_lidar/5_Revision_1st/Sentinel1_ExampleData2/"

## Common settings
outputDirectory <- file.path(inputDirectory, "expected_output/")
dir.create(outputDirectory)
scaleFactor  <- 1e4      ## Scale factor used for compressing rasters
dataType     <- "INT4S"  ## Data type used to store raster data (integers between -2,147,483,647 and 2,147,483,647)
samplingSite <- "ALB"    ## Name of study site (forest) used in this example
cOpts  <- c("COMPRESS=DEFLATE", "INTERLEAVE=BAND") ## Compression for GeoTiff files and efficient file structure


## Required input files for process 1 (Calculate the seasonal median) ##
# Data description:
# There are two multi-layer GeoTiff per polarisation (VV, VH), containing the gamma_0
# backscatter coefficients derived via ESA SNAP toolbox. Each stack contains 107 layers 
# which represent all acquisitions between 2016-01-01 and 2016-12-31.
# 1) Stack for VV polarisation:
paste0(inputDirectory, "ALB_gamma0_VV.tif")
# 2) Stack for VH polarisation:
paste0(inputDirectory, "ALB_gamma0_VH.tif")
# 3) Data acquisition dates:
paste0(inputDirectory, "ALB_acquisition_dates.csv")

## Required input files for process 2 (Calculate GLCM) ##
# 1) Calculated mean values of all our 1-ha plots
paste0(inputDirectory, "S1_mean.csv")
# 2) output files from process 1

## Required input files for process 3 (Create 9x9 neighborhood Map) ##
# 1) output files from process 1

## Required input files for process 4 (Extract the metrics from the example plots) ##
# 1) Shapefile of center of example plots
paste0(inputDirectory,"ALB_example_plots.shp")
# 2) output files from process 2 and 3


############# Process 1. Calculate Yearly, Winter, Summer Median and difference and ratio between VV and VH ############

## Acquisition dates of example data
dates <- read.csv(paste0(inputDirectory, "ALB_acquisition_dates.csv"), header = F, col.names = "Dates")
dates <- as.character(dates$Dates)
months <- substr(dates, 5, 6)
winterScenes <- which(months %in% c("12", "01", "02"))
summerScenes <- which(months %in% c("07", "08", "09"))
outBasename  <- paste0(outputDirectory, "/L2_", samplingSite, "_")

for (pol in c("VH","VV")){
	
	# Import a Geotiff file processed in SNAP
	s1stack <- stack(paste0(inputDirectory, "ALB_gamma0_",pol,".tif"))
	
	# Stack summer and winter scenes
	summer <- s1stack[[summerScenes]]
	winter <- s1stack[[winterScenes]]
	
	# Calculate Yearly, Winter, Summer Median
	
	med_summer <- calc(summer, fun = median, na.rm = TRUE, 
			filename = paste0(outBasename, pol, "_summer.tif"),
			datatype = dataType, options = cOpts, overwrite = TRUE)
	med_winter <- calc(winter, fun = median, na.rm = TRUE,
			filename = paste0(outBasename, pol, "_winter.tif"), 
			datatype = dataType, options = cOpts, overwrite = TRUE)
	median_all <- calc(s1stack, fun = median, na.rm = TRUE , 
			filename = paste0(outBasename, pol, "_year.tif"),  
			datatype = dataType, options = cOpts, overwrite = TRUE)
	
	# Calculate difference between winter and summer median
	diff_season <- calc(stack(med_summer, med_winter), fun = function(x) {x[,1]-x[,2]}, forcefun = TRUE,
			filename = paste0(outBasename, pol, "_diff_summer_winter.tif"),
			datatype = dataType, options = cOpts, overwrite = TRUE)
	
}

## Import yearly meadian backscatter values
vv_med <- raster(paste0(outBasename, "VV_year.tif"))
vH_med <- raster(paste0(outBasename, "VH_year.tif"))

## Calculate difference between median polarisations
overlay(vv_med, vH_med, fun = function(vv,vh) {(vv - vh)}, 
		forcefun = TRUE, filename = paste0(outBasename,"diff_VH_VV.tif"),
		datatype = dataType, options = cOpts, overwrite = TRUE)

## Calculate polarisation ratios
overlay(vv_med, vH_med,fun =  function(vv,vh) {vv/vh*10000}, # Scale the ratio values due to data type of integers
		forcefun = TRUE, filename = paste0(outBasename,"ratio_VV_VH.tif"),
		datatype = dataType, options = cOpts, overwrite = TRUE)


############# Process 2. Calculate 9x9 GLCM #############################

## Retrieve the 5% and 95% percentiles from all plots (as opposed to min/max, which are too sensible for extrema inherent in SAR data)
s1_all    <- read.csv(paste0(inputDirectory, "S1_mean.csv")); names(s1_all)
s1_all$S1_ratio_VV_VH <- s1_all$S1_ratio_VV_VH*10000 # Scale the ratio values due to data type of integers
s1_mean <- s1_all[,-c(1:2)]; names(s1_mean)

min5  <- sapply(s1_mean, quantile, probs = 0.01, na.rm = TRUE)
max95 <- sapply(s1_mean, quantile, probs = 0.99, na.rm = TRUE)

minmax_tab <- data.frame(pol = gsub("S1_","", names(s1_mean)), totMin = min5, totMax = max95)
minmax_tab$layers <- list.files(outputDirectory, pattern = "^L2.*tif$", full.names = TRUE)

## Calculate GLCM Texture 
glcm_list <- lapply(1:nrow(minmax_tab), function(k) {
			img_k   <- raster(minmax_tab$layers[k])
			img_k_c <- clamp(img_k, lower = minmax_tab$totMin[k], 
					         upper = minmax_tab$totMax[k], useValues = TRUE)
			
			texture_glcm <- glcm(img_k_c,
					n_grey = 32, 
					window = c(9,9),
					min_x = minmax_tab$totMin[k], 
					max_x = minmax_tab$totMax[k], 
					statistics = c("dissimilarity", "entropy"),
					scale_factor = scaleFactor,
					asinteger = TRUE)
		})
texture_stack <- stack(glcm_list)
names(texture_stack) <-  paste0(rep(minmax_tab$pol, each = 2), c("_GLCM_9x9_DIS", "_GLCM_9x9_ENT"))
outfile    <- paste0(outputDirectory, "/L3_texture_",samplingSite,"_N9x9_glcm.tif")
writeRaster(texture_stack, filename = outfile, datatype = dataType, options = cOpts, overwrite = TRUE)


############# Process 3. Create 9x9 neighborhood Map ##################################

## Calculate mean in 9x9 neighborhood
neighborhood_list <- lapply(1:nrow(minmax_tab), function(k) {
			img_k <- raster(minmax_tab$layers[k])
			focal_mean <- focal(img_k, w = matrix(nrow=9,ncol=9, data=1), fun = mean, na.rm = TRUE )
		})
neighborhood_stack <- stack(neighborhood_list)
names(neighborhood_stack) <-  paste0(minmax_tab$pol, "_9x9_mean")
outfile    <- paste0(outputDirectory, "/L3_neighborhood_",samplingSite,"_N9x9_mean.tif")
writeRaster(neighborhood_stack, filename = outfile, datatype = dataType, options = cOpts, overwrite = TRUE)

## Calculate standard deviation in 9x9 neighborhood
heterogeneity_list <- lapply(1:nrow(minmax_tab), function(k) {
			img_k <- raster(minmax_tab$layers[k])
			focal_sd <- focal(img_k, w = matrix(nrow=9,ncol=9, data=1), fun = sd, na.rm = TRUE )
		})
heterogeneity_stack <- stack(heterogeneity_list)
names(heterogeneity_stack) <-  paste0(minmax_tab$pol, "_9x9_SD")
outfile    <- paste0(outputDirectory, "/L3_heterogeneity_",samplingSite,"_N9x9_SD.tif")
writeRaster(heterogeneity_stack, filename = outfile, datatype = dataType, options = cOpts, overwrite = TRUE)


############# Process 4. Extract the metrics from the example plots ##################################

# Import all variables as a single raster stack
grids <- list.files(outputDirectory, pattern = "^L3.*tif$", full.names = TRUE)
s     <- stack(grids)

# Create the variable names
grids_name <- c(unlist(sapply(paste0(lapply(strsplit(basename(grids),"_"),function(x){x[2]}),"_stack"), function(x){names(get(x))})))

# Read the polygon shapefile of study plots
plots <- readOGR(dsn = paste0(inputDirectory,samplingSite,"_example_plots.shp"),layer = paste0(samplingSite,"_example_plots"))

# Extract raster mean and standard deviation of all pixels falling within each polygon area (plots)

ext <- lapply(as.list(s), extract, y = plots)
ex  <- matrix(data = unlist(ext), nrow = length(plots))
colnames(ex) <- grids_name

# Write to a data frame
df <- data.frame(Explrtr = plots$EP,ex)
df$ratio_VV_VH_9x9_SD <- df$ratio_VV_VH_9x9_SD/10000 # Rescale the ratio values
df$ratio_VV_VH_9x9_mean <- df$ratio_VV_VH_9x9_mean/10000 # Rescale the ratio values

#write to a CSV file
write.table(df, file = paste0(outputDirectory,"/S1_",samplingSite,"_plots_stat.csv"), sep=",", row.names = F, col.names = T)



