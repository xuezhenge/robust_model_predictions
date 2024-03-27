rm()
rm(list=ls())
library(rgbif)
library(tictoc)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(BiodiversityR)
library(spThin)
library(raster)
library(biomod2)
library(usdm)
library(usmap)
library(terra)
# library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale='medium', returnclass= 'sf')

setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM')
#setwd("/home/xge/scratch/correlativeSDM")
species_list <- read.csv(paste0('species_names/','species_list.csv'),header = TRUE)
species_list <- as.list(species_list)$spec_list

log_con <- file(paste0("log_files/","data_cleanning_all_species.txt"), open="a")
len <- length(species_list)

all_records_thin <- {}
error_species <- {}
i <- 0
for (i in 0:(len-1)){
  tryCatch({
    time0 = Sys.time()
    i <- i + 1
    skip_to_next <- FALSE
    species <- species_list[i]
    cat(paste0("===============Running ",i," species: ",species,"==============="), file = log_con,sep="\n")
  #-----Step 1: obtain and clean occurrence data------------------------------
  occs <- read.csv(file.path('outputs',species,'data_cleaning/raw_data_global.csv'))
  no.raw_records <- nrow(occs)
  if(!file.exists(file.path('outputs',species,'data_cleaning/data_cleaned_coord_global.csv'))){
    #--------Checking species' coordinates-----------
    geo.clean <- clean_coordinates(x = occs, 
                                   lon = "decimalLongitude",
                                   lat = "decimalLatitude",
                                   species = "species", 
                                   countries = "countryCode",
                                   value = "clean",
                                   tests = c("capitals","centroids","equal","gbif","institutions","seas","zeros"))
    # export the global data after coordinate check
    occs.out<- geo.clean %>% dplyr::select(scientificName,species, decimalLongitude,decimalLatitude, year) 
    # remove spatial duplicates
    occs.out<- occs.out[!duplicated(occs.out),]
    no.cleaned_global_records <- nrow((occs.out))
    # export the global data after coordinate check
    write.csv(occs.out, 
              file.path('outputs',species,'data_cleaning/data_cleaned_coord_global.csv'),row.names = FALSE)}else{
    cat(paste0("---file exists!----"), file = log_con,sep="\n")   
              }
  no.raw_records <- nrow(occs)
  time1 = Sys.time()
  #-----Step 2: data thinning------------------------------
  occs <- read.csv(file.path('outputs',species,'data_cleaning/data_cleaned_coord_global.csv'))
  no.cleaned_records <- nrow(occs)
  colnames(occs) <- c("scientificname","spec","lon",'lat','year')
  ## Run spatial thinning
  output.dir = file.path('outputs',species,'data_thinning')
  if(!exists(output.dir)) dir.create(output.dir,recursive = TRUE)
  occs <- data.frame(occs$lon,occs$lat)
  colnames(occs) <- c('lon','lat')
  if(!file.exists(file.path('outputs',species,'data_thinning/data_thinned.csv'))){
  #remove duplicates
  occs <- distinct(occs)
  no.records <- nrow(occs)
  if (no.records > 100000){
    thinned_data <- ensemble.spatialThin.quant(occs, thin.km = 10, 
                                               runs = 1, silent = FALSE, verbose = FALSE, 
                                               LON.length = 21, LAT.length = 41)
    thinned_data <- ensemble.spatialThin(thinned_data, thin.km = 10, runs = 1)
  }else if(no.records > 5000 & no.records <=100000){
    thinned_data <- ensemble.spatialThin.quant(occs, thin.km = 10, 
                                      runs = 1, silent = FALSE, verbose = FALSE, 
                                      LON.length = 21, LAT.length = 21)
    thinned_data <- ensemble.spatialThin(thinned_data, thin.km = 10, runs = 1)
  }else if(no.records > 1000 & no.records <=5000){
    thinned_data <- ensemble.spatialThin.quant(occs, thin.km = 10, 
                                      runs = 1, silent = FALSE, verbose = FALSE, 
                                      LON.length = 6, LAT.length = 6)
    thinned_data <- ensemble.spatialThin(thinned_data, thin.km = 10, runs = 1)
  }else{
    thinned_data <- ensemble.spatialThin(occs, thin.km = 10, runs = 1)
  }
  write.csv(thinned_data,
            file.path('outputs',species,'data_thinning/data_thinned.csv'),
            row.names = FALSE)
  }else{
    cat(paste0("---file exists!----"), file = log_con,sep="\n")   
  }
  thinned_data <- read.csv(file.path('outputs',species,'data_thinning/data_thinned.csv'))
  no.3records <- c(species, no.raw_records, no.cleaned_records,nrow(thinned_data))
  all_records_thin <- rbind(all_records_thin, no.3records)
  cat(paste0("---------------have thinned species occurrence data-------------"), file = log_con,sep="\n")
  cat(no.3records, file = log_con,sep="\n")
  time2 = Sys.time()
  #-----Step 3: plot thinned occurrence records---------------
  ## save image of input data summary
  summary.dir <- file.path("outputs",species,'summary_doc')
  if(!exists(summary.dir)) dir.create(summary.dir,recursive = TRUE)
  occs <- read.csv(file.path('outputs',species,'data_thinning/data_thinned.csv'))
  no.thin_records <- nrow(occs)
  theme_set(theme_bw())
  p<-ggplot(data=world) +
    geom_sf()+
    geom_point(
      data = occs,
      aes(lon,lat,color = "red", alpha = 0.25)
    ) +
    theme(legend.position = "right") + 
    theme(panel.background = element_rect(colour = "black")) + 
    labs(title = paste0("Global occurrence records:",species, "(", no.thin_records," records)")) +
    theme(legend.position = "right")
  ggsave(paste0(summary.dir,'/',species,"_thinned_occ_map.png"),plot=p)
  time3 = Sys.time()
  time.clean <- time1 - time0
  time.thin <- time2 - time1
  total.time <- time3 - time0
  cat(paste0("**************Time of running ",i," species: ",species,"********"),  file = log_con,sep="\n")
  cat(paste0("Time of cleaning occurrence data: ", time.clean), file = log_con,sep="\n")
  cat(paste0("Time of thinning occurrence data: ",time.thin), file = log_con,sep="\n")
  cat(paste0("Total time of processing occurrence data: ",total.time), file = log_con,sep="\n")
  cat(paste0("****************************************************************"), file = log_con,sep="\n")
  }, 
  error=function(e){
    error_species <- rbind(error_species,species)
    skip_to_next <<- TRUE
    cat(paste0("Oops! --> Error in Loop ",i), file = log_con,sep="\n")
    cat("ERROR :",conditionMessage(e), "\n")})
  if(skip_to_next) {stop}
}

all_records_thin <- as.data.frame(all_records_thin)
names(all_records_thin) <- c('species', 'no.raw_records','no.cleaned_records', 'no.thin_records')

if(!exists("summary_tables")) dir.create("summary_tables")
write.csv(all_records_thin,paste0("summary_tables/all_records_thin.csv"))
cat(paste0("====================Summary CSV saved!!! Done!!! ================="), file = log_con,sep="\n")
cat(error_species, file = log_con,sep="\n")
close(log_con)

