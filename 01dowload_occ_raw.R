
rm(list=ls())
library(rgbif)
library(dplyr)
library(argparse)
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--speciesid", default='6', type="integer",
                    help="input the species id, e.g., 1")
args <- parser$parse_args()
species.id <- args$speciesid

setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM')
#setwd("/home/xge/scratch/correlativeSDM")
species_list <- read.csv(paste0('species_names/','species_list.csv'),header = TRUE)
species_list <- as.list(species_list)$spec_list
species <- species_list[species.id]
get_occ_raw <- function(species){
  species0 <- gsub("_", " ", species)
  occs0 <- occ_search(scientificName = species0, hasCoordinate = TRUE,limit = 1000000)
  occs <- occs0$data
  no.raw_records <- nrow(occs) #number of records with coordinates
  colnames(occs) #Column names returned from gbif follow the DarwinCore standard
  #exporting raw data
  clean.dir <- file.path('outputs',species,'data_cleaning')
  if(!exists(clean.dir)) dir.create(clean.dir, recursive = TRUE)
  write.csv(occs, 
            file.path('outputs',species,'data_cleaning/raw_data_global.csv'), 
            row.names = FALSE)
}
if(!exists(paste0('log_files/'))) dir.create(paste0('log_files/'))
log_con <- file(paste0('log_files/',species,'.txt'), open="a")
get_occ_raw(species)
cat(paste0("===============Have downloaded ",species.id," species: ",species,"==============="), file = log_con,sep="\n")

