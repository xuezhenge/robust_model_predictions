rm()
rm(list=ls())
library(rgbif)
#library(Taxonstand)
library(CoordinateCleaner)
# library(maps)
library(dplyr)
# library(devtools)
library(usmap)
library(ggplot2)
library(spThin)
library(raster)
library(biomod2)
library(usdm)
library(terra)
library(argparse)
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--speciesid", default='5', type="integer",
                    help="input the species id, e.g., 1")
args <- parser$parse_args()
species.id <- args$speciesid

#setwd("/home/xge/scratch/correlativeSDM")
setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM')
species_list <- read.csv(paste0('species_names/','species_list.csv'),header = TRUE)
species_list <- as.list(species_list)$spec_list
species <- species_list[species.id]

log_con <- file(paste0("log_files/","generatePA_",species,".txt"), open="a")
len <- length(species_list)

cat(paste0("===============Running ",species.id," species: ",species,"==============="), file = log_con,sep="\n")


## function to get PA dataset
get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA.table
  )
}

# function to get parameters for BIOMOD_FormatingData
get_PA_paras <- function(model.group,no.records){
  #Group A: 'MaxEnt','BIOCLIM'/'SRE','GLM','GAM','GBM',
  if (model.group == 'GroupA'){
    PAstrategy = 'random'
    PArep = 20
    if (no.records > 10000){
      PAno = no.records
    }else{PAno = 10000}
  }
  #Group B: 'MARS'
  else if(model.group == 'GroupB'){
    PAstrategy = 'random'
    PArep = 20
    PAno = no.records}
  #Group C: 'CTA','FDA', 'ANN', 'RF'
  else if(model.group == 'GroupC'){
    if(no.records<=200){
      PAstrategy = 'disk'
    }else{
      PAstrategy = 'sre'}
    PArep = 20
    PAno = no.records}
  paras = c(PAstrategy,PArep,PAno)
  return(paras)
}
# generate pseudo absence data for different models (which have 20 replicates)
get_pre_pa <- function(model.group,no.records,reps,myResp,myRespXY,bios.hist,sp){
  PA.paras <- get_PA_paras(model.group,no.records)
  PA.data <- get_PA_data(model.group,no.records,myResp,myRespXY,bios.hist,sp)
  csv.dir <- file.path('outputs',sp,'pre_pa_csv')
  if(!exists(csv.dir)) dir.create(csv.dir)
  output.dir <- file.path(csv.dir,model.group)
  if(!exists(output.dir)) dir.create(output.dir)
  
  # generate presence data
  pres.data <- myRespXY %>% 
    rename(x=lon,y=lat)
  pres.data$status <- 1
  # generate pseudo absence data (1 or 20 replicates) 
  for (i in 1:reps){
    pa.data.table <- get_PAtab(PA.data) %>% 
      filter(is.na(status))
    pa.data <- filter(pa.data.table,pa.data.table[i+3] == TRUE)
    pa.data <- pa.data[,c("x", "y","status")]
    pa.data$status <- 0
    # generate presence - pseudo absence data
    PA <- rbind(pres.data,pa.data)
    PA <- PA %>% rename(Species = status)
    # save presence - pseudo absence data
    file.dir <- paste0(output.dir,"/pre_pa_",model.group,i,".csv")
    write.csv(PA,file.dir, row.names = TRUE)
  }
}

# Generate pseudo-absence based on the PA parameters
get_PA_data <- function(model.group,no.records,myResp,myRespXY,bios.hist,sp){
  PA.paras <- get_PA_paras(model.group,no.records)
  if (PA.paras[1] == 'disk'){
    PA.data <- BIOMOD_FormatingData(
      resp.var = myResp,# presence data (all 1's)
      resp.xy = myRespXY,# coordinates of presences
      expl.var = bios.hist,# RasterStack
      resp.name = sp,# name of species
      PA.strategy = PA.paras[1], 
      PA.dist.min = 220000,
      PA.dist.max = NULL,
      PA.nb.rep = as.numeric(PA.paras[2]), 
      PA.nb.absences = as.numeric(PA.paras[3])
    )} else {
      PA.data <- BIOMOD_FormatingData(
        resp.var = myResp,# presence data (all 1's)
        resp.xy = myRespXY,# coordinates of presences
        expl.var = bios.hist,# RasterStack
        resp.name = sp,# name of species
        PA.strategy = PA.paras[1], 
        PA.nb.rep = as.numeric(PA.paras[2]), 
        PA.nb.absences = as.numeric(PA.paras[3])
      )}
  return(PA.data)
}

sel_ev_hist<- function(bios.hist){
  ex <- raster::extract(bios.hist,spg)
  head(ex)
  ex <- na.omit(ex)
  head(ex)
  #identify the variables which have collinearity problem
  v <- vifstep(ex, th=4)
  v
  #exclude these variables
  bios.hist <- exclude(bios.hist,v)
  return(bios.hist)
}

#-----Step 1: select climate variables--------------------------------------
unlink(file.path('outputs',species,'sp_bios_names'), recursive = TRUE)
# 1.load occurrence data
spg <- read.csv(file.path('outputs',species,'data_thinning/data_thinned.csv'))
spg$species.id <- 1
spg.cc <- spg
no.records <- nrow(spg)
coordinates(spg) <- c('lon','lat')
# 2.prepare environment data-------
#worldclim 2 data
hist.dir <- "inputs/world_clim2/world/1971-2000"
sp.bios.dir <- file.path('outputs',species,"sp_bios_names")
if(!exists(sp.bios.dir)) dir.create(sp.bios.dir)
fns<- list.files(hist.dir)
climate.dir <- file.path(hist.dir,fns)
bios.hist <- raster::stack(climate.dir)
# 3.select historical environment variables based on VIF (multi-collinearity)
bios.hist <- sel_ev_hist(bios.hist)
sel.evnames <- names(bios.hist)
#save the select climate variables names
write.csv(sel.evnames,file.path(sp.bios.dir,'hist_bio_names.csv'))
cat(paste0("------------------Have output the climate rasters---------------"), file = log_con,sep="\n")
#----Step 2: generate dataset that contains presence and pseudo-absence data----
#1. load climate variables
# import the selected climate variable names
hist.dir <- "inputs/world_clim2/world/1971-2000"
sel.evnames <- read.csv(file.path('outputs',species,'sp_bios_names/hist_bio_names.csv'))
sel.evnames <- as.vector(sel.evnames$x)
fns<- paste0(sel.evnames,".tif")
climate.dir <- paste0(hist.dir,'/',fns)
bios.hist <- raster::stack(climate.dir)
#2.format occurrence data----------
# Get corresponding presence/absence data
# load occurrence data
spg.cc <- read.csv(file.path('outputs',species,'data_thinning/data_thinned.csv'))
ex_env <- raster::extract(bios.hist,spg.cc)
spg.env <-cbind(spg.cc,ex_env)
spg.cc <- na.omit(spg.env)
spg.cc <- spg.cc[,c('lon','lat')]
spg.cc$species.id <- 1
no.records <- nrow(spg.cc)
myResp <- spg.cc['species.id']
# Get corresponding XY coordinates
myRespXY <- spg.cc[, c('lon', 'lat')]
no.records <- nrow(spg.cc)
#3. generate presence pseudo-absence data
#4. generate presence pseudo-absence data and write pre-pa data into csv format
# Get pseudo -absence for different models when generating pseudo-absence
#Group A: 'MaxEnt','BIOCLIM'/'SRE','GLM','GAM','GBM',
#GroupB: 'MARS'
#Group C: 'CTA','FDA','ANN', 'RF'
get_pre_pa('GroupA',no.records,reps=20,myResp,myRespXY,bios.hist,species)
get_pre_pa('GroupB',no.records,reps=20,myResp,myRespXY,bios.hist,species)
get_pre_pa('GroupC',no.records,reps=20,myResp,myRespXY,bios.hist,species)
cat(paste0("------------------Have generated PA data---------------"), file = log_con,sep="\n")
close(log_con)
