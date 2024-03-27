rm()
rm(list=ls())
library(raster)
library(sp)
library(usdm)
library(biomod2)
library(dplyr)
library(tidyverse)
library(mapview)
library(maps)
library(graphics)
library(blockCV)
library(sf)
library(tidyr)
library(argparse)
#https://cran.r-project.org/web/packages/blockCV/vignettes/BlockCV_for_SDM.html
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--speciesid", default='6', type="integer",
                    help="input the species id, e.g., 1")
args <- parser$parse_args()
species.id <- args$speciesid

setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM')
species_list <- read.csv(paste0('species_names/','species_list.csv'),header = TRUE)
species_list <- as.list(species_list)$spec_list
species <- species_list[species.id]

log_con <- file(paste0("log_files/generateDST",species,".txt"), open="a")
#get data split table (blockCV) and preparation for null model 
get_DST <- function(pa.id,model.group){
  #create folders
  DST.dir <- file.path('outputs',species,'DataSplitTable',model.group)
  if(!exists(DST.dir)) dir.create(DST.dir,recursive = TRUE)
  Null.train_pre.dir <- file.path('outputs',species,'nullmodel/train_pre',model.group)
  if(!exists(Null.train_pre.dir)) dir.create(Null.train_pre.dir,recursive = TRUE)
  Null.rest.dir <- file.path('outputs',species,'nullmodel/rest_data',model.group)
  if(!exists(Null.rest.dir)) dir.create(Null.rest.dir,recursive = TRUE)
  #--------Step1: Load pre-pa data---------
  pre_pa.dir <- file.path('outputs',species,"pre_pa_csv",model.group)
  #---------Step2: Load historical climate data -------
  hist.dir <- "inputs/world_clim2/world/1971-2000"
  sel.evnames <- read.csv(file.path('outputs',species,'sp_bios_names/hist_bio_names.csv'))
  sel.evnames <- as.vector(sel.evnames$x)
  fns<- paste0(sel.evnames,".tif")
  climate.dir <- file.path(hist.dir,fns)
  bios.hist <- raster::stack(climate.dir)
  #---------Step3: Load pesudo-absence and presence data -------
  pre_pa.dir <- file.path('outputs',species,"pre_pa_csv",model.group)
  file.list <- list.files(pre_pa.dir)
  file <- paste0('pre_pa_',model.group,pa.id,'.csv')
  pre_pa <- read.csv(file.path(pre_pa.dir,file))
  pre.records <- nrow(pre_pa[which(pre_pa$Species==1), ])
  #---------Step4: Splitting between training and test records (blockCV)-------
  pa_data <- st_as_sf(pre_pa, coords = c("x", "y"), crs = crs(bios.hist))
  # spatial blocking by rows with systematic assignment
  sb2 <- spatialBlock(speciesData = pa_data, # presence-background data
                      species = "Species",
                      rasterLayer = bios.hist,
                      rows = 3,
                      cols = 5,
                      k = 5,
                      selection = "systematic",
                      biomod2Format = TRUE)
  
  # spatial blocking by rows and columns with checkerboard assignment
  sb3 <- spatialBlock(speciesData = pa_data,
                      species = "Species",
                      rasterLayer = bios.hist,
                      rows = 10,
                      cols = 15,
                      selection = "checkerboard",
                      biomod2Format = TRUE)
  
  # environmental clustering
  eb <- envBlock(rasterLayer = bios.hist,
                 speciesData = pa_data,
                 species = "Species",
                 k = 5,
                 standardization = "standard", # rescale variables between 0 and 1
                 rasterBlock = FALSE,
                 numLimit = 1)
  
  # 2. Defining the folds for DataSplitTable
  # note that biomodTable should be used here not folds
  # use generated folds from spatialBlock in previous section
  # DST.sb1 <- sb1$biomodTable
  DST.sb2 <- sb2$biomodTable
  DST.sb3 <- sb3$biomodTable
  DST.eb <- eb$biomodTable
  #select the runs which meets the following requirements
  # cols.sb1 <- which(sb1$records[2] > 0.5*pre.records & sb1$records[4] > 0.1*pre.records)
  cols.sb2 <- which(sb2$records[2] > 0.5*pre.records & sb2$records[4] > 0.1*pre.records)
  cols.sb3 <- which(sb3$records[2] > 0.5*pre.records & sb3$records[4] > 0.1*pre.records)
  cols.eb <- which(eb$records[2] > 0.5*pre.records & eb$records[4] > 0.1*pre.records)
  # DST2.sb1 <- DST.sb1[,cols.sb1]
  DST2.sb2 <- DST.sb2[,cols.sb2]
  DST2.sb3 <- DST.sb3[,cols.sb3]
  DST2.eb <- DST.eb[,cols.eb]
  DST.all <- as.data.frame(cbind(DST2.sb2,DST2.sb3,DST2.eb))
  col.DST <- ncol(DST.all)
  for (i in 1:col.DST){
    colnames(DST.all)[i] <- paste0('RUN',i)
    #save data split table
    DST <- DST.all[,i]
    DST <- cbind(pre_pa,DST)
    filename <- paste0(model.group,'_pa',pa.id,'_DST',i,'.csv')
    write.csv(DST,file.path(DST.dir,filename))
    #training presence data (for null model)
    run.pre_pa <- cbind(pre_pa,DST.all[,i])
    train.pre <- filter(run.pre_pa,Species==1 & DST.all[,i] == 'TRUE')
    filename <- paste0(model.group,'_pa',pa.id,'_DST',i,'.csv')
    write.csv(train.pre,paste0(Null.train_pre.dir,'/',filename))
    #test presence data and all pseduo-absence data (for null model)
    test.pre <- filter(run.pre_pa,Species==1 & DST.all[,i] == 'FALSE')
    pa <- filter(run.pre_pa,Species==0)
    rest_data <- rbind(test.pre,pa)
    filename <- paste0(model.group,'_pa',pa.id,'_DST',i,'.csv')
    write.csv(rest_data,paste0(Null.rest.dir,'/',filename))
  }
}

model.groups <- c('GroupA','GroupB','GroupC')
unlink(file.path('outputs',species,'DataSplitTable'), recursive = TRUE)
cat(paste0("===============Running ",species.id," species: ",species,"==============="), file = log_con,sep="\n")
for (model.group in model.groups){
  for (pa.id in 1:20){
    print(paste0(species, model.group,pa.id))
    DST.all <- get_DST(pa.id,model.group)
  }
}
cat(paste0("------------------Have generated DST---------------"), file = log_con,sep="\n")
close(log_con)

