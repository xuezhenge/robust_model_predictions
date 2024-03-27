library(terra)
library(raster)
library(rgdal)
library(dismo)
library(tools)
library(dplyr)
library(biomod2)
library(tictoc)
library(parallel)
library(doParallel)
library(argparse)

# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--speciesid", default='1', type="integer",
                    help="input the species id, e.g., 1")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
species.id <- args$speciesid
n.cores <- 24

setwd("/home/xge/scratch/correlativeSDM")
#setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM/')
species_list <- list.files('outputs')
species <- species_list[species.id]
speciesV2 <- gsub("_", ".", species)

#---------Step 1: get climate------
get_climate <- function(time.period,scenario){
  setwd("/home/xge/scratch/correlativeSDM")
  #setwd("/Users/xuezhenge/Desktop/application_paper/correlativeSDM")
  data.dir <- file.path('inputs/world_clim2/world',time.period)
  sel.evnames <- read.csv(file.path('outputs',species,'sp_bios_names/hist_bio_names.csv'))
  sel.evnames <- as.vector(sel.evnames$x)
  fns<- paste0(sel.evnames,".tif")
  climate.dir <- paste0(data.dir,'/',fns)
  bios <- raster::stack(climate.dir)
  return(bios)
}

#---------Step 2: do predictions------
get_projections <- function(i){
  setwd("/home/xge/scratch/correlativeSDM/BiomodOutput")
  #setwd("/Volumes/XG_data/state_species/BiomodOutput")
  #setwd("/Users/xuezhenge/Desktop/state_species/BiomodOutput")
  fn.id <- om_sel$fnid[i]
  model.name <- om_sel$model[i]
  biomodoutput.name <- paste0(speciesV2,'.',time.period,scenario,model.name,fn.id,'.projection.out')
  biomodoutput.dir <- paste0('proj_',time.period,scenario,model.name,fn.id)
  biomodoutput.path <- file.path(speciesV2,biomodoutput.dir,biomodoutput.name)
  if (file.exists(biomodoutput.path)==TRUE){
     cat(paste0(biomodoutput.name,' file exists!!'), file = log_con,sep="\n") 
  }else{
  bios <- get_climate(time.period,scenario)
  setwd("/home/xge/scratch/correlativeSDM/BiomodOutput")
  if (model.name == 'MAXENT'){
    model.name = 'MAXENT.Phillips.2'
  }
  biomodoutput.name <- paste0(speciesV2,'.',model.name,fn.id,'.models.out')
  biomodoutput.path <- file.path(speciesV2,biomodoutput.name)
  load(biomodoutput.path)
  myBiomodModelOut <- get(load(biomodoutput.path))
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    proj.name = paste0(time.period,scenario,model.name,fn.id),
                                    new.env = bios,
                                    build.clamping.mask = FALSE)
  cat(paste0(biomodoutput.name,' file just done!!'), file = log_con,sep="\n") 
  }
}

om_alls <- read.csv(file.path('outputs',species,'summary_doc','omrates.csv'))
# om_sel <- filter(om_alls,om10 <= 0.1 | om5 <= 0.05 | om0 == 0)
om_sel <- filter(om_alls,om10 <= 0.1)
len <- nrow(om_sel)
reps = seq(1:len)

if(!exists("log_files/model_projection/")) dir.create("log_files/model_projection/",recursive = TRUE)
log_con <- file(paste0("log_files/model_projection/",species.id,species,"fut.txt"), open="a")
cat(paste0("---",species.id, species," Start Future Projection:----"), file = log_con,sep="\n") 

time.periods <- c('2021-2040','2081-2100')
scenarios <- c('ssp585')
for (time.period in time.periods){
  for (scenario in scenarios){
    cat(paste0("---",time.period,' ', scenario," Running!----"), file = log_con,sep="\n") 
    mclapply(reps,get_projections, mc.cores = n.cores)
    toc()
    cat(paste0("---",time.period,' ', scenario," Done!----"), file = log_con,sep="\n")
  }
}
close(log_con)
