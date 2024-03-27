library(terra)
library(rgdal)
library(dismo)
library(tools)
library(dplyr)
library(biomod2)
library(RColorBrewer)
library(argparse)

# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--speciesid", default='1', type="integer",
                    help="input the species id, e.g., 1")

args <- parser$parse_args()
species.id <- args$speciesid

setwd("/home/xge/scratch/correlativeSDM")
species_list <- list.files('outputs')
species <- species_list[species.id]
speciesV2 <- gsub("_", ".", species)

time.periods <- c('2021-2040','2081-2100')
scenarios <- c('ssp585')
om_idxs <- c(10)
# get model projections under different conditions
get_projections <- function(species,time.period,scenario,om_sel,i){
  setwd("/home/xge/scratch/correlativeSDM/BiomodOutput")
  fn.id <- om_sel$fnid[i]
  model.name <- om_sel$model[i]
  if (model.name == 'MAXENT'){
    model.name = 'MAXENT.Phillips.2'
  }
  if (scenario == 'hist'){
    biomodoutput.name <- paste0(speciesV2,'.',scenario,model.name,fn.id,'.projection.out')
    biomodoutput.dir <- paste0('proj_',scenario,model.name,fn.id)
  }else{
    biomodoutput.name <- paste0(speciesV2,'.',time.period,scenario,model.name,fn.id,'.projection.out')
    biomodoutput.dir <- paste0('proj_',time.period,scenario,model.name,fn.id) 
  }
  biomodoutput.path <- file.path(speciesV2,biomodoutput.dir,biomodoutput.name)
  load(biomodoutput.path)
  myBiomodModelProj <- get(load(biomodoutput.path))
  myBiomodModelProj <- myBiomodModelProj@proj.out@val
  setwd("/home/xge/scratch/correlativeSDM")
  return(myBiomodModelProj)
}

get_ensplot <- function(species,time.period,scenario,om_idx){
  setwd("/home/xge/scratch/correlativeSDM")
  #create outraster folder
  outraster.dir = paste0('outrasters/',species,'/',time.period,'_',scenario,'_',om_idx)
  if(!exists(outraster.dir)) dir.create(outraster.dir,recursive = TRUE)
  fn = paste0('omth',om_idx,'.tif')
  outraster <- file.path(outraster.dir,fn)
  plot.dir <- paste0('outplots/',species,'/',time.period,'_',scenario,'/omth',om_idx,'.png')
  if (file.exists(plot.dir)==TRUE){
    cat(paste0(time.period,' ', scenario, ' ',om_idx,' Plot existed!'), file = log_con,sep="\n") 
  }else{
    #select the runs that meet the om requirement
    om_alls <- read.csv(file.path('outputs',species,'summary_doc','omrates.csv'))
    if (om_idx == 10){
      om_sel <- filter(om_alls,om10 <= 0.1)
      }else if(om_idx == 5){
      om_sel <- filter(om_alls,om5 <= 0.05)
      }else if(om_idx == 0){
        om_sel <- filter(om_alls,om0 == 0)
      }else if(om_idx == 'pod'){
        om_sel <- filter(om_alls,ompod == 0)
      }else if(om_idx == 'all'){
        om_sel <- om_alls
      }
    #stack these runs
    len <- nrow(om_sel)
    myBiomodModelProjs <- get_projections(species,time.period,scenario,om_sel,1)
    if (len > 1){
      for (i in 2:len){
        myBiomodModelProj <- get_projections(species,time.period,scenario,om_sel,i)
        myBiomodModelProjs <- stack(myBiomodModelProjs,myBiomodModelProj)
      }
    }
    #rescale the projections and get the mean
    ens <- mean(myBiomodModelProjs)/1000
    #save the ensemble projection
    # save raster
    writeRaster(ens, filename=outraster, overwrite=TRUE)
    # save plot
    outplot.dir = paste0('outplots/',species,'/',time.period,'_',scenario)
    if(!exists(outplot.dir)) dir.create(outplot.dir,recursive = TRUE)
    fn = paste0('omth',om_idx,'.png')
    outplot = file.path(outplot.dir,fn)
    cols <- brewer.pal(9, "OrRd")
    pal <- colorRampPalette(cols)
    png(file=outplot, width=600*5, height=300*5,res=300)
    plot(ens,zlim = c(0,1),col=pal(50),main=paste0(species,' ',time.period,' ',scenario,' omth',om_idx))
    dev.off()
    cat(paste0(time.period,' ', scenario, ' ',om_idx,'Plot saved!'), file = log_con,sep="\n") 
  }
}

# run these functions
if(!exists("log_files/model_projection/")) dir.create("log_files/model_projection/",recursive = TRUE)
log_con <- file(paste0("log_files/model_projection/",species.id,species,".txt"), open="a")
cat(paste0('------------Species: ',species,' ensemble-------------'), file = log_con,sep="\n")
get_ensplot(species,'ensmean','hist','all')
get_ensplot(species,'ensmean','hist','pod')
for (om_idx in om_idxs){
  print(paste0('Select the runs based on om',om_idx))
  get_ensplot(species,'ensmean','hist',om_idx)
  print('1971-2000 Plot saved!')
  for (time.period in time.periods){
    for (scenario in scenarios){
      get_ensplot(species,time.period,scenario,om_idx)
    }
  }
}



