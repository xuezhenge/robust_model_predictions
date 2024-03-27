library(terra)
library(rgdal)
library(dismo)
library(tools)
library(dplyr)
library(biomod2)
library(tictoc)
library(parallel)
library(doParallel)
library(argparse)
library(raster)

# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--modelgroup", default='GroupA', type="character",
                    help="input the group id, e.g., GroupA")
parser$add_argument("--speciesid", default='5', type="integer",
                    help="input the species id, e.g., 1")
parser$add_argument("--modelname", default='MAXENT', type="character",
                    help="input the model name, e.g., GLM")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
model.group <- args$modelgroup
species.id <- args$speciesid
model.name <- args$modelname
n.cores <- 24

setwd("/home/xge/scratch/correlativeSDM")
#setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM/')
species_list <- list.files('outputs')
species <- species_list[species.id]
#---------Step1: Load climate data -------
#-- historical climate data
hist.dir <- "inputs/world_clim2/world/1971-2000"
sel.evnames <- read.csv(file.path('outputs',species,'sp_bios_names/hist_bio_names.csv'))
sel.evnames <- as.vector(sel.evnames$x)
fns<- paste0(sel.evnames,".tif")
climate.dir <- paste0(hist.dir,'/',fns)
bios.hist <- raster::stack(climate.dir)

#---------Step2: Load presence - pseudo absence data -------
get_random.prepa.data <- function(species,model.group,fn.id,random.id){
  prepa.dir<- file.path('outputs',species,"pre_pa_rp",model.group)
  filelist <- list.files(prepa.dir)
  filename <- filelist[fn.id]
  prepa.all <- read.csv(file.path(prepa.dir,filename))
  random.pre <- filter(prepa.all,random == random.id)
  rp <- data.frame(random.pre$x,random.pre$y)
  pa <- filter(prepa.all,random == 0)
  random.pre.pa <- rbind(random.pre,pa)
  return(random.pre.pa)
}

get_real.prepa.data <- function(species,model.group,fn.id){
  setwd("/home/xge/scratch/correlativeSDM")
  #setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM/')
  prepa.dir<- file.path('outputs',species,'DataSplitTable',model.group)
  filelist <- list.files(prepa.dir)
  filename <- filelist[fn.id]
  prepa <- read.csv(file.path(prepa.dir,filename))
  return(prepa)
}

get_BiomodData <- function(pre.pa,bios.hist){
  myResp <- as.numeric(pre.pa$Species)
  myRespXY <- pre.pa[, c('x', 'y')]
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = bios.hist, # explanatory raster data
                                       resp.xy = myRespXY,
                                       resp.name = species,
                                       na.rm = TRUE)
  return(myBiomodData)
}

# 3. define individual models options ---- 
Sp_opt <- 
  BIOMOD_ModelingOptions(
  )

get_om.test <- function(species,pre.pa,pre.test,model.name,bios.hist,fn.id){
  myBiomodData <- get_BiomodData(pre.pa,bios.hist)
  DST <- as.matrix(pre.pa$DST)
  if (model.name == 'MAXENT'){
    model.name = 'MAXENT.Phillips.2'
  }
  setwd("/home/xge/scratch/correlativeSDM/BiomodOutput")
  #setwd("/Users/xuezhenge/Desktop/application_paper/correlativeSDM/BiomodOutput")
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                      modeling.id = paste0(model.name,fn.id),
                                      models = c(model.name),
                                      bm.options = Sp_opt,
                                      data.split.table = DST,
                                      metric.eval = c('POD'),
                                      var.import = 0,
                                      prevalence = 0.5,
                                      do.full.models = FALSE,
                                      scale.models = TRUE,
                                      save.output = FALSE,
                                      nb.cpu = 1)
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    proj.name = paste0('hist',model.name,fn.id),
                                    new.env = bios.hist,
                                    models.chosen = 'all',
                                    metric.binary = 'all',
                                    metric.filter = 'all',
                                    build.clamping.mask = FALSE)
  pre.train <- filter(pre.pa, DST == 'TRUE' & Species == '1')
  #get coordinates
  pre.train.coord <- data.frame(pre.train$x,pre.train$y)
  pre.test.coord <- data.frame(pre.test$x,pre.test$y)
  #get projected value for each train and test data
  proj.value <- myBiomodProj@proj.out@val[[1]]
  pre.train.proj.value <- raster::extract(proj.value,pre.train.coord)
  pre.test.proj.value <- raster::extract(proj.value,pre.test.coord)
  pre.train.om10 <- quantile(pre.train.proj.value,0.1)
  pre.train.om5 <- quantile(pre.train.proj.value,0.05)
  pre.train.om0 <- quantile(pre.train.proj.value,0)
  t <- length(pre.test.proj.value)
  om10.test = length(pre.test.proj.value[pre.test.proj.value < pre.train.om10])/t
  om5.test = length(pre.test.proj.value[pre.test.proj.value < pre.train.om5])/t
  om0.test = length(pre.test.proj.value[pre.test.proj.value < pre.train.om0])/t
  POD <- get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  om_POD <- 1 - POD$Testing.data
  Th_POD <- POD$Cutoff
  oms <- data.frame(fn.id,om_POD,om0.test,om5.test,om10.test,Th_POD, pre.train.om0, pre.train.om5, pre.train.om10)
  return(oms)
}

get_real_oms <- function(fn.id){
  fn_om <- filelist[fn.id]
  om_outfile.dir <- file.path(real.all.dir,fn_om)
  if(file.exists(om_outfile.dir)==FALSE){
    real.pre.pa <- get_real.prepa.data(species,model.group,fn.id)
    real.pre.test <- filter(real.pre.pa, Species == '1' & DST == 'FALSE')
    oms <- get_om.test(species,real.pre.pa,real.pre.test,model.name,bios.hist,fn.id)
    names(oms) <- c('fnid','ompod','om0','om5','om10','thpod','th0','th5','th10')
    setwd("/home/xge/scratch/correlativeSDM")
    #setwd('/Users/xuezhenge/Desktop/application_paper/correlativeSDM')
    write.csv(oms,om_outfile.dir)
    cat(paste0("---",model.name,fn.id," Done!----"), file = log_con,sep="\n")
  }else{cat(paste0("---",model.name,fn.id," exists!----"), file = log_con,sep="\n")}
}

# output folder
if(!exists("log_files/model_projection/")) dir.create("log_files/model_projection/",recursive = TRUE)
log_con <- file(paste0("log_files/model_projection/",species.id,species,".txt"), open="a")
cat(paste0("---",species.id,species,"-",model.name, "----"), file = log_con,sep="\n") 
real.all.dir <- file.path('outputs',species,'omission_rates',model.name)
if(!exists(real.all.dir)) dir.create(real.all.dir,recursive = TRUE)
# do parallel
prepa.dir<- file.path('outputs',species,'DataSplitTable',model.group)
filelist <- list.files(prepa.dir)
len <- length(filelist)
reps = seq(1:len)
mclapply(reps,get_real_oms, mc.cores = n.cores)
cat(paste0("---",model.name," model hist Done!----"), file = log_con,sep="\n")
close(log_con)






