
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(argparse)
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--speciesid", default='3', type="integer",
                    help="input the species id, e.g., 1")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
species.id <- args$speciesid

setwd("/home/xge/scratch/correlativeSDM")
species_list <- list.files('outputs')
species <- species_list[species.id]
modelnames = c("GLM", "GBM", "GAM", "CTA","FDA", "MARS", "RF","MAXENT")
#modelnames = c("GBM", "GAM", "CTA","FDA", "MARS", "RF","MAXENT")

if(!exists("log_files/model_projection/")) dir.create("log_files/model_projection/",recursive = TRUE)
log_con <- file(paste0("log_files/model_projection/",species.id,species,".txt"), open="a")

om_alls <- {}
for (model.name in modelnames){
  print(model.name)
  om.dir<- file.path('outputs',species,"omission_rates",model.name)
  fns <- list.files(om.dir)
  len <- length(fns)
  for (i in 1:len){
    om_data <- read.csv(file.path(om.dir, fns[i]))
    om_data$model <- model.name
    om_alls <- rbind(om_alls,om_data)
  }
  ompods <- filter(om_alls, ompod == 0)
  om0s <- filter(om_alls, om0 == 0)
  om5s <- filter(om_alls, om5 <= 0.05)
  om10s <- filter(om_alls,om10 <= 0.1)
  
  #write summary table
  out.dir<- file.path('outputs',species,'summary_doc','omrates.csv')
  write.csv(om_alls,out.dir)
  
  # check how many runs were selected from each model category
  model_no <- count(om_alls,'model')
  select_om10_no <- count(om10s, 'model')
  comparisons <- merge(model_no, select_om10_no, all = TRUE,by = "model")
  names(comparisons) <- c('modelname','freq_all','freq_om10')
  comparisons$ratio <- comparisons$freq_om10/comparisons$freq_all
  out.dir<- file.path('outputs',species,'summary_doc','selected_runs_permodel.csv')
  write.csv(comparisons,out.dir)
  cat(paste0("---",species," sum omission rate tables Done!----"), file = log_con,sep="\n") 
}

