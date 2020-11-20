library(aws.s3)
library(tidyr)
library(glmnet)
library(argparse)
library(data.table)
library(compiler)
source("utility/gxe_util.R")
source("utility/query.R")
source("utility/val_util.R")
source("utility/gws_util.R")

# set up random seed
iseed = 623041
set.seed(iseed)

# read argument
parser <- ArgumentParser()
parser$add_argument("--bucket_dir", required=T)
parser$add_argument("--crop", required=T)
parser$add_argument("--fix_effect", required=T)
parser$add_argument("--env_col", required=T)
parser$add_argument("--pipeline_dir", required = T)
parser$add_argument("--model_def", required = F)
parser$add_argument("--trait_list", required = F)
parser$add_argument("--repeats", required = T)
parser$add_argument("--repeat", required = T)
args <- parser$parse_args()
print(args)

bucket_dir      = args[['bucket_dir']]
crop            = args[['crop']]
env_col         = args[['env_col']]
pipeline_dir    = args[['pipeline_dir']]
traitList       = args[['trait_list']]
fix_effect      = args[['fix_effect']]
modelDef        = args[['model_def']]
JobNum          = as.numeric(args[['repeats']])
iFold           = as.numeric(args[['repeat']])

output_dir      = paste0(pipeline_dir,"/gxe")
consolidate_dir = paste0(pipeline_dir,"/consolidated_raw")
origin_dat_dir  = paste0(consolidate_dir, "/GenoPhenoRaw.Rdata")
year_dat_dir    = paste0(consolidate_dir,"/yearOutGenoPhenoRaw.Rdata")
local_tmp_dir   = "./tmp"

#####
# load data and dependencies

# if traitList is given then use it, otherwise use default
if (is.null(traitList)){
  traits=scan(paste("../resources/", str_to_title(crop),".gws.traits.txt",sep=''), what='character')
} else if (tolower(traitList) == "default"){#default
  traits=scan(paste("../resources/", str_to_title(crop),".gws.traits.txt",sep=''), what='character')
} else {
  reTry('save_object(file = paste(local_tmp_dir,"/",basename(traitList), sep = ""),
        object=traitList, bucket=bucket_dir)')
  traits=scan(paste(local_tmp_dir,"/",basename(traitList), sep = ""), what='character')
}
print("trait list in making training sets")
print(traits)

# load train configuration
trainConf=NULL
if(tolower(modelDef)=='default'){ # default model configureations
  print("default model configurations in resources")
  trainConf = read.csv(paste("../resources/",crop,"_training_set_config.csv",sep=""),  stringsAsFactors=F, colClasses="character")
} else if(grepl("\\.csv$",modelDef)){ # custom model configurations
  print("custom model configurations (csv)")
  reTry('save_object(file = paste(local_tmp_dir,"/",basename(modelDef), sep = ""),
        object=modelDef, bucket=bucket_dir)')
  trainConf=fread(paste(local_tmp_dir,"/",basename(modelDef), sep = ""), stringsAsFactors=F, colClasses="character") %>% as.data.frame()
} else {
  print("custom model configurations (text file, population lists)")
  reTry('save_object(file = paste(local_tmp_dir,"/",basename(modelDef), sep = ""),
        object=modelDef, bucket=bucket_dir)')
  fileList=scan( paste(local_tmp_dir,"/",basename(modelDef), sep = ""), what='character')
  inv=subset(inv, FILE%in%fileList) # limit inventory to the supplied file list
  trainConf=data.frame(NAME="CUSTOM_MODEL",REGION="",COUNTRY='',RM='',YEAR='',GENDER='',SUB_GROUP='',PROG='',EXOTIC_COUNTRY='',stringsAsFactors=F) #make empty training configuration entry
}


# load udr design
########

markerMap$MARKER = as.character(markerMap$MARKER)

logFold=rep(F,JobNum)
logFold[iFold]=T


# set up model hyperparameters
Alpha = 0.001
maxRAM <- 400e6
lambdaSeq=rev((exp(seq(log(0.0001), log(300),(log(300)-log(0.0001))/100))))


# end of load data and dependencies
#####

for(i in (1:nrow(trainConf))){
  
  trainSetCond = trainConf[i,,drop=F]
  modelName = trainSetCond[,'NAME']
  cat("\n========= ",modelName," ==========\n\n")
  
  conslidate_dat_dir  = paste0(consolidate_dir,"/",modelName,".consolidatedGenoPhenoRaw.Rdata")
  origin_dat_dir      = paste0(consolidate_dir,"/",modelName,".GenoPhenoRaw.Rdata")
  year_dat_dir        = paste0(consolidate_dir,"/",modelName,".yearOutGenoPhenoRaw.Rdata")
  
  logFold=rep(F,JobNum)
  logFold[iFold]=T
  
  # running origin out
  GXEprocess(modelName, model_out = paste0(output_dir, "/originOut"),
             dat_dir = origin_dat_dir, bucket_dir = bucket_dir)
  
  
  # running year out
  GXEprocess(modelName, model_out = paste0(output_dir, "/yearOut"),
             dat_dir = year_dat_dir, bucket_dir = bucket_dir)
  
  # running 5-fold
  for (j in 1:5) {
    fold_dat_dir = paste0(consolidate_dir,"/",modelName,'.5foldGenoPhenoRaw.',j,'.Rdata')
    GXEprocess(modelName = modelName, model_out = paste0(output_dir, "/crossVal_fold_",j),
               dat_dir = fold_dat_dir, bucket_dir = bucket_dir)
  }
}


