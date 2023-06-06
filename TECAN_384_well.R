#!/usr/bin/env Rscript

## Example TECAN experiment analysis
## 06/06/2023
## maximilian.stammnitz@crg.eu
## andre.faure@crg.eu

###########################
### CHECK DEPENDENCIES
###########################

#R packages
required_packages <- c(
  "argparser",
  "readxl",
  "growthrates",
  "beeswarm",
  "scales",
  "data.table")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)!=0){
  stop(paste0("Required R packages not installed. Please install the following packages: ", paste(missing_packages, sep = ", ")), call. = FALSE)
}

###########################
### FUNCTIONS
###########################

#Import TECAN plate design
import_tecan_plate_design <- function(
  input_path){
  ## import 384-well plate design matrix from summary excel sheet
  suppressMessages(suppressWarnings(design <- as.matrix(read_xlsx(input_path, sheet = 1))[2:18,1:25]))
  colnames(design) <- design[1,]
  design <- design[-1,]
  rownames(design) <- design[,1]
  design <- design[,-1]
  design_dt <- data.table(
    id = paste0(rep(rownames(design), each = dim(design)[2]), as.integer(colnames(design))),
    name = as.character(t(design)))

  #Check design matches that of 384-well plate 
  if(length(as.character(design))!=384 | nrow(design_dt)!=384){
    stop("Excel file format invalid for 384-well plate.", call. = FALSE)
  }

  return(design_dt)
}

#Import TECAN plate data
import_tecan_plate_data <- function(
  input_path){
  ## import TECAN results (all wells) from summary excel sheet
  suppressMessages(suppressWarnings(tecan <- as.matrix(read_xlsx(excel_path, sheet = 1))))
  tecan <- tecan[22:nrow(tecan),]
  rownames(tecan) <- tecan[,1]
  tecan <- tecan[,-1]
  return(tecan)
}

#Calculate growth rate results
calculate_growth_rates <- function(
  plate_design,
  well_ids,
  tecan_data,
  h_parameter
  ){

  # tecan_select <- tecan_data[c('D2', 'G23', 'N18'),]
  tecan_select <- tecan_data[well_ids,]
  colnames(tecan_select) <- as.numeric(tecan_data[1,])/3600 ### convert time to hours
  class(tecan_select) <- 'numeric' ### convert table to numeric format

  ## determine each well's max. exponential rates by scanning for the log-linear growth range
  ## h-parameter: here set to 15 consecutive time points, i.e. growth looking at the best window of 3-4h, feel free to vary/reduce
  results_list <- list()
  for (i in 1:length(well_ids)){
    
    ## summarise well data and run fit
    tmp <- rbind(as.numeric(colnames(tecan_select)),
                 as.numeric(tecan_select[i,]))
    tmp <- tmp[,!is.na(tmp[2,])]
    tmp <- fit_easylinear(time = tmp[1,], y = tmp[2,], h = h_parameter) ### adjust h paramater

    ### you can run plot(tmp) to check if the slope calculation looks good
    results_list[[well_ids[i]]] <- data.table(
      'well_id' = well_ids[i],
      'well_name' = unlist(plate_design[id==well_ids[i],name]),
      'maxGR' = tmp@par[['mumax']],
      'lag' = tmp@par[['lag']])
  }
  return(rbindlist(results_list))
}


###########################
### COMMAND-LINE OPTIONS
###########################

library(argparser)

#Create Argument Parser
parser <- arg_parser(description = "Your tool description")

#Add Positional Argument for File Path
parser <- add_argument(parser, "excel_path", help = "Path to the Excel file")

#Add Optional Argument for Integer
parser <- add_argument(parser, "--parameter", type = "integer", default=15, help = "h-parameter; number of consecutive time points to evaluate maximum growth rate (default:15)")
parser <- add_argument(parser, "--wells", default='all', help = "Comma-separated list of well ids or well names (default:'all')")
parser <- add_argument(parser, "--outputPath", help = "Output file path (default: print results to stdout)")

#Parse the Command Line Arguments
args <- parse_args(parser)

###########################
### INPUT VALIDATION
###########################

#Check if Excel file path is provided
if(is.null(args[['excel_path']])){
  stop("Excel file path is required.", call. = FALSE)
}else if(!file.exists(args[["excel_path"]])){
  stop("Excel file path does not exist", call. = FALSE)
}

###########################
### DEPENDENCIES
###########################

suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(growthrates)))
suppressMessages(suppressWarnings(library(beeswarm)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(data.table)))

###########################
### GLOBALS
###########################

#Excel file path
excel_path <- args[['excel_path']]
hParameter <- args[['parameter']]
wells <- unlist(strsplit(args[['wells']], ','))
outputPath <- args[['outputPath']]
well_names <- FALSE

###########################
### IMPORT PLATE DESIGN
###########################

design_dt <- import_tecan_plate_design(input_path = excel_path)

#Default well ids (all non-NA names)
if(length(wells)==1 & wells[1] == "all"){
  wells <- design_dt[!is.na(name),id]
}

#Wells ids or names supplied?
if(sum(wells %in% design_dt[,name])>sum(wells %in% design_dt[,id])){
  well_names <- TRUE
}

#Check well ids are valid 
if(sum(wells %in% design_dt[,id])!=length(wells) & !well_names){
  stop("Invalid well ids specified.", call. = FALSE)
}

#Check well names are valid 
if(sum(wells %in% design_dt[,name])!=length(wells) & well_names){
  stop("Invalid well names specified.", call. = FALSE)
}

#Translate well names to well ids
if(well_names){
  wells <- design_dt[name %in% wells,id]
}

###########################
### IMPORT PLATE DATA
###########################

tecan <- import_tecan_plate_data(input_path = excel_path)

###########################
### CALCULATE GROWTH RATES
###########################

result_dt <- calculate_growth_rates(
  plate_design = design_dt,
  well_ids = wells,
  tecan_data = tecan,
  h_parameter = hParameter)

###########################
### PRINT OR SAVE RESULTS
###########################

#Print or save
if(is.na(outputPath)){
  print(result_dt)
}else{
  write.table(result_dt, file = outputPath, sep = "\t", quote = F, row.names = F)
}




