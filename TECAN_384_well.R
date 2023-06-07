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
  input_path,
  well_ids
  ){
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

  #Default well ids (all non-NA names)
  if(length(well_ids)==1 & well_ids[1] == "all"){
    well_ids <- design_dt[!is.na(name),id]
  }

  #Wells ids or names supplied?
  well_names <- FALSE
  if(sum(well_ids %in% design_dt[,name])>sum(well_ids %in% design_dt[,id])){
    well_names <- TRUE
  }

  #Check well ids are valid 
  if(sum(well_ids %in% design_dt[,id])!=length(well_ids) & !well_names){
    stop("Invalid well ids specified.", call. = FALSE)
  }

  #Check well names are valid 
  if(sum(well_ids %in% design_dt[,name])!=length(well_ids) & well_names){
    stop("Invalid well names specified.", call. = FALSE)
  }

  #Translate well names to well ids
  if(well_names){
    well_ids <- design_dt[name %in% well_ids,id]
  }

  return(list(
    'design' = design_dt,
    'well_ids' = well_ids))
}

#Import TECAN plate data
import_tecan_plate_data <- function(
  input_path,
  well_ids
  ){
  ## import TECAN results (all wells) from summary excel sheet
  suppressMessages(suppressWarnings(tecan <- as.matrix(read_xlsx(excel_path, sheet = 1))))
  tecan <- tecan[22:nrow(tecan),]
  rownames(tecan) <- tecan[,1]
  tecan <- tecan[,-1]

  # tecan_select <- tecan_data[c('D2', 'G23', 'N18'),]
  tecan_select <- tecan[well_ids,]
  colnames(tecan_select) <- as.numeric(tecan[1,])/3600 ### convert time to hours
  class(tecan_select) <- 'numeric' ### convert table to numeric format

  return(tecan_select)
}

#Calculate growth rate results using heuristic approach from growthrates package i.e. similar to the "growth rates made easy"-method of Hall et al. (2013)
calculate_growth_rates_heuristic <- function(
  tecan_data,
  plate_design,
  h_parameter
  ){

  ## determine each well's max. exponential rates by scanning for the log-linear growth range
  ## h-parameter: here set to 15 consecutive time points, i.e. growth looking at the best window of 3-4h, feel free to vary/reduce
  well_ids <- rownames(tecan_data)
  results_list <- list()
  for (i in 1:length(well_ids)){
    
    ## summarise well data and run fit
    tmp <- rbind(as.numeric(colnames(tecan_data)),
                 as.numeric(tecan_data[i,]))
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
parser <- add_argument(parser, "--method", default='heuristic', help = "Maximum growth rate method")
parser <- add_argument(parser, "--parameter", type = "integer", default=15, help = "h-parameter; number of consecutive time points to evaluate maximum growth rate")
parser <- add_argument(parser, "--wells", default='all', help = "Comma-separated list of well ids or well names")
parser <- add_argument(parser, "--outputPath", help = "Output file path (default: no output file; print results to stdout)")

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
method <- args[['method']]
hParameter <- args[['parameter']]
wells <- unlist(strsplit(args[['wells']], ','))
outputPath <- args[['outputPath']]

###########################
### IMPORT PLATE DESIGN
###########################

tecan_design <- import_tecan_plate_design(
  input_path = excel_path,
  well_ids = wells)
design_dt <- tecan_design[['design']]
wells <- tecan_design[['well_ids']]

###########################
### IMPORT PLATE DATA MATRIX
###########################

tecan_mat <- import_tecan_plate_data(
  input_path = excel_path,
  well_ids = wells)

###########################
### CALCULATE GROWTH RATES
###########################

if(method == "heuristic"){
  result_dt <- calculate_growth_rates_heuristic(
    tecan_data = tecan_mat,
    plate_design = design_dt,
    h_parameter = hParameter)
}

###########################
### PRINT OR SAVE RESULTS
###########################

#Print or save
if(is.na(outputPath)){
  print(result_dt)
}else{
  write.table(result_dt, file = outputPath, sep = "\t", quote = F, row.names = F)
}




