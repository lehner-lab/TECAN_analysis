#!/usr/bin/env Rscript

## TECAN experiment analysis
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
  "data.table",
  "ggplot2",
  "reshape2")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)!=0){
  stop(paste0("Required R packages not installed. Please install the following packages: ", paste(missing_packages, sep = ", ")), call. = FALSE)
}

###########################
### FUNCTIONS
###########################

# #Import TECAN plate design - obsolete
# import_tecan_plate_design <- function(
#   input_path,
#   well_ids
#   ){
#   ## import 384-well plate design matrix from summary excel sheet
#   suppressMessages(suppressWarnings(design <- as.matrix(read_xlsx(input_path, sheet = 1))[2:18,1:25]))
#   colnames(design) <- design[1,]
#   design <- design[-1,]
#   rownames(design) <- design[,1]
#   design <- design[,-1]
#   design_dt <- data.table(
#     id = paste0(rep(rownames(design), each = dim(design)[2]), as.integer(colnames(design))),
#     name = as.character(t(design)))

#   #Check design matches that of 384-well plate 
#   if(length(as.character(design))!=384 | nrow(design_dt)!=384){
#     stop("Excel file format invalid for 384-well plate.", call. = FALSE)
#   }

#   #Default well ids (all non-NA names)
#   if(length(well_ids)==1 & well_ids[1] == "all"){
#     well_ids <- design_dt[!is.na(name),id]
#   }

#   #Wells ids or names supplied?
#   well_names <- FALSE
#   if(sum(well_ids %in% design_dt[,name])>sum(well_ids %in% design_dt[,id])){
#     well_names <- TRUE
#   }

#   #Check well ids are valid 
#   if(sum(well_ids %in% design_dt[,id])!=length(well_ids) & !well_names){
#     stop("Invalid well ids specified.", call. = FALSE)
#   }

#   #Check well names are valid 
#   if(sum(well_ids %in% design_dt[,name])!=length(well_ids) & well_names){
#     stop("Invalid well names specified.", call. = FALSE)
#   }

#   #Translate well names to well ids
#   if(well_names){
#     well_ids <- design_dt[name %in% well_ids,id]
#   }

#   return(list(
#     'design' = design_dt,
#     'well_ids' = well_ids))
# }

#Import TECAN plate data
import_tecan_plate_data <- function(
  input_path,
  well_ids
  ){
  ## import TECAN results (all wells) from summary excel sheet
  suppressMessages(suppressWarnings(tecan <- as.matrix(read_xlsx(excel_path, sheet = 1))))
  row1 <- grep("Cycle Nr.", tecan[,1])[1]+1
  tecan <- tecan[row1:nrow(tecan),]
  rowN <- which(is.na(tecan[,1]))[1]-1
  if(!is.na(rowN)){
    tecan <- tecan[1:rowN,]
  }
  rownames(tecan) <- tecan[,1]
  tecan <- tecan[,-1]

  #Default well ids (all non-NA names)
  if(length(well_ids)==1 & well_ids[1] == "all"){
    well_ids <- rownames(tecan)[!grepl("\\[", rownames(tecan))]
  }

  #Check well ids are valid 
  if(sum(well_ids %in% rownames(tecan))!=length(well_ids)){
    stop("Invalid well ids specified.", call. = FALSE)
  }

  # tecan_select <- tecan_data[c('D2', 'G23', 'N18'),]
  tecan_select <- tecan[well_ids,]
  colnames(tecan_select) <- as.numeric(tecan[1,])/3600 ### convert time to hours
  class(tecan_select) <- 'numeric' ### convert table to numeric format

  return(tecan_select)
}

#Calculate growth rate results using heuristic approach from growthrates package i.e. similar to the "growth rates made easy"-method of Hall et al. (2013)
calculate_growth_rates_heuristic <- function(
  tecan_data,
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
      'maxGR' = tmp@par[['mumax']],
      'lag' = tmp@par[['lag']])
  }
  return(rbindlist(results_list))
}

#Set well status based on dead and lag thresholds
set_well_status <- function(
  input_dt,
  dead_threshold,
  lag_threshold
  ){
  input_dt[, status := 'PASS']
  input_dt[maxGR<dead_threshold, status := paste0('FAIL (maxGR<', dead_threshold, ')')]
  input_dt[lag>lag_threshold, status := paste0('FAIL (lag>', lag_threshold, ')')]
  return(input_dt)
}

#ggplot2-like colour scale in HCL space.
gg_color_hue <- function(
  n
  ){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Plot growth curves
plot_growth_curves <- function(
  tecan_data,
  input_dt){
  #Format plot data
  plot_dt <- as.data.table(t(tecan_data))
  plot_dt[, time := as.numeric(colnames(tecan_data))]
  plot_dt <- reshape2::melt(plot_dt, id = c('time'))
  plot_dt <- as.data.table(plot_dt)[!is.na(value)]
  colnames(plot_dt)[2:3] <- c('well_id', 'OD')
  #Well status
  plot_dt <- merge(plot_dt, input_dt[,.(well_id, status)], by = 'well_id')
  plot_dt[, status_id := paste0(status, ": ", well_id)]
  plot_dt <- plot_dt[order(status_id, decreasing = T)]
  status_id_levels <- plot_dt[!duplicated(status_id),status_id]
  plot_dt[, status_id := factor(status_id, levels = status_id_levels)]
  #Plot colours
  plot_cols <- gg_color_hue(3)
  plot_dt[,status_col := plot_cols[2]]
  plot_dt[grep("maxGR", status), status_col := plot_cols[1]]
  plot_dt[grep("lag", status), status_col := plot_cols[3]]
  plot_dt[, status_col := factor(status_col, levels = c(plot_cols[2], plot_cols[1], plot_cols[3]))]
  plot_cols <- plot_dt[,as.character(status_col)]
  names(plot_cols) <- plot_dt[, status_id]
  plot_dt[, status := factor(status, levels = unique(status))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(time, OD, color = status_id)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::xlab("Time (h)") +
    ggplot2::ylab("Optical density") +
    ggplot2::facet_wrap(~status, ncol = 1) +
    ggplot2::theme_bw()
  if(input_dt[status!="PASS",.N]>20){
    d <- d + ggplot2::scale_colour_manual(values = plot_cols, breaks = names(plot_cols)[!grepl("PASS", names(plot_cols))], guide="none")
  }else{
    d <- d + ggplot2::scale_colour_manual(values = plot_cols, breaks = names(plot_cols)[!grepl("PASS", names(plot_cols))])
  }
  ggplot2::ggsave(paste0(outputPrefix, '.pdf'), d, width = 8, height = 8, useDingbats=FALSE)
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
parser <- add_argument(parser, "--wells", default='all', help = "Comma-separated list of well ids")
parser <- add_argument(parser, "--deadThreshold", type = "double", default=0.05, help = "Growth rate threshold for dead variants")
parser <- add_argument(parser, "--lagThreshold", type = "double", default=48.0, help = "Lag time threshold for problematic variants")
parser <- add_argument(parser, "--outputPrefix", help = "Output path prefix (default: no output file; print results to stdout)")

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
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(reshape2)))

###########################
### GLOBALS
###########################

#Excel file path
excel_path <- args[['excel_path']]
method <- args[['method']]
hParameter <- args[['parameter']]
wells <- unlist(strsplit(args[['wells']], ','))
dead_threshold <- args[['deadThreshold']]
lag_threshold <- args[['lagThreshold']]
outputPrefix <- args[['outputPrefix']]

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
    h_parameter = hParameter)
}

###########################
### SET WELL STATUS
###########################

result_dt <- set_well_status(
    input_dt = result_dt,
    dead_threshold = dead_threshold,
    lag_threshold = lag_threshold)

###########################
### PLOT GROWTH CURVES
###########################

plot_growth_curves(
  tecan_data = tecan_mat,
  input_dt = result_dt)

###########################
### PRINT OR SAVE RESULTS
###########################

#Print or save
if(is.na(outputPrefix)){
  print(result_dt)
}else{
  write.table(result_dt, file = paste0(outputPrefix, '.txt'), sep = "\t", quote = F, row.names = F)
}




