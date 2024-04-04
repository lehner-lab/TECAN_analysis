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
    tmp_fit <- fit_easylinear(time = tmp[1,], y = tmp[2,], h = h_parameter) ### adjust h paramater

    ### you can run plot(tmp) to check if the slope calculation looks good
    results_list[[well_ids[i]]] <- data.table(
      'well_id' = well_ids[i],
      'maxGR' = tmp_fit@par[['mumax']],
      'lag' = tmp_fit@par[['lag']],
      'maxOD' = max(tmp[2,]))
  }
  return(rbindlist(results_list))
}

#Set well status based on dead and lag thresholds
set_well_status <- function(
  input_dt,
  dead_threshold,
  lag_threshold,
  design_path
  ){
  input_dt[, status := 'PASS']
  input_dt[maxGR<dead_threshold | maxOD<od_threshold, status := paste0('FAIL (maxGR<', dead_threshold, ')')]
  input_dt[lag>lag_threshold, status := paste0('FAIL (lag>', lag_threshold, ')')]

  #Check if design path exists
  if(is.null(design_path)){
    return(input_dt)
  }else if(!file.exists(design_path)){
    return(input_dt)
  }

  #Load design
  design_dt <- fread(design_path)
  #Check if Well column exists
  if(!"Well" %in% names(design_dt)){
    stop("Well column missing in plate design file ('--designPath').", call. = FALSE)
  }
  #Check all well_ids exist in Well column
  if(sum(!input_dt[,well_id] %in% design_dt[,Well])!=0){
    stop("One of more wells missing in plate design file ('--designPath').", call. = FALSE)
  }
  #Check if Plasmid column exists
  if(!"Plasmid" %in% names(design_dt)){
    print("Warning: Plasmid column missing in plate design file ('--designPath'). Plasmid-specific growth curves and boxplots will not be produced.", call. = FALSE)
  }else if(!"Well class" %in% names(design_dt)){
    print("Warning: 'Well class' column missing in plate design file ('--designPath'). Plasmid-specific boxplots will not be produced.", call. = FALSE)
  }
  design_dt <- design_dt[,.SD,,.SDcols = names(design_dt)[names(design_dt) %in% c("Well", "Plasmid", "Well class")]]

  #Merge results with design
  well_ids_sorted <- input_dt[,well_id]
  names(design_dt)[names(design_dt)=='Well'] <- 'well_id'
  names(design_dt) <- gsub(" ", "_", names(design_dt))
  input_dt <- merge(input_dt, design_dt[], by = 'well_id', all.x = T)

  #Sort by well ids
  input_dt <- input_dt[, well_id_factor := factor(well_id, levels = well_ids_sorted)][order(well_id_factor)]
  input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="well_id_factor"]]

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
    ggplot2::geom_line(size = 0.5) +
    ggplot2::xlab("Time (h)") +
    ggplot2::ylab("Optical density") +
    ggplot2::facet_wrap(~status, ncol = 1) +
    ggplot2::theme_bw()
  if(input_dt[status!="PASS",.N]>20){
    d <- d + ggplot2::scale_colour_manual(values = plot_cols, breaks = names(plot_cols)[!grepl("PASS", names(plot_cols))], guide="none")
  }else{
    d <- d + ggplot2::scale_colour_manual(values = plot_cols, breaks = names(plot_cols)[!grepl("PASS", names(plot_cols))])
  }
  ggplot2::ggsave(paste0(outputPrefix, '.pdf'), d, width = plot_width, height = plot_height, useDingbats=FALSE)

  #Individual traces per plasmid
  plasmid_dir <- paste0(outputPrefix, "_plasmids")
  if('Plasmid' %in% names(input_dt)){
    dir.create(plasmid_dir, showWarnings = FALSE)
    for(i in unique(input_dt[,Plasmid])){
      plot_dt_plasmid <- plot_dt[well_id %in% input_dt[Plasmid==i,unlist(well_id)]]
      d <- ggplot2::ggplot(plot_dt_plasmid,ggplot2::aes(time, OD, color = status_id, linetype = well_id)) +
        ggplot2::geom_line(size = 0.5) +
        ggplot2::xlab("Time (h)") +
        ggplot2::ylab("Optical density") +
        ggplot2::coord_cartesian(ylim = range(plot_dt[,OD])) +
        ggplot2::theme_bw()
      if(plot_dt_plasmid[status!="PASS",.N]>20){
        d <- d + ggplot2::scale_colour_manual(values = plot_cols, breaks = names(plot_cols)[!grepl("PASS", names(plot_cols))], guide="none")
      }else{
        d <- d + ggplot2::scale_colour_manual(values = plot_cols, breaks = names(plot_cols)[!grepl("PASS", names(plot_cols))])
      }
      ggplot2::ggsave(file.path(plasmid_dir, paste0(basename(outputPrefix), "_", i, '.pdf')), d, width = plot_width, height = plot_height, useDingbats=FALSE)

    }
  }
}

#Plot growth boxplots
plot_growth_boxplots <- function(
  tecan_data,
  input_dt,
  pass_proportion_threshold = 0.5){

  #Check if required columns present
  if(!'Plasmid' %in% names(result_dt) | !'Well_class' %in% names(result_dt)){
    return()
  }

  #Format plot data
  plot_dt <- copy(input_dt)
  plot_dt[, plot_class := 'sample']
  plot_dt[Well_class!="Sample", plot_class := 'remainder']
  #Sort plasmids by median maxGR and set plasmid pass status
  sort_dt <- plot_dt[,.(median_maxGR = median(maxGR), plasmid_pass = (sum(status=='PASS'))/.N>pass_proportion_threshold),.(Plasmid, plot_class)]
  plot_dt[, Plasmid_factor := factor(Plasmid, levels = sort_dt[order(plot_class, median_maxGR),Plasmid])]
  plot_dt <- merge(plot_dt, sort_dt[,.(Plasmid, plasmid_pass)], by = "Plasmid")
  #Box status
  plot_dt[, box_status := Well_class]
  plot_dt[Well_class=="Sample", box_status := 'FAIL']
  plot_dt[Well_class=="Sample" & plasmid_pass==TRUE, box_status := 'PASS']
  plot_dt[!Well_class %in% c("Positive control", "Negative control", "Sample"), box_status := 'Remainder']
  #Box status colours
  gg_cols = gg_color_hue(3)
  plot_cols <- c("dark grey", gg_cols[1], gg_cols[2], "light grey", gg_cols[3])
  names(plot_cols) <- c("Remainder", "Negative control", "Positive control", "FAIL", "PASS")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Plasmid_factor, maxGR, fill = box_status)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(color = 'black', ggplot2::aes(shape = status)) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Maximum growth rate") +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = plot_cols) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(paste0(outputPrefix, '_boxplots.pdf'), d, width = plot_width, height = plot_height, useDingbats=FALSE)
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
parser <- add_argument(parser, "--method", default = 'heuristic', help = "Maximum growth rate method")
parser <- add_argument(parser, "--parameter", type = "integer", default = 15, help = "h-parameter; number of consecutive time points to evaluate maximum growth rate")
parser <- add_argument(parser, "--wells", default = 'all', help = "Comma-separated list of well ids")
parser <- add_argument(parser, "--deadThreshold", type = "double", default = 0.05, help = "Growth rate threshold for dead variants")
parser <- add_argument(parser, "--lagThreshold", type = "double", default = 48.0, help = "Lag time threshold for problematic variants")
parser <- add_argument(parser, "--outputPrefix", help = "Output path prefix (default: no output file; print results to stdout)")
parser <- add_argument(parser, "--designPath", default = NULL, help = "Path to the plain text file with Well, Plasmid and Well class columns (optional)")
parser <- add_argument(parser, "--plotWidth", type = "integer", default = 8, help = "Plot width in inches (default:8)")
parser <- add_argument(parser, "--plotHeight", type = "integer", default = 8, help = "Plot height in inches (default:8)")
parser <- add_argument(parser, "--ODThreshold", type = "double", default = 0, help = "Minimum optical density required to escape 'deadThreshold' (default:0)")

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
design_path <- args[['designPath']]
plot_width <- args[['plotWidth']]
plot_height <- args[['plotHeight']]
od_threshold <- args[['ODThreshold']]

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
    lag_threshold = lag_threshold,
    design_path = design_path)

###########################
### PLOT GROWTH CURVES
###########################

plot_growth_curves(
  tecan_data = tecan_mat,
  input_dt = result_dt)

###########################
### PLOT GROWTH BOXPLOTS
###########################

plot_growth_boxplots(
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




