##########################################
#### recovery data from TECAN MACHINE ####
#### Aina M August 2020               ####
##########################################

#####
# This script will generate an excel table with all the OD values for each cycle and well giving as input the log files generated by the TECAN program

##### FILES and FORMAT
# This script needs the log temporary files, which are stored locally in the TECAN computer:
#   On Windows XP computers, the log files (txt files) are stored in: 
#         \Documents\Tecan\LogFiles\icontrol\1.10\<Instrument_Serial_Number>\.
#   Their name will contain the word "measurement" on it. Here is an example: 
#         "RNwork_Measurement-log.txt.2020-08-25.2"
#   Your data may be distributed among different log files, identify all of them and put them together in the same folder.
#   Name all the files in a similar way. For example, adding 001, 002, 003 at the beginning of the file names and in the order the data is stored, being 001 the file containing the information of the first cycles, and so on.
#   You can see which information is stored in the files by identifying the rows like the following one:
#         15:13:11,405 [Measurement Thread  ] INFO  - SCAN for Well=E11 Row=4/Col=10  MeasNr=1 Cycle=7=11:300   [SystemTicker=1766229126]
#         Here we have the time of the measurement, the well to which it corresponds (E11), and the cycle (11)
#         Some rows below you will find the OD measurement: 
#         15:13:11,420 [Measurement Thread  ] INFO  - OD Value = 0,0985153565730153


#### TECAN 384 well plate
# this might not work for 384 well plate tecans since the storage exceeds the maximum allowed for the machine, thus it rewrites the information every ~12 hours and the previous data is lost in the log files.
# this will likely produce gaps in the recovered data
####

#### CYCLES 
# If you have some cycle numbers repeated in the working files (e.g you have restarted the tecan machine after some cycles)
#   the information from the older cycles won't be stored in the new Excel file. The new cycles with the same number will overwrite the previous data.
#   if you want to store the data from those cycles as well, generate a separate file with the first TECAN running cycles and run this script separately
####


####### The following parameters may need to be changed
        # Select the log files that contain the information about the OD measurements. In this case, all the files starting by "00" in this folder
files_path <- "/Users/Aina/Desktop/tecan_data_recovery/00*"
        # Specify the number of cycles done. This will be the number of columns in the new table that will be created
num_cycles <- 89 
        # Specify the name of the wells for which you have measurements. In this case, those are all the wells except the border ones (e.g 60 wells out of the total 96 wells)
wells=c(paste0("B", 2:11), paste0("C", 2:11), paste0("D", 2:11), paste0("E", 2:11), paste0("F", 2:11), paste0("G", 2:11))
        # Specify the name of the output excel file.It will be stored in the same folder where the script is located
output_file='recovered_TECAN_data.xlsx'
#######


#load libraries
if (!require("readtext")) install.packages("readtext")
if (!require("stringr")) install.packages("stringr")
if (!require("data.table")) install.packages("data.table")
if (!require("openxlsx")) install.packages("openxlsx")

library(readtext)
library(stringr)
library(data.table)
library(openxlsx)

#read the files
files <- readtext(paste0(files_path))
all_text<-files[2]

#create an empty data table
dt<- data.table( well_id = wells)
for(i in 1:num_cycles){
  dt[, paste(i):=0]
}

#start reading all files and storing the OD data
well<-""
for(i in 1:nrow(all_text)){
  txt<-all_text[i,]
  lines<- strsplit(txt, "\n")
  for (l in 1:length(lines[[1]])){
    if(str_detect(lines[[1]][l], "SCAN")){
      line_parts<- strsplit(lines[[1]][l], " ")
      well=str_sub(line_parts[[1]][11], 6, -1)
      cycle<-as.numeric(strsplit(strsplit(line_parts[[1]][15], "=")[[1]][3], ":")[[1]][1])
    }
    else if(str_detect(lines[[1]][l], "OD Value")){
      if (well!=""){
        OD=as.numeric(sub(",", ".", str_sub(lines[[1]][l], 56, 72), fixed = TRUE))
        dt[ well_id==well, cycle+1] <- OD
        well<-""
      }
    }
    
  }
}

#save the data into an excel file
write.xlsx(dt, output_file)

