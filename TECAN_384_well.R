## Example TECAN experiment analysis
## 02/06/2023
## maximilian.stammnitz@crg.eu

library(readxl)
library(growthrates)
library(beeswarm)
library(scales)

## import 384-well plate design matrix from summary excel sheet
setup <- as.matrix(read_xlsx('TECAN_growths.xlsx', sheet = 1))[2:18,1:25]
colnames(setup) <- setup[1,]
setup <- setup[-1,]
rownames(setup) <- setup[,1]
setup <- setup[,-1]

## import TECAN results (all wells) from summary excel sheet
tecan <- as.matrix(read_xlsx('TECAN_growths.xlsx', sheet = 1))[22:331,1:301]
rownames(tecan) <- tecan[,1]
tecan <- tecan[,-1]

## extract results for three example wells
results <- tecan[c('D2', 'G23', 'N18'),]
colnames(results) <- as.numeric(tecan[1,])/3600 ### convert time to hours
class(results) <- 'numeric' ### convert table to numeric format

## determine each well's max. exponential rates by scanning for the log-linear growth range
## h-parameter: here set to 15 consecutive time points, i.e. growth looking at the best window of 3-4h, feel free to vary/reduce
results.rates <- rep(NA, length = nrow(results))
names(results.rates) <- rownames(results)
for (i in 1:length(results.rates)){
  
  ## summarise well data and run fit
  tmp <- rbind(as.numeric(colnames(results)),
               as.numeric(results[i,]))
  tmp <- tmp[,!is.na(tmp[2,])]
  tmp <- fit_easylinear(time = tmp[1,], y = tmp[2,], h = 15) ### adjust h paramater

  ### you can run plot(tmp) to check if the slope calculation looks good
  ### you can run tmp@par[['lag']] to also calculate the lag phase up to max. ∆OD
  
  ## extract max. ∆OD
  tmp <- tmp@par[['mumax']]
  results.rates[i] <- tmp
  
}

## results.rates
## D2          G23         N18 
## 0.006208704 0.171462067 0.108460533 
