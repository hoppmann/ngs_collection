
options(warn = -1)

## retrieve commandline options
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  print(paste("This script need exactly 2 input parameter. First fileIn and then FileOut\n ", length(args), " parameters were provided"))
  stop()
}

fileIn <- args[1]
fileOut <- args[2]


if (!file.exists(fileIn)) {
  print(paste("ERROR: File ", fileIn, " not found!"))
  stop()
}


## give out filename currently running on
print (fileIn)


## read in results table data
data <- read.table(fileIn, sep = "\t", header = T)

## create table for collecting new annotated data
outData <- data.frame()


## go throu file line wise and add modifications

for (i in 1:nrow(data)) {
# for (a in 1:1) {
  
  ###########################################
  ######## get number of predictions ########
  ###########################################

  ## get row
  row <- data[i,]
  
  ## reset caller counter
  damCaller <- 0
  allCaller <- 0
  
  
  #######################
  ######## siftPrediction
  siftscore = as.numeric(as.character(row$sift_score))
  if(! is.na(siftscore)) {
    
    allCaller = allCaller + 1
    
    if (siftscore < 0.05){
      damCaller = damCaller + 1
    }
    
  }
  
  ################
  ####### polyphen
  polyscore = as.numeric(as.character(row$polyphen_score))
  
  if(! is.na(polyscore)) {
    
    allCaller = allCaller + 1
    
    if (polyscore > 0.446){
      damCaller = damCaller + 1
    }
    
  }
  
  
  
  
  ############
  ###### mealr
  metaLRscore = as.numeric(as.character(row$vep_metalr_score))
  
  if(! is.na(metaLRscore)) {
    
    allCaller = allCaller + 1
    
    if (metaLRscore >= 0.5){
      damCaller = damCaller + 1
    }
    
  }
  
  
  
  #############
  ###### meaSVM
  metaSVMscore = as.numeric(as.character(row$vep_metasvm_score))
  
  if(! is.na(metaSVMscore)) {
    
    allCaller = allCaller + 1
    
    if (metaSVMscore >= 0.5){
      damCaller = damCaller + 1
    }
    
  }
  
  
  
  ##################
  ###### REVEL score
  revelScore = as.numeric(as.character(row$vep_revel_score))
  revelPred = ""
  
  if(! is.na(revelScore)) {
    
    allCaller = allCaller + 1
        
    if (revelScore > 0.5){
      damCaller = damCaller + 1
      revelPred = "D"
    } else {
      revelPred = "T"
    }
    
  }
  
  
  
  #################
  ###### CADD score
  caddScore = as.numeric(as.character(row$cadd_scaled))
  caddPred = ""
  if(! is.na(caddScore)) {
    
    allCaller = allCaller + 1
    
    if (caddScore > 15){
      damCaller = damCaller + 1
      caddPred = "D"
    } else {
      caddPred = "T"
    }
    
  }
  
  
  
  
  
  ##############
  ####### fathmm
  allFat = as.numeric(unlist(strsplit(as.character(row$vep_fathmm_score), "&")))
  if (length(allFat) < 1) {
    fathmmScore = NA
  } else {
    fathmmScore = max(allFat)
    fathmmScore = as.numeric(as.character(row$vep_fathmm_score))
  }
  if(! is.na(fathmmScore) ) {
    
    allCaller = allCaller + 1
    
    if (fathmmScore <= -2.5){
      damCaller = damCaller + 1
    }
    
  }
  
  #########################
  ######## write as decimal 
  
  percentDamagingPred <- damCaller / allCaller
  
  ##########################
  #### add additional column

  row["revel_score_pred"] <- revelPred
  row["cadd_pred"] <- caddPred  
  row["damagingPredictions"] <- damCaller
  row["totalPredictions"] <- allCaller
  row["percentDamagingPredictions"] <- percentDamagingPred

  
  ## clean variables
  rm(damCaller, allCaller, percentDamagingPred, fathmmScore, metaLRscore, metaSVMscore,
     allFat, caddScore, polyscore, siftscore, revelScore)
  
  
  ##########################################
  ######## activity of splice sites ########
  ##########################################

  ## init/reset variable
  allSSPredValues = c()

  ###############
  ####  EntScore

  ## calculate max differnece
  varMaxEntScore <- as.numeric(unlist(strsplit(as.character(row$alamut_varMaxEntScore), ",")))
  wtMaxEntScore <- as.numeric(unlist(strsplit(as.character(row$alamut_wtMaxEntScore), ",")))
  percent <- varMaxEntScore / wtMaxEntScore

  allSSPredValues <- append(allSSPredValues, min(percent))
  # minPercentEntScore <- min(percent)
  rm(varMaxEntScore, wtMaxEntScore, percent)


  ############
  #### NNScore

  varNNScore <- as.numeric(unlist(strsplit(as.character(row$alamut_varNNSScore), ",")))
  wtNNScore <- as.numeric(unlist(strsplit(as.character(row$alamut_wtNNSScore), ",")))
  percent <- varNNScore / wtNNScore

  allSSPredValues <- append(allSSPredValues, min(percent))
  # minPercentNNScore <- min(percent)
  rm(varNNScore, wtNNScore, percent)


  ############
  #### GSScore

  varGSScore <- as.numeric(unlist(strsplit(as.character(row$alamut_varGSScore), ",")))
  wtGSScore <- as.numeric(unlist(strsplit(as.character(row$alamut_wtGSScore), ",")))
  percent <- varGSScore / wtGSScore

  allSSPredValues <- append(allSSPredValues, min(percent))
  # minPercentGSScore <- min(percent)
  rm(varGSScore, wtGSScore, percent)



  #############
  #### SSFSCore

  varSSFScore <- as.numeric(unlist(strsplit(as.character(row$alamut_varSSFScore), ",")))
  wtSSFScore <- as.numeric(unlist(strsplit(as.character(row$alamut_wtSSFScore), ",")))
  percent <- varSSFScore / wtSSFScore

  allSSPredValues <- append(allSSPredValues, min(percent))
  # minPercentSSFScore <- min(percent)
  rm(varSSFScore, wtSSFScore, percent)


  #####################################
  #### calculate majority for splicing

  # ## add values to vector to work with vector later on
  # allSSPredValues <- c(minPercentEntScore, minPercentGSScore, minPercentNNScore, minPercentSSFScore)

  ## get number of predictions not na
  allSSPred <- sum(!is.na(allSSPredValues))
  highSSPred <- 0
  medSSPred <- 0


  ## check for each prediction if loss is >45% -> high catagory
  ## or (<45% AND >15%) -> med catagory
  ## check if value is na if so no catagory
  for (j in 1:length(allSSPredValues)){
    if (!is.na(allSSPredValues[j]) && allSSPredValues[j] < 0.55){
      highSSPred = highSSPred +1
    } else if (!is.na(allSSPredValues[j]) && allSSPredValues[j] >= 0.45 && allSSPredValues[j] < 0.85){
    } else if (!is.na(allSSPredValues[j]) && allSSPredValues[j] >= 0.45 && allSSPredValues[j] < 0.85){
      medSSPred = medSSPred +1
    }
  }


  #### add additional column for SS predictions
  row["SSPred15"] <- medSSPred
  row["SSPred45"] <- highSSPred
  row["allSSPred"] <- allSSPred
  row["SSPred15Dec"] <- medSSPred/allSSPred
  row["SSPred45Dec"] <- highSSPred/allSSPred


  ####################
  #### clean variables
  rm(j, allSSPred, allSSPredValues, highSSPred, medSSPred, minPercentEntScore, minPercentGSScore, minPercentNNScore, minPercentSSFScore)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### add annotated row to new Dataframe
  outData <- rbind(outData, row)
   
} # end of iteration over each line


######## remove unwanted columns
# outData <- outData[,-c("alamut_wtMaxEntScore", "alamut_varMaxEntScore", "alamut_wtNNSScore", 
#                  "alamut_varNNSScore", "alamut_wtGSScore", "alamut_varGSScore", "alamut_wtSSFScore",
#                  "alamut_varSSFScore", "sift_score", "polyphen_score", "vep_fathmm_score",
#                  "vep_metalr_score", "vep_metasvm_score", "vep_revel_score")]


## splice site prediction
outData$alamut_wtGSScore <- NULL
outData$alamut_varGSScore <- NULL
outData$alamut_wtMaxEntScore <- NULL
outData$alamut_varMaxEntScore <- NULL
outData$alamut_wtNNSScore <- NULL
outData$alamut_varNNSScore <- NULL
outData$alamut_wtSSFScore <- NULL
outData$alamut_varSSFScore <- NULL

## prediction tools
outData$sift_score <- NULL
outData$polyphen_score <- NULL
outData$vep_fathmm_score <- NULL
outData$vep_metalr_score <- NULL
outData$vep_metasvm_score <- NULL
outData$vep_revel_score <- NULL
outData$cadd_scaled <- NULL

# index <- which(colnames(outData)=="")
# outData <- outData[,]


#### write out new file
write.table(outData, fileOut, col.names = T, row.names = F, quote = F, sep = "\t")























