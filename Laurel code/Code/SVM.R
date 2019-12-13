######################################################
######################################################
## Laurel Joannes - laureljoannes@hotmail.com       ##
## November 2018                                    ##
######################################################
######################################################

######################################################
## Label Powerset Method
######################################################
rm(list = setdiff(ls(), lsf.str()))
dir <- getwd()
source(paste0(dir,"/Code/Import.R"))
modelsave <- "SVMLabelPowerSet"
mkdir(paste0(dir, "/Data/", modelsave))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

SVMFunction <- function(trainingSet, testSet, kernelPara, name, gammaPara, costPara){
  name <- paste0("kernel: ", kernelPara, ", gamma: ", gammaPara, ", cost: ", costPara, ", ", name)
  trainADs <- makeListOfADsPerPatientFunction(unique(trainingSet$user_id))
  replaceFrame <- as.data.frame(cbind(trainADs$AD, trainADs$AD))
  colnames(replaceFrame) <- c("prediction", "tempPred")
  replaceFrame$tempPred <- as.character(trainADs$AD)
  testSetADS <- makeListOfADsPerPatientFunction(unique(testSet$user_id))
  
  trainingSet <- merge(trainADs, trainingSet, by = "user_id")
  trainingSet <- trainingSet[order(trainingSet$user_id),]
  trainingSet$AD <- as.character(trainingSet$AD)
  
  set.seed(100)
  model_svm <- svm(AD~. , data  = trainingSet[,-1], scale = FALSE, type = "C-classification", gamma = gammaPara, 
                   kernel = kernelPara, cost = costPara)
  saveRDS(model_svm, paste0(dir, "/Data/", modelsave, "/Models/",name,".rda"))
  predicted <- predict(model_svm, testSet[,-1])
  print(paste0(name, " predicts ", length(unique(predicted)), " unique diseases" ))
  
  predicted <- cbind(predicted, testSetADS)
  colnames(predicted) <- c("tempPred", "user_id", "real")
  predicted <- predicted[,c(2,1,3)]
  
  prediction <- merge(predicted, replaceFrame, by = "tempPred")
  prediction <- prediction[!duplicated(prediction), ]
  prediction <- prediction[,c(2,4,3)]
  
  sum <- 0
  disease <- c()
  for(i in 1:nrow(prediction)){
    sum <- sum + length(prediction$prediction[i][[1]])
    disease <- as.vector(sort(unique(c(disease, unlist(prediction$prediction[i][[1]])))))
  }
  print(paste0(name," predicts ", length(disease), "  unique diseases"))
  print(paste0(name, " on average predicts " , sum / nrow(prediction), " diseases per patient"))
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name)
  performances <- calculatePerformanceFunction(prediction, name, "no")
  outputs <- list(performances, scores)
  return(outputs)
}

# 5-Fold Cross Validation
set.seed(5)
flds <- createFolds(unique(train$user_id), k = 5, list = TRUE, returnTrain = FALSE)
costPara <- c(0.1, 1, 10, 100, 500, 1000, 10000)
performances <- scores <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(costPara)){
    sampleTestIDs <- train[flds[[j]],]$user_id
    OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
    tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
    tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
    outputs <- SVMFunction(tempTrainset, tempTestset, "linear", "5-fold", 0, costPara[i])
    performances <- rbind(performances, as.data.frame(outputs[1]))
    scores <- rbind(scores, as.data.frame(outputs[2]))
  }
}


gammaPara <- c(2^(-(length(train)-1)), 2^(-(length(train)-1)/2), 2^-100, 2^-10, 1/(length(train)-1), 1/((length(train)-1)/2), 1)
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(costPara)){
    for(k in 1:length(gammaPara)){
      sampleTestIDs <- train[flds[[j]],]$user_id
      OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
      tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
      tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
      outputs <- SVMFunction(tempTrainset, tempTestset, "radial", "5-fold", gammaPara[k], costPara[i])
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

endScores <- endPerformances <- c()
temp1 <- dplyr::filter(scores, (grepl("linear", scores$Name)))
temp2 <- dplyr::filter(performances, (grepl("linear", performances$Name)))
for(i in 1:length(costPara)){
  temp3 <- dplyr::filter(temp1, (grepl(paste0("cost: ", costPara[i],","), temp1$Name)))
  temp4 <- dplyr::filter(temp2, (grepl(paste0("cost: ", costPara[i],","), temp2$Name)))
  endScores <- rbind(endScores, c(temp3$Name[1], round(colSums(temp3[,-1])/nrow(temp3),4)))
  endPerformances <- rbind(endPerformances, c(temp4$Name[1], round(colSums(temp4[,-1])/nrow(temp4),4)))
}

temp1 <- dplyr::filter(scores, (grepl("radial", scores$Name)))
temp2 <- dplyr::filter(performances, (grepl("radial", performances$Name)))
for(i in 1:length(costPara)){
  temp3 <- dplyr::filter(temp1, (grepl(paste0("cost: ", costPara[i],","), temp1$Name)))
  temp4 <- dplyr::filter(temp2, (grepl(paste0("cost: ", costPara[i],","), temp2$Name)))
  for(j in 1:length(gammaPara)){
    temp5 <- dplyr::filter(temp3, (grepl(paste0("gamma: ", gammaPara[j], ","), temp3$Name)))
    temp6 <- dplyr::filter(temp4, (grepl(paste0("gamma: ", gammaPara[j],","), temp4$Name)))
    endScores <- rbind(endScores, c(temp5$Name[1], round(colSums(temp5[,-1])/nrow(temp5),4)))
    endPerformances <- rbind(endPerformances, c(temp6$Name[1], round(colSums(temp6[,-1])/nrow(temp6),4)))
  }
}

write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.txt"), sep=" & ")
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.txt"), sep=" & ")

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/end_5Fold.txt"), sep=" & ", row.names = FALSE)

# Validation
endScores <- endPerformances <- c()

gammaPara <- c(2^-10, 1/(length(train)-1), 1/((length(train)-1)/2))
for(i in 1: length(gammaPara)){
  outputs <- SVMFunction(train, validation, "radial", "Validation", gammaPara[i], 100)
  endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
  endScores <- rbind(endScores, as.data.frame(outputs[2]))
}

endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))


write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Validation.txt"), sep=" & ")
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Validation.txt"), sep=" & ")


# Test 
endScores2 <- endPerformances2 <- c()

trainValidation <- rbind(train, validation)
trainValidation <- trainValidation[order(trainValidation$user_id),]

gammaPara <- c(2^-10, 1/(length(train)-1), 1/((length(train)-1)/2))
for(i in 1:length(gammaPara)){
  outputs <- SVMFunction(trainValidation, test, "radial", "test", gammaPara[i], 100)
  endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
  endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
}


endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))



write.csv(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.txt"), sep=" & ")
write.csv(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.txt"), sep=" & ")

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_ValTest.txt"), sep=" & ")
write.table(endPerformances3, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_ValTest.txt"), sep=" & ")


rm(list = setdiff(ls(), lsf.str()))
######################################################
## Binary Relevance Method
######################################################
dir <- getwd()
source(paste0(dir,"/Code/Import.R"))
modelsave <- "SVMBinary"
mkdir(paste0(dir, "/Data/", modelsave, "/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

BinarySVMFunction <- function(trainingSet, testSet, testUsers, kernelPara, name, gammaPara, costPara, threshold, j){
  name2 <- paste0(" Kernel = ", kernelPara, " gamma = ", gammaPara , " Cost = ", costPara, " threshold = ", threshold, " ", name)
  if(j == 0 ){
    name3 <- paste0(dir,"/Data/", modelsave, "/Models/",name2, ".rda")
  }else{
    map <- paste0("5fold",j)
    name3 <- paste0(dir,"/Data/", modelsave, "/Models/",map,"/",name2, ".rda")
  }
  
  testADs <- makeListOfADsPerPatientFunction(testUsers)
  
  output.model <- ebr(trainingSet, "SVM", gamma = gammaPara, type = "C-classification", 
                      kernel = kernelPara, cost = costPara, scale=FALSE, cores = 4, seed = 2018, attr.space = 1, subsample = 1, m = 10)
  prediction  <- predict(output.model, testSet, cores = 4)
  
  saveRDS(output.model, file = name3)
  
  prediction <- as.data.frame(as.matrix(prediction))
  prediction$prediction <- NA
  
  for(i in 1:nrow(prediction)){
    temp <- c()
    for(j in 1:(ncol(prediction)-2)){
      if(prediction[i,j] > threshold){
        temp <- c(temp, j)
      }
    }
    if(length(temp) == 0){
      temp <- 0
    }
    prediction$prediction[i] <- list(temp)
  }
  
  prediction <- cbind(prediction, testADs)
  prediction <- prediction[,c(19,18,20)]
  colnames(prediction) <- c("user_id", "prediction", "real")
  
  sum <- 0
  disease <- c()
  for(i in 1:nrow(prediction)){
    sum <- sum + length(prediction$prediction[i][[1]])
    disease <- as.vector(sort(unique(c(disease, unlist(prediction$prediction[i][[1]])))))
  }
  print(paste0(name2," predicts ", length(disease), "  unique diseases"))
  print(paste0(name2, " on average predicts " , sum / nrow(prediction), " diseases per patient"))
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name2)
  performances <- calculatePerformanceFunction(prediction, name2, "no")
  
  
  outputs <- list(performances, scores)
  return(outputs)
}


# 5-Fold Cross Validation
set.seed(5)
flds <- createFolds(unique(trainMultiLabel$user_id), k = 5, list = TRUE, returnTrain = FALSE)
costPara <- c(0.1, 1, 10, 100, 500, 1000)
thresholdPara <- c(0.3,0.4,0.5)

performances <- scores <- c()

for(j in 1:5){
  print(paste0(j, " of 5, linear"))
  for(i in 1:length(costPara)){
    for(l in 1:length(thresholdPara)){
      sampleTestIDs <- train[flds[[j]],]$user_id
      OtherTrainIDs <- setdiff(unique(trainMultiLabel$user_id), sampleTestIDs)
      tempTempTestset <- dplyr::filter(trainMultiLabel, user_id %in% sampleTestIDs)
      tempTempTrainset <- dplyr::filter(trainMultiLabel, user_id %in% OtherTrainIDs)
      
      cols <- sapply(tempTempTrainset, is.logical)
      tempTempTrainset[,cols] <- lapply(tempTempTrainset[,cols], as.numeric)
      cols <- sapply(tempTempTestset, is.logical)
      tempTempTestset[,cols] <- lapply(tempTempTestset[,cols], as.numeric)
      
      tempTrainset <- mldr_from_dataframe(tempTempTrainset[,-1], labelIndices =  c(477:493))
      tempTestset <- mldr_from_dataframe(tempTempTestset[,-1], labelIndices =  c(477:493))
      
      outputs <- BinarySVMFunction(tempTrainset, tempTestset, sampleTestIDs, "linear", "5-fold", 0, costPara[i], thresholdPara[l], j)
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

scores <- readRDS(paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
performances <- readRDS(paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))


gammaPara <- c(2^-10, 1/(length(trainMultiLabel)-1), 1/((length(trainMultiLabel)-1)/2), 1)
for(j in 1:5){
  print(paste0(j, " of 5, radial"))
  for(i in 1:length(costPara)){
    for(l in 1:length(thresholdPara)){
      for(k in 1:length(gammaPara)){
        sampleTestIDs <- train[flds[[j]],]$user_id
        OtherTrainIDs <- setdiff(unique(trainMultiLabel$user_id), sampleTestIDs)
        tempTempTestset <- dplyr::filter(trainMultiLabel, user_id %in% sampleTestIDs)
        tempTempTrainset <- dplyr::filter(trainMultiLabel, user_id %in% OtherTrainIDs)
        
        cols <- sapply(tempTempTrainset, is.logical)
        tempTempTrainset[,cols] <- lapply(tempTempTrainset[,cols], as.numeric)
        cols <- sapply(tempTempTestset, is.logical)
        tempTempTestset[,cols] <- lapply(tempTempTestset[,cols], as.numeric)
        
        tempTrainset <- mldr_from_dataframe(tempTempTrainset[,-1], labelIndices =  c(477:493))
        tempTestset <- mldr_from_dataframe(tempTempTestset[,-1], labelIndices =  c(477:493))
        
        outputs <- BinarySVMFunction(tempTrainset, tempTestset, sampleTestIDs, "radial", "5-fold", gammaPara[k], costPara[i], thresholdPara[l], j)
        performances <- rbind(performances, as.data.frame(outputs[1]))
        scores <- rbind(scores, as.data.frame(outputs[2]))
      }
    }
  }
}




saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

endScores <- endPerformances <- c()
temp1 <- dplyr::filter(scores, (grepl("linear", scores$Name)))
temp2 <- dplyr::filter(performances, (grepl("linear", performances$Name)))
for(i in 1:length(costPara)){
  temp3 <- dplyr::filter(temp1, (grepl(paste0("Cost = ", costPara[i], " t"), temp1$Name)))
  temp4 <- dplyr::filter(temp2, (grepl(paste0("Cost = ", costPara[i], " t"), temp2$Name)))
  for(j in 1:length(thresholdPara)){
    temp5 <- dplyr::filter(temp3, (grepl(paste0("threshold = ", thresholdPara[j], " 5"), temp3$Name)))
    temp6 <- dplyr::filter(temp4, (grepl(paste0("threshold = ", thresholdPara[j], " 5"), temp4$Name)))
    endScores <- rbind(endScores, c(temp5$Name[1], round(colSums(temp5[,-1])/nrow(temp5),4)))
    endPerformances <- rbind(endPerformances, c(temp6$Name[1], round(colSums(temp6[,-1])/nrow(temp6),4)))
  }
}

temp1 <- dplyr::filter(scores, (grepl("radial", scores$Name)))
temp2 <- dplyr::filter(performances, (grepl("radial", performances$Name)))
for(i in 1:length(costPara)){
  temp3 <- dplyr::filter(temp1, (grepl(paste0("Cost = ", costPara[i], " t"), temp1$Name)))
  temp4 <- dplyr::filter(temp2, (grepl(paste0("Cost = ", costPara[i], " t"), temp2$Name)))
  for(j in 1:length(gammaPara)){
    temp5 <- dplyr::filter(temp3, (grepl(paste0("gamma = ", gammaPara[j], " C"), temp3$Name)))
    temp6 <- dplyr::filter(temp4, (grepl(paste0("gamma = ", gammaPara[j], " C"), temp4$Name)))
    for(k in 1:length(thresholdPara)){
      temp7 <- dplyr::filter(temp5, (grepl(paste0("threshold = ", thresholdPara[k], " 5"), temp5$Name)))
      temp8 <- dplyr::filter(temp6, (grepl(paste0("threshold = ", thresholdPara[k], " 5"), temp6$Name)))
      endScores <- rbind(endScores, c(temp7$Name[1], round(colSums(temp7[,-1])/nrow(temp7),4)))
      endPerformances <- rbind(endPerformances, c(temp8$Name[1], round(colSums(temp8[,-1])/nrow(temp8),4)))
    }
  }
}

write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.txt"), sep=" & ", row.names=FALSE)

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/end_5Fold.txt"), sep=" & ", row.names=FALSE)


#### Validation

cols <- sapply(trainMultiLabel, is.logical)
trainMultiLabel[,cols] <- lapply(trainMultiLabel[,cols], as.numeric)
cols <- sapply(validationMultiLabel, is.logical)
validationMultiLabel[,cols] <- lapply(validationMultiLabel[,cols], as.numeric)

trainBinary <- mldr_from_dataframe(trainMultiLabel[,-1], labelIndices =  c(477:493))
validationBinary <- mldr_from_dataframe(validationMultiLabel[,-1], labelIndices =  c(477:493))

endScores <- endPerformances <- c()

outputs <- BinarySVMFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), "radial", "Validation", 1/(length(trainMultiLabel)-1), 10, 0.3, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- BinarySVMFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), "radial", "Validation", 1/((length(trainMultiLabel)-1)/2), 10, 0.3, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- BinarySVMFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), "radial", "Validation", 2^(-10), 100, 0.3, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))

endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))


write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Validation.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Validation.txt"), sep=" & ", row.names=FALSE)

# Test 
endScores2 <- endPerformances2 <- c()
trainValidation <- rbind(trainMultiLabel, validationMultiLabel)
trainValidation <- trainValidation[order(trainValidation$user_id),]

cols <- sapply(trainValidation, is.logical)
trainValidation[,cols] <- lapply(trainValidation[,cols], as.numeric)
trainValidationBinary <- mldr_from_dataframe(trainValidation[,-1], labelIndices = c(477:493))

cols <- sapply(testMultiLabel, is.logical)
testMultiLabel[,cols] <- lapply(testMultiLabel[,cols], as.numeric)
testBinary <- mldr_from_dataframe(testMultiLabel[,-1], labelIndices = c(477:493))

endScores2 <- endPerformances2 <- c()

outputs <- BinarySVMFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id), "radial", "Test", 1/(length(trainMultiLabel)-1), 10, 0.3, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- BinarySVMFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id), "radial", "Test",1/((length(trainMultiLabel)-1)/2), 10, 0.3, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- BinarySVMFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id),  "radial", "Test", 2^(-10), 100, 0.3, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))

endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))

write.csv(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.txt"), sep=" & ", row.names=FALSE)

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_ValTest.txt"), sep=" & ", row.names=FALSE)
write.table(endPerformances3, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_ValTest.txt"), sep=" & ", row.names=FALSE)



rm(list = setdiff(ls(), lsf.str()))
######################################################
## Algorithm Adaptation Method
######################################################
dir <- getwd()
source(paste0(dir,"/Code/Import.R"))
modelsave <- "SVMMultilabel"
mkdir(paste0(dir, "/Data/", modelsave, "/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))


multiLabelSVMFunction <- function(trainingSet, testSet, kernelPara, name, gammaPara, costPara, thresholdPara, j){
  if(kernelPara == "Linear"){
    kernelPara2 <- 0
  }else if(kernelPara == "Radial"){
    kernelPara2 <- 2
  }
  
  name2 <- paste0("kernel: ", kernelPara, ", gamma: ", gammaPara, ", cost: ", costPara, ", ", "threshold: ", thresholdPara, ", ", name)[1]
  
  if(j == 0 ){
    name3 <- paste0(dir,"/Data/", modelsave, "/Models/",name2, ".rda")
  }else{
    map <- paste0("5fold",j)
    name3 <- paste0(dir, "/Data/", modelsave, "/Models/",map,"/",name2, ".rda")
  }
  
  trainAD <- dplyr::filter(usersWithAD, user_id %in% unique(trainingSet$user_id))
  newTrain <- c()
  temp <- trainingSet
  for(i in 2:9){
    temp$AD <- as.vector(unlist(trainAD[,i]))
    newTrain <- rbind(newTrain, temp, make.row.names = FALSE)
  }
  newTrain <- newTrain[order(newTrain$user_id),]
  
  newTrain <- na.omit(newTrain)
  
  svmOptions <- paste0("-j ", costPara, " -t ", kernelPara2, " -g ", gammaPara)
  set.seed(100)
  model_svm <- svmlight(AD ~. ,data = newTrain[,-1], svm.options = svmOptions ,pathsvm = dir, type = "C")
  saveRDS(model_svm, file = name3)
  
  testADs <- makeListOfADsPerPatientFunction(unique(testSet$user_id))
  prediction <- as.data.frame(predict(model_svm, testSet[,-1], type = "C")$posterior)
  prediction <- prediction[, (as.character(c(1:16,0)))]
  prediction$prediction <- NA
  for(i in 1:nrow(prediction)){
    temp <- c()
    for(j in 1:(ncol(prediction)-2)){
      if(prediction[i,j] >= thresholdPara){
        temp <- c(temp, j)
      }
    }
    if(length(temp) == 0){
      temp <- 0
    }
    prediction$prediction[i] <- list(temp)
  }
  prediction <- cbind(prediction, testADs)
  prediction <- prediction[,c(19,18,20)]
  colnames(prediction) <- c("user_id", "prediction", "real")
  
  sum <- 0
  disease <- c()
  for(i in 1:nrow(prediction)){
    sum <- sum + length(prediction$prediction[i][[1]])
    disease <- as.vector(sort(unique(c(disease, unlist(prediction$prediction[i][[1]])))))
  }
  print(paste0(name2," predicts ", length(disease), "  unique diseases"))
  print(paste0(name2, " on average predicts " , sum / nrow(prediction), " diseases per patient"))
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name2)
  performances <- calculatePerformanceFunction(prediction, name2, "no")
  
  outputs <- list(performances, scores)
  return(outputs)
}


set.seed(5)
flds <- createFolds(unique(train$user_id), k = 5, list = TRUE, returnTrain = FALSE)
costPara <- c(0.1, 1, 10, 25, 50, 100, 500, 1000, 10000)
thresholdPara <- c(0.1, 0.2, 0.3)
performances <- scores <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(costPara)){
    for(k in 1:length(thresholdPara)){
      sampleTestIDs <- train[flds[[j]],]$user_id
      OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
      tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
      tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
      outputs <- multiLabelSVMFunction(tempTrainset, tempTestset, "Linear", "5-fold", 0, costPara[i], thresholdPara[k], j)
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold1.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold1.rda"))

gammaPara <- c(2^-10, 1/(length(train)-1), 1/((length(train)-1)/2), 1)
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(costPara)){
    for(k in 1:length(gammaPara)){
      for(l in 1:length(thresholdPara)){
        sampleTestIDs <- train[flds[[j]],]$user_id
        OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
        tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
        tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
        outputs <- multiLabelSVMFunction(tempTrainset, tempTestset, "Radial", "5-fold", gammaPara[k], costPara[i], thresholdPara[l], j)
        performances <- rbind(performances, as.data.frame(outputs[1]))
        scores <- rbind(scores, as.data.frame(outputs[2]))
      }
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

endScores <- endPerformances <- c()
temp1 <- dplyr::filter(scores, (grepl("Linear", scores$Name)))
temp2 <- dplyr::filter(performances, (grepl("Linear", performances$Name)))
for(i in 1:length(costPara)){
  temp3 <- dplyr::filter(temp1, (grepl(paste0("cost: ", paste0(costPara[i],",")), temp1$Name)))
  temp4 <- dplyr::filter(temp2, (grepl(paste0("cost: ", paste0(costPara[i],",")), temp2$Name)))
  for(j in 1:length(thresholdPara)){
    temp5 <- dplyr::filter(temp3, grepl(paste0("threshold: ", thresholdPara[j]), temp3$Name))
    temp6 <- dplyr::filter(temp4, grepl(paste0("threshold: ", thresholdPara[j]), temp4$Name))
    endScores <- rbind(endScores, c(temp5$Name[1], round(colSums(temp5[,-1])/nrow(temp5),4)))
    endPerformances <- rbind(endPerformances, c(temp6$Name[1], round(colSums(temp6[,-1])/nrow(temp6),4)))
  }
  
}

temp1 <- dplyr::filter(scores, (grepl("Radial", scores$Name)))
temp2 <- dplyr::filter(performances, (grepl("Radial", performances$Name)))
for(i in 1:length(costPara)){
  temp3 <- dplyr::filter(temp1, (grepl(paste0("cost: ", paste0(costPara[i],",")), temp1$Name)))
  temp4 <- dplyr::filter(temp2, (grepl(paste0("cost: ", paste0(costPara[i],",")), temp2$Name)))
  for(j in 1:length(gammaPara)){
    temp5 <- dplyr::filter(temp3, (grepl(paste0("gamma: ", gammaPara[j]), temp3$Name)))
    temp6 <- dplyr::filter(temp4, (grepl(paste0("gamma: ", gammaPara[j]), temp4$Name)))
    for(k in 1:length(thresholdPara)){
      temp7 <- dplyr::filter(temp5, grepl(paste0("threshold: ", thresholdPara[k]), temp5$Name))
      temp8 <- dplyr::filter(temp6, grepl(paste0("threshold: ", thresholdPara[k]), temp6$Name))
      endScores <- rbind(endScores, c(temp7$Name[1], round(colSums(temp7[,-1])/nrow(temp7),4)))
      endPerformances <- rbind(endPerformances, c(temp8$Name[1], round(colSums(temp8[,-1])/nrow(temp8),4)))
    }
  }
}

write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.txt"), sep=" & ")
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.txt"), sep=" & ")

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/end_5Fold.txt"), sep=" & ", row.names = FALSE)


### Validation
endScores <- endPerformances <- c()

outputs <- multiLabelSVMFunction(train, validation, "Radial", "Validation", 2^(-10), 10, 0.2, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- multiLabelSVMFunction(train, validation, "Radial", "Validation", 1/((length(train)-1)/2), 10, 0.2, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- multiLabelSVMFunction(train, validation, "Radial", "Validation", 1/((length(train)-1)/2), 10, 0.1, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))

endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))

write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Validation.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Validation.txt"), sep=" & ", row.names=FALSE)


### Test
endScores2 <- endPerformances2 <- c()
trainValidation <- rbind(train, validation)
trainValidation <- trainValidation[order(trainValidation$user_id),]

outputs <- multiLabelSVMFunction(trainValidation, test, "Radial", "Test", 2^(-10), 10, 0.2, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- multiLabelSVMFunction(trainValidation, test, "Radial", "Test", 1/((length(train)-1)/2), 10, 0.2, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- multiLabelSVMFunction(trainValidation, test, "Radial", "Test", 1/((length(train)-1)/2), 10, 0.1, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))

endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))


write.csv(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.txt"), sep=" & ", row.names=FALSE)

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_ValTest.txt"), sep=" & ", row.names=FALSE)
write.table(endPerformances3, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_ValTest.txt"), sep=" & ", row.names=FALSE)






rm(list = setdiff(ls(), lsf.str()))

