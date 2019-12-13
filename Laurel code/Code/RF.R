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
modelsave <- "RandomForestLabelPowerset"
mkdir(paste0(dir, "/Data/", modelsave,"/Models/FINAL"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))


randomForestFunction<- function(trainingSet, testset, name, ntree, mtry){
  name2 <- paste0(" ntree = ", ntree, ", mtry = ", mtry, name)
  name3 <- paste0(dir, "/Data/", modelsave, "/Models/",name2, ".rda")
  name4 <- paste0(dir, "/Data/", modelsave, "/Models/FINAL ",name2, ".rda")

  trainADs <- makeListOfADsPerPatientFunction(unique(trainingSet$user_id))
  trainingSet <- merge(trainADs, trainingSet, by = "user_id")
  trainingSet <- trainingSet[order(trainingSet$user_id),]
  trainingSet$AD <- as.factor(as.character(trainingSet$AD))
  set.seed(50)
  output.forest <- randomForest(x = trainingSet[,-which(names(trainingSet) %in% c('user_id', 'AD'))], y = trainingSet[,'AD'], ntree = ntree, mtry = mtry , importance = TRUE)
  saveRDS(output.forest, file = name3)
  
  trainADs <- makeListOfADsPerPatientFunction(unique(testset$user_id))
  predicted <- as.data.frame(predict(output.forest, testset, type = "class"))
  predicted <- cbind(predicted, trainADs)
  colnames(predicted) <- c("tempPred", "user_id", "real")
  predicted <- predicted[,c(2,1,3)]
  
  replaceFrame <- as.data.frame(cbind(trainADs$AD, trainADs$AD))
  colnames(replaceFrame) <- c("prediction", "tempPred")
  replaceFrame$tempPred <- as.factor(as.character(trainADs$AD))
  prediction <- merge(predicted, replaceFrame, by = "tempPred")
  prediction <- prediction[!duplicated(prediction), ]
  prediction <- prediction[,c(2,4,3)]
  
  sum <- 0
  disease <- c()
  for(i in 1:nrow(prediction)){
    sum <- sum + length(prediction$prediction[i][[1]])
    disease <- as.vector(sort(unique(c(disease, unlist(prediction$prediction[i][[1]])))))
  }
  print(paste0(name2," predicts ", length(disease), "  unique diseases"))
  print(paste0(name2, " on average predicts " , sum / nrow(prediction), " diseases per patient"))
  
  prediction <- scoreCalculationFunction(prediction)
  saveRDS(prediction, file = name4)
  scores <- returnScoreFunction(prediction, name2)
  performances <- calculatePerformanceFunction(prediction, name2, "no")
  outputs <- list(performances, scores)
  return(outputs)
}

set.seed(5)
flds <- createFolds(unique(train$user_id), k = 5, list = TRUE, returnTrain = FALSE)
ntreePara <- c(250, 500, 1000)
mtryPara <- c(15, sqrt((length(train)-1)), 30, 60, 100)
scores <- performances <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(ntreePara)){
    for(k in 1:length(mtryPara)){
      sampleTestIDs <- train[flds[[j]],]$user_id
      OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
      tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
      tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
      outputs <- randomForestFunction(tempTrainset, tempTestset, "5-fold", ntreePara[i], mtryPara[k])
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

endScores <- endPerformances <- c()
for(j in 1:length(ntreePara)){
  temp1 <- dplyr::filter(scores, (grepl(paste0("ntree = ", ntreePara[j]), scores$Name)))
  temp2 <- dplyr::filter(performances, (grepl(paste0("ntree = ", ntreePara[j]), performances$Name)))
  for(i in 1:length(mtryPara)){
    temp3 <- dplyr::filter(temp1, (grepl(paste0("mtry = ", mtryPara[i]), temp1$Name)))
    temp4 <- dplyr::filter(temp2, (grepl(paste0("mtry = ", mtryPara[i]), temp2$Name)))
    endScores <- rbind(endScores, c(temp3$Name[1], round(colSums(temp3[,-1])/nrow(temp3),4)))
    endPerformances <- rbind(endPerformances, c(temp4$Name[1], round(colSums(temp4[,-1])/nrow(temp4),4)))
  }
}

write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_5Fold.txt"), sep=" & ")
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_5Fold.txt"), sep=" & ")

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/end_5Fold.txt"), sep=" & ")


# Validation
endScores <- endPerformances <- c()
ntreePara <- c(250, 500)
mtryPara <- c(60)

for(i in 1: length(ntreePara)){
  for(j in 1:length(mtryPara)){
    outputs <- randomForestFunction(train, validation, "Validation", ntreePara[i], mtryPara[j])
    endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
    endScores <- rbind(endScores, as.data.frame(outputs[2]))
  }
}

outputs <- randomForestFunction(train, validation, "Validation", 1000, 30)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))

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
ntreePara <- c(250, 500)
mtryPara <- c(60)

for(i in 1: length(ntreePara)){
  for(j in 1:length(mtryPara)){
    outputs <- randomForestFunction(trainValidation, test, "Test", ntreePara[i], mtryPara[j])
    endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
    endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
  }
}

outputs <- randomForestFunction(trainValidation, test, "Test", 1000, 30)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- randomForestFunction(trainValidation, test, "Test", 500, 60)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))

endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))


write.csv(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endScores_Test.txt"), sep=" & ")
write.csv(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_Test.txt"), sep=" & ")

endScores <- cbind(endScores, endScores2)
endPerformances <- cbind(endPerformances, endPerformances2)
endScores <- endScores[,-6]
endPerformances <- endPerformances[,-5]
write.table(endScores, paste0(dir, "/Data/", modelsave, "mForest/EndScores-EndPerformances/endScores_ValTest.txt"), sep=" & ")
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances/endPerformances_ValTest.txt"), sep=" & ")

rm(list = setdiff(ls(), lsf.str()))


######################################################
## Binary Relevance Method
######################################################
dir <- getwd()
source(paste0(dir,"/Code/Import.R"))
modelsave <- "RandomForestBinary"
mkdir(paste0(dir, "/Data/", modelsave,"/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

BinaryRFFunction <- function(trainingSet, testSet, testUsers, ntree, name, mtry, threshold){
  name2 <- paste0("ntree = ", ntree , ", mtry =  ", mtry, ", threshold = ", threshold, " ", name)
  name3 <- paste0(dir, "Data/",modelsave, "/Models/",name2, ".rda")
  
  testADs <- makeListOfADsPerPatientFunction(testUsers)
  
  output.model <- ebr(trainingSet, "RF", ntree = ntree, mtry  = mtry, importance = TRUE, scale=FALSE, cores = 4, seed = 2018, attr.space = 1, subsample = 1, m = 10)
  prediction  <- as.probability(predict(output.model, testSet, cores = 4))
  
  
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
ntreePara <- c(250)
mtryPara <- c(15, sqrt((length(train)-1)), 30, 60, 100)
thresholdPara <- c(0.3,0.4,0.5)

performances <- scores <- c()
for(j in 2:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(ntreePara)){
    for(k in 1:length(mtryPara)){
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
        
        outputs <- BinaryRFFunction(tempTrainset, tempTestset, sampleTestIDs, ntreePara[i], "5-fold", mtryPara[k], thresholdPara[l])
        performances <- rbind(performances, as.data.frame(outputs[1]))
        scores <- rbind(scores, as.data.frame(outputs[2]))
      }
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

scores <- readRDS(paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/scores.5fold.rda"))
performances <- readRDS(paste0(dir, "/Data/", modelsave, "/Endscores-EndPerformances/performances.5fold.rda"))

endScores <- endPerformances <- c()
for(i in 1:length(ntreePara)){
  temp1 <- dplyr::filter(scores, (grepl(paste0("ntree = ", ntreePara[i]), scores$Name)))
  temp2 <- dplyr::filter(performances, (grepl(paste0("ntree = ", ntreePara[i]), performances$Name)))
  for(k in 1:length(mtryPara)){
    temp3 <- dplyr::filter(temp1, (grepl(paste0("mtry =  ", mtryPara[k]), temp1$Name)))
    temp4 <- dplyr::filter(temp2, (grepl(paste0("mtry =  ", mtryPara[k]), temp2$Name)))
    for(j in 1:length(thresholdPara)){
      temp5 <- dplyr::filter(temp3, (grepl(paste0("threshold = ", thresholdPara[j]), temp3$Name)))
      temp6 <- dplyr::filter(temp4, (grepl(paste0("threshold = ", thresholdPara[j]), temp4$Name)))
      endScores <- rbind(endScores, c(temp5$Name[1], round(colSums(temp5[,-1])/nrow(temp5),4)))
      endPerformances <- rbind(endPerformances, c(temp6$Name[1], round(colSums(temp6[,-1])/nrow(temp6),4)))
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

outputs <- BinaryRFFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), 250, "Validation", sqrt((length(train)-1)), 0.3)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- BinaryRFFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), 250, "Validation", 60, 0.4)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- BinaryRFFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), 250, "Validation", 100, 0.3)
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

outputs <- BinaryRFFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id), 250, "Test", sqrt((length(train)-1)), 0.3)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- BinaryRFFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id), 250, "Test", 60, 0.4)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- BinaryRFFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id), 250, "Test", 100, 0.3)
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
modelsave <- "RandomForestMultiLabel"
mkdir(paste0(dir, "/Data/", modelsave,"/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

multiLabelRandomForestFunction <- function(trainingSet, testset, name, ntree, mtry, labels, j){
  name2 <- paste0(" ntree = ", ntree, ", mtry = ", mtry, " ", name)
  if(j == 0 ){
    name3 <- paste0(dir,"/Data/", modelsave, "/Models/",name2, ".rda")
  }else{
    map <- paste0("5fold",j)
    name3 <- paste0(dir,"/Data/", modelsave, "/Models/",map,"/",name2, ".rda")
  }
  
  
  data.task <- makeMultilabelTask(id='multi', data = trainingSet[,-1], target = labels)
  learn.rfsrc <- makeLearner("multilabel.randomForestSRC", ntree = ntree, mtry= mtry, importance = TRUE)
  set.seed(2018)
  output.forest <- train(learn.rfsrc, data.task)
  saveRDS(output.forest, file = name3)
  
  testADs <- makeListOfADsPerPatientFunction(unique(testset$user_id))
  
  predicted <- predict(output.forest, newdata = testset)
  
  prediction <- as.data.frame(predicted)[,18:34]
  
  prediction$prediction <- NA
  for(j in 1:nrow(prediction)){
    temp <- c()
    if(prediction[j,17] == TRUE){
      temp <- c(temp, 0)
    }
    for(i in 1:16){
      if(prediction[j,i] == TRUE){
        temp <- c(temp, i)
      }
    }
    if(length(temp) == 0){
      temp <- 0
    }
    prediction$prediction[j] <- list(temp)
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

labels <- c(paste0("Disease.", seq(1:16)), paste0("Disease.", 0))

set.seed(5)
flds <- createFolds(unique(trainMultiLabel$user_id), k = 5, list = TRUE, returnTrain = FALSE)
ntreePara <- c(100, 250, 500, 1000)
mtryPara <- c(15, round(sqrt((length(train)-1))), 30, 60, 100, 250, 500)


scores <- performances <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(ntreePara)){
    for(k in 1:length(mtryPara)){
      sampleTestIDs <- train[flds[[j]],]$user_id
      OtherTrainIDs <- setdiff(unique(trainMultiLabel$user_id), sampleTestIDs)
      tempTestset <- dplyr::filter(trainMultiLabel, user_id %in% sampleTestIDs)
      tempTrainset <- dplyr::filter(trainMultiLabel, user_id %in% OtherTrainIDs)
      outputs <- multiLabelRandomForestFunction(tempTrainset, tempTestset, "5-fold", ntreePara[i], mtryPara[k], labels, j)
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}

scores <- readRDS(paste0(dir, "/Data/", modelsave, "/Endscores-Performances/scores.5fold.rda"))
performances <- readRDS(paste0(dir, "/Data/", modelsave, "/Endscores-Performances/performances.5fold.rda"))

saveRDS(scores, paste0(dir, "/Data/", modelsave, "/Endscores-Performances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/", modelsave, "/Endscores-Performances/performances.5fold.rda"))

endScores <- endPerformances <- c()
for(j in 1:length(ntreePara)){
  temp1 <- dplyr::filter(scores, (grepl(paste0("ntree = ", ntreePara[j], ','), scores$Name)))
  temp2 <- dplyr::filter(performances, (grepl(paste0("ntree = ", ntreePara[j], ","), performances$Name)))
  for(i in 1:length(mtryPara)){
    temp3 <- dplyr::filter(temp1, (grepl(paste0("mtry = ", mtryPara[i], " 5"), temp1$Name)))
    temp4 <- dplyr::filter(temp2, (grepl(paste0("mtry = ", mtryPara[i], " 5"), temp2$Name)))
    endScores <- rbind(endScores, c(temp3$Name[1], round(colSums(temp3[,-1])/nrow(temp3),4)))
    endPerformances <- rbind(endPerformances, c(temp4$Name[1], round(colSums(temp4[,-1])/nrow(temp4),4)))
  }
}

write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_5Fold.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_5Fold.txt"), sep=" & ", row.names=FALSE)

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/end_5Fold.txt"), sep=" & ", row.names=FALSE)


# Validation
labels <- c(paste0("Disease.", seq(1:16)), paste0("Disease.", 0))
endScores <- endPerformances <- c()

outputs <- multiLabelRandomForestFunction(trainMultiLabel, validationMultiLabel, "Validation", 100, 500, labels, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- multiLabelRandomForestFunction(trainMultiLabel, validationMultiLabel, "Validation", 250, 500, labels, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- multiLabelRandomForestFunction(trainMultiLabel, validationMultiLabel, "Validation", 500, 500, labels, 0)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))

endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))


write.csv(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_Validation.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_Validation.txt"), sep=" & ", row.names=FALSE)


# Test 
endScores2 <- endPerformances2 <- c()
trainValidation <- rbind(trainMultiLabel, validationMultiLabel)
trainValidation <- trainValidation[order(trainValidation$user_id),]

outputs <- multiLabelRandomForestFunction(trainValidation, testMultiLabel, "Test", 100, 500, labels, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- multiLabelRandomForestFunction(trainValidation, testMultiLabel, "Test", 250, 500, labels, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- multiLabelRandomForestFunction(trainValidation, testMultiLabel, "Test", 500, 500, labels, 0)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))

endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))


write.csv(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_Test.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_Test.txt"), sep=" & ", row.names=FALSE)

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endScores_ValTest.txt"), sep=" & ", row.names=FALSE)
write.table(endPerformances3, paste0(dir, "/Data/", modelsave, "/EndScores-Performances/endPerformances_ValTest.txt"), sep=" & ", row.names=FALSE)



rm(list = setdiff(ls(), lsf.str()))