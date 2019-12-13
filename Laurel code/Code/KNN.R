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
modelsave <- "KNNLabelPowerset"
mkdir(paste0(dir, "/Data/", modelsave,"/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

KNNFunction <- function(traindata, testdata, method, name, K, frequency){
  userList <- makeListOfADsPerPatientFunction(unique(traindata$user_id))
  dataset <- rbind(traindata, testdata)
  row.names(dataset) <- c(traindata$user_id, testdata$user_id)
  name <- paste0(method, " , ", K, "-nearest neighbor, frequency = ", frequency, " ,", name)
  name2 <- paste0(dir,"/Models/", name, ".rda")
  
  
  if(method %in% c("euclidean")){    # For computational time reasons!
    euc_dist <- as.matrix(dist(dataset[1:nrow(dataset),-1], method = method), labels = TRUE)
    colnames(euc_dist) <- rownames(euc_dist) <- dataset[['user_id']]
  }else{
    euc_dist <- as.matrix(distance(as.data.frame(dataset[1:nrow(dataset),-1]), method = method), labels = TRUE)
    colnames(euc_dist) <- rownames(euc_dist) <- dataset[['user_id']]
  }
  
  euc_dist <- subset(euc_dist, row.names(euc_dist) %in% unique(testdata$user_id))
  euc_dist <- euc_dist[, !colnames(euc_dist) %in% unique(testdata$user_id)]
  
  colNames <- c("user_id", "similar", "prediction", "real")
  prediction <- as.data.frame(setNames(replicate(4,numeric(0), simplify = T), colNames))
  prediction[1:nrow(euc_dist),] <- NA
  
  for(i in 1:nrow(euc_dist)){
    patient <- as.numeric(rownames(euc_dist)[i])
    tempMin <- as.vector(as.numeric(names(apply(euc_dist,1,function(x) which(x==min(x)))[i][[1]])))
    if(K <= length(tempMin)){
      if(length(tempMin)>K){
        set.seed(1)
        similar_patients <- as.vector(sample(tempMin,K))
      }else{
        similar_patients <- as.vector(tempMin)
      }
    }else{
      similar_patients <- as.vector(as.numeric(names(sort(euc_dist[i,], decreasing = FALSE)[1:K])))
    }
    
    prediction$user_id[i] <- patient
    prediction$similar[i] <- list(similar_patients)
    
    if(K != 1){
      prediction$prediction[i] <- list(unique(Filter(function (elem) length(which(
        as.vector(as.numeric(unlist(dplyr::filter(userList, 
                                                  user_id %in% similar_patients)[,-1]))) == elem)) > frequency, 
        as.vector(as.numeric(unlist(dplyr::filter(userList, 
                                                  user_id %in% similar_patients)[,-1]))) 
      )))
    }else{
      prediction$prediction[i] <- list(unique(as.numeric(unlist(dplyr::filter(userList, 
                                                                              user_id %in% similar_patients)[,-1]))))
    }
    prediction$real[i] <- list(na.omit(as.numeric(dplyr::filter(usersWithAD, user_id == prediction$user_id[i])[,-1])))
  }
  
  saveRDS(prediction, file = name2)
  prediction <- scoreCalculationFunction(prediction[,c('user_id', 'prediction', 'real')])
  scores2 <- returnScoreFunction(prediction, name)
  performances <- calculatePerformanceFunction(prediction, name, "no")
  
  sum <- 0
  disease <- c()
  for(i in 1:nrow(prediction)){
    sum <- sum + length(prediction$prediction[i][[1]])
    disease <- as.vector(sort(unique(c(disease, unlist(prediction$prediction[i][[1]])))))
  }
  print(paste0(name," predicts ", length(disease), "  unique diseases"))
  print(paste0(name, " on average predicts " , sum / nrow(prediction), " diseases per patient"))
  
  output <- list(performances, scores2, prediction)
  return(output)
}

distanceMethods = c("euclidean", "jaccard", "pearson")
KPara <- c(1,3,5,7,10)
freq <- c(1,2)

# 5-Fold Cross Validation
set.seed(5)
flds <-createFolds(unique(train$user_id), k = 5, list = TRUE, returnTrain = FALSE)
performances <- scores <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(k in 1:length(distanceMethods)){
    print(paste0(k, " of 3"))
    for(i in 1:length(KPara)){
      for(l in 1:length(freq)){
        sampleTestIDs <- train[flds[[j]],]$user_id
        OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
        tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
        tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
        outputs <- KNNFunction(tempTrainset, tempTestset, distanceMethods[k], "5-fold", KPara[i], freq[l])
        performances <- rbind(performances, as.data.frame(outputs[1]))
        scores <- rbind(scores, as.data.frame(outputs[2]))
      }
    }
  }
}

saveRDS(scores, paste0(dir,"/Data/",modelsave,"/EndScores-Performances/scores.5Fold.rda"))
saveRDS(performances, paste0(dir,"/Data/",modelsave,"/EndScores-Performances/performances.5Fold.rda"))

scores <- readRDS(paste0(dir,"/Data/",modelsave,"/EndScores-Performances/scores.5Fold.rda"))
performances <- readRDS(paste0(dir,"/Data/",modelsave,"/EndScores-Performances/performances.5Fold.rda"))

endPerformances <- endScores <- c()
for(k in 1:length(distanceMethods)){
  temp1 <- dplyr::filter(scores, (grepl(distanceMethods[k], scores$Name)))
  temp3 <- dplyr::filter(performances, (grepl(distanceMethods[k], performances$Name)))
  for(j in 1:length(KPara)){
    temp2 <- dplyr::filter(temp1, (grepl(paste0(KPara[j], "-nearest"), temp1$Name)))
    temp4 <- dplyr::filter(temp3, (grepl(paste0(KPara[j], "-nearest"), temp3$Name)))
    for(i in 1:length(freq)){
      temp5 <- dplyr::filter(temp2, (grepl(paste0("frequency = ", freq[i]), temp2$Name)))
      temp6 <- dplyr::filter(temp4, (grepl(paste0("frequency = ", freq[i]), temp4$Name)))
      endScores <- rbind(endScores, c(temp5$Name, round(colSums(temp5[,-1])/nrow(temp5),4)))
      endPerformances <- rbind(endPerformances, c(temp6$Name, round(colSums(temp6[,-1])/nrow(temp6),4)))
    }
  }
}

endScores <- endScores[,5:9]
endPerformances <- endPerformances[,5:8]
write.csv(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_5Fold.txt"), sep=" & ")
write.csv(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_5Fold.txt"), sep=" & ")
                            

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir,"/Data/",modelsave,"/EndScores-Performances/end_5Fold.txt"), sep=" & ")

# Validation
distanceMethods = c("euclidean", "jaccard")
KPara <- c(1,10)
freq <- 1

endScores <- endPerformances <- c()
for(k in 1:length(distanceMethods)){
  print(paste0(k, " of 4"))
  for(i in 1:length(KPara)){
    outputs <- KNNFunction(train, validation, distanceMethods[k], "Validation", KPara[i], freq)
    endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
    endScores <- rbind(endScores, as.data.frame(outputs[2]))
  }
}

endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))

write.csv(endScores, paste0("/Data/",modelsave,"/EndScores-Performances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0("/Data/",modelsave,"/EndScores-Performances/endScores_Validation.txt"), sep=" & ")
write.csv(endPerformances, paste0("/Data/",modelsave,"/EndScores-Performances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0("/Data/",modelsave,"/EndScores-Performances/endPerformances_Validation.txt"), sep=" & ")


# Test
distanceMethods = c("euclidean", "jaccard")
KPara <- c(1,10)
freq <- 1

endScores2 <- endPerformances2 <- c()
trainValidation <- rbind(train, validation)
trainValidation <- trainValidation[order(trainValidation$user_id),]
for(k in 1:length(distanceMethods)){
  print(paste0(k, " of 4"))
  for(i in 1:length(KPara)){
    outputs <- KNNFunction(trainValidation, test, distanceMethods[k], "Test", KPara[i], freq)
    endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
    endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
  }
}



endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))


write.csv(endScores2, paste0(dir,"/Data/",modelsave,"/EndScores-Performances/endScores_Test.csv", row.names=FALSE))
write.table(endScores2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_Test.txt"), sep=" & ")
write.csv(endPerformances2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_Test.txt"), sep=" & ")


rm(list = setdiff(ls(), lsf.str()))

######################################################
## Binary Relevance Method
######################################################
dir <- getwd()
source(paste0(dir,"/Code/Import.R"))
modelsave <- "KNNBinary"
mkdir(paste0(dir, "/Data/", modelsave,"/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

BinaryKNNFunction <- function(trainingSet, testSet, testUsers, name, K, threshold){
  name2 <- paste0(" Method = Euclidean", ", K = ", K , ", threshold = ", threshold, " ", name)
  name3 <- paste0(dir,"/Data/",modelsave,"/Models/",name2, ".rda")
  
  testADs <- makeListOfADsPerPatientFunction(testUsers)
  
  output.model <- ebr(trainingSet, "KNN", k = K, scale=FALSE, cores = 4, seed = 2018, attr.space = 1, subsample = 1, m = 10)
  prediction  <- predict(output.model, testSet, cores = 4)
  
  saveRDS(output.model, file = name3)
  
  prediction <- as.data.frame(as.matrix(as.probability(prediction)))
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

tempDisease <- dplyr::filter(usersWithAD, user_id %in% validation$user_id)
# 5-Fold Cross Validation
set.seed(5)
flds <- createFolds(unique(trainMultiLabel$user_id), k = 5, list = TRUE, returnTrain = FALSE)
KPara <- c(1,3,5,7,10, 12, 15, 20, 25, 30, 35, 40)
thresholdPara <- c(0.3,0.4,0.5)

performances <- scores <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(KPara)){
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
      
      outputs <- BinaryKNNFunction(tempTrainset, tempTestset, sampleTestIDs, "5-fold", KPara[i], thresholdPara[l])
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}

saveRDS(scores, paste0(dir, "/Data/",modelsave,"/Endscores-EndPerformances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/",modelsave,"/Endscores-EndPerformances/performances.5fold.rda"))

scores <- readRDS(paste0(dir, "/Data/",modelsave,"/Endscores-EndPerformances/scores.5fold.rda"))
performances <- readRDS(paste0(dir, "/Data/",modelsave,"/Endscores-EndPerformances/performances.5fold.rda"))

endScores <- endPerformances <- c()
for(i in 1:length(KPara)){
  temp3 <- dplyr::filter(scores, (grepl(paste0("K = ", KPara[i], ","), scores$Name)))
  temp4 <- dplyr::filter(performances, (grepl(paste0("K = ", KPara[i], ","), performances$Name)))
  for(j in 1:length(thresholdPara)){
    temp5 <- dplyr::filter(temp3, (grepl(paste0("threshold = ", thresholdPara[j], " 5"), temp3$Name)))
    temp6 <- dplyr::filter(temp4, (grepl(paste0("threshold = ", thresholdPara[j], " 5"), temp4$Name)))
    endScores <- rbind(endScores, c(temp5$Name[1], round(colSums(temp5[,-1])/nrow(temp5),4)))
    endPerformances <- rbind(endPerformances, c(temp6$Name[1], round(colSums(temp6[,-1])/nrow(temp6),4)))
  }
}

write.csv(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_5Fold.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_5Fold.txt"), sep=" & ", row.names=FALSE)

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/end_5Fold.txt"), sep=" & ", row.names=FALSE)

#### Validation

cols <- sapply(trainMultiLabel, is.logical)
trainMultiLabel[,cols] <- lapply(trainMultiLabel[,cols], as.numeric)
cols <- sapply(validationMultiLabel, is.logical)
validationMultiLabel[,cols] <- lapply(validationMultiLabel[,cols], as.numeric)

trainBinary <- mldr_from_dataframe(trainMultiLabel[,-1], labelIndices =  c(477:493))
validationBinary <- mldr_from_dataframe(validationMultiLabel[,-1], labelIndices =  c(477:493))

endScores <- endPerformances <- c()

outputs <- BinaryKNNFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id), "Validation", 7, 0.4)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- BinaryKNNFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id),  "Validation", 12, 0.4)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- BinaryKNNFunction(trainBinary, validationBinary, unique(validationMultiLabel$user_id),  "Validation", 30, 0.3)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))

endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))


write.csv(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_Validation.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_Validation.txt"), sep=" & ", row.names=FALSE)

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

outputs <- BinaryKNNFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id), "Test", 7, 0.4)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- BinaryKNNFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id),  "Test", 12, 0.4)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- BinaryKNNFunction(trainValidationBinary, testBinary, unique(testMultiLabel$user_id),  "Test", 30, 0.3)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))

endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))

write.csv(endScores2, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_Test.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances2, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_Test.txt"), sep=" & ", row.names=FALSE)

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endScores_ValTest.txt"), sep=" & ", row.names=FALSE)
write.table(endPerformances3, paste0(dir, "/Data/",modelsave,"/EndScores-EndPerformances/endPerformances_ValTest.txt"), sep=" & ", row.names=FALSE)

rm(list = setdiff(ls(), lsf.str()))

######################################################
## Algorithm Adaptation Method
######################################################
dir <- getwd()
source(paste0(dir,"/Code/Import.R"))
modelsave <- "KNNMultiLabel"
mkdir(paste0(dir, "/Data/", modelsave, "/Models"))
mkdir(paste0(dir, "/Data/", modelsave, "/EndScores-EndPerformances"))

multiLabelKNNFunction <- function(trainingSet, testset, testUsers, name, K, threshold){
  
  name2 <- paste0(" K = ", K, ", threshold = ", threshold, " ", name)
  name3 <- paste0(dir,"/Data/",modelsave,"/Models/",name2, ".rda")
  
  testADs <- makeListOfADsPerPatientFunction(testUsers)
  
  output.model <- mlknn(trainingSet, k = K, seed = 2018)
  prediction <- predict(output.model, testset)
  
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

set.seed(5)
flds <- createFolds(unique(trainMultiLabel$user_id), k = 5, list = TRUE, returnTrain = FALSE)
KPara <- c(1,3,5,7,10, 12, 15, 20, 25, 30, 35, 40)
thresholdPara <- c(0.1,0.2,0.25,0.3, 0.35)
scores <- performances <- c()
for(j in 1:5){
  print(paste0(j, " of 5"))
  for(i in 1:length(KPara)){
    for(k in 1:length(thresholdPara)){
      sampleTestIDs <- train[flds[[j]],]$user_id
      OtherTrainIDs <- setdiff(unique(trainMultiLabel$user_id), sampleTestIDs)
      tempTempTestset <- as.data.frame(dplyr::filter(trainMultiLabel, user_id %in% sampleTestIDs))
      tempTempTrainset <- as.data.frame(dplyr::filter(trainMultiLabel, user_id %in% OtherTrainIDs))
      tempTrainset <- mldr_from_dataframe(tempTempTrainset[,-1], labelIndices = c(477:493))
      tempTestset <- mldr_from_dataframe(tempTempTestset[,-1], labelIndices = c(477:493))
      
      outputs <- multiLabelKNNFunction(tempTrainset, tempTestset, unique(tempTempTestset$user_id), "5-fold", KPara[i], thresholdPara[k])
      performances <- rbind(performances, as.data.frame(outputs[1]))
      scores <- rbind(scores, as.data.frame(outputs[2]))
    }
  }
}


saveRDS(scores, paste0(dir, "/Data/",modelsave,"/Endscores-Performances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/",modelsave,"/Endscores-Performances/performances.5fold.rda"))

scores <- readRDS(paste0(dir, "/Data/",modelsave,"/Endscores-Performances/scores.5fold.rda"))
performances <- readRDS(paste0(dir, "/Data/",modelsave,"/Endscores-Performances/performances.5fold.rda"))

endScores <- endPerformances <- c()
for(j in 1:length(KPara)){
  temp1 <- dplyr::filter(scores, (grepl(paste0("K = ", KPara[j], ","), scores$Name)))
  temp2 <- dplyr::filter(performances, (grepl(paste0("K = ", KPara[j], ","), performances$Name)))
  for(i in 1:length(thresholdPara)){
    temp3 <- dplyr::filter(temp1, (grepl(paste0("threshold = ", thresholdPara[i], " 5"), temp1$Name)))
    temp4 <- dplyr::filter(temp2, (grepl(paste0("threshold = ", thresholdPara[i], " 5"), temp2$Name)))
    endScores <- rbind(endScores, c(temp3$Name[1], round(colSums(temp3[,-1])/nrow(temp3),4)))
    endPerformances <- rbind(endPerformances, c(temp4$Name[1], round(colSums(temp4[,-1])/nrow(temp4),4)))
  }
}

write.csv(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_5Fold.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_5Fold.txt"), sep=" & ",  row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_5Fold.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_5Fold.txt"), sep=" & ",  row.names=FALSE)

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/end_5Fold.txt"), sep=" & ", row.names=FALSE)

# Validation
endScores <- endPerformances <- c()
trainML <- mldr_from_dataframe(trainMultiLabel[,-1], labelIndices = c(477:493))
validationML <- mldr_from_dataframe(validationMultiLabel[,-1], labelIndices = c(477:493))

outputs <- multiLabelKNNFunction(trainML, validationML, unique(validationMultiLabel$user_id), "Validation", 15, 0.25)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- multiLabelKNNFunction(trainML, validationML, unique(validationMultiLabel$user_id), "Validation", 25, 0.25)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))
outputs <- multiLabelKNNFunction(trainML, validationML, unique(validationMultiLabel$user_id), "Validation", 35, 0.1)
endPerformances <- rbind(endPerformances, as.data.frame(outputs[1]))
endScores <- rbind(endScores, as.data.frame(outputs[2]))


endPerformances <- cbind(endPerformances$Name, round(endPerformances[,-1],4))
colnames(endPerformances) <- c("Name", colnames(endPerformances[,-1]))
endScores <- cbind(endScores$Name, round(endScores[,-1],4))
colnames(endScores) <- c("Name", colnames(endScores[,-1]))


write.csv(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_Validation.csv"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_Validation.txt"), sep=" & ",  row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_Validation.csv"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_Validation.txt"), sep=" & ",  row.names=FALSE)


# Test 
tempTrainValidationML <- rbind(trainMultiLabel, validationMultiLabel)
tempTrainValidationML <- tempTrainValidationML[order(tempTrainValidationML$user_id),]
trainValidationML <- mldr_from_dataframe(tempTrainValidationML[,-1], labelIndices = c(477:493))
testML <- mldr_from_dataframe(testMultiLabel[,-1], labelIndices = c(477:493))

endScores2 <- endPerformances2 <- c()

outputs <- multiLabelKNNFunction(trainValidationML, testML, unique(testMultiLabel$user_id), "Test", 15, 0.25)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- multiLabelKNNFunction(trainValidationML, testML, unique(testMultiLabel$user_id), "Test", 25, 0.25)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))
outputs <- multiLabelKNNFunction(trainValidationML, testML, unique(testMultiLabel$user_id), "Test", 35, 0.1)
endPerformances2 <- rbind(endPerformances2, as.data.frame(outputs[1]))
endScores2 <- rbind(endScores2, as.data.frame(outputs[2]))

endPerformances2 <- cbind(endPerformances2$Name, round(endPerformances2[,-1],4))
colnames(endPerformances2) <- c("Name", colnames(endPerformances2[,-1]))
endScores2 <- cbind(endScores2$Name, round(endScores2[,-1],4))
colnames(endScores2) <- c("Name", colnames(endScores2[,-1]))

write.csv(endScores2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_Test.csv"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_Test.txt"), sep=" & ",  row.names=FALSE)
write.csv(endPerformances2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_Test.csv"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_Test.txt"), sep=" & ",  row.names=FALSE)

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endScores_ValTest.txt"), sep=" & ",  row.names=FALSE)
write.table(endPerformances3, paste0(dir, "/Data/",modelsave,"/EndScores-Performances/endPerformances_ValTest.txt"), sep=" & ",  row.names=FALSE)

rm(list = setdiff(ls(), lsf.str()))

