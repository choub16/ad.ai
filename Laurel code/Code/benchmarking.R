dir <- getwd()
source(paste0(dir,"/Code/functions.R"))

mkdir(paste0(dir, "/Data/Benchmark/Endscores-Performances/"))

makeProbADs <- function(trainSet){
  tempUserADS <- dplyr::filter(usersWithAD, user_id %in% unique(trainSet$user_id))
  distributionADs <- as.data.frame(na.omit(c(as.vector(as.numeric(unlist(tempUserADS[,2:9]))))))
  distributionADs <- as.data.frame(table(distributionADs))
  return(distributionADs)
}

calculateMeanNoADs <- function(trainSet){
  meanNoADs <- sum(dplyr::filter(usersWithAD, user_id %in% unique(trainSet$user_id))[,10])/nrow(trainSet)
  return(meanNoADs)
}

makeBenchmark1 <- function(trainSet, testSet, name){
  real <- makeListOfADsPerPatientFunction(unique(testSet$user_id))
  prediction <- as.data.frame(real)
  prediction$prediction <- NA
  for(i in 1:nrow(prediction)){
    prediction$prediction[i] <- list(as.vector(unique(sample(0:16,sample(1:8, 1, replace=TRUE), replace = TRUE))))
  }
  prediction <- prediction[,c(1,3,2)]
  colnames(prediction) <- c("user_id", "prediction", "real")
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name)
  performances <- calculatePerformanceFunction(prediction, name, "no")
  
  outputs <- list(performances, scores)
  return(outputs)
}

makeBenchmark2 <- function(trainSet, testSet, name){
  real <- makeListOfADsPerPatientFunction(unique(testSet$user_id))
  prediction <- as.data.frame(real)
  prediction$prediction <- NA
  for(i in 1:nrow(prediction)){
    prediction$prediction[i] <- 1
  }
  prediction <- prediction[,c(1,3,2)]
  colnames(prediction) <- c("user_id", "prediction", "real")
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name)
  performances <- calculatePerformanceFunction(prediction, name, "no")
  
  outputs <- list(performances, scores)
  return(outputs)
}

makeBenchmark3 <- function(trainSet, testSet, distribution, name){
  real <- makeListOfADsPerPatientFunction(unique(testSet$user_id))
  prediction <- as.data.frame(real)
  prediction$prediction <- NA
  for(i in 1:nrow(prediction)){
    for(i in 1:nrow(prediction)){
      prediction$prediction[i] <- list(as.vector(unique(sample(0:16, size = sample(1:8,1,replace = TRUE), replace = TRUE, prob = distribution$prob))))
    }
  }
  prediction <- prediction[,c(1,3,2)]
  colnames(prediction) <- c("user_id", "prediction", "real")
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name)
  performances <- calculatePerformanceFunction(prediction, name, "no")
  
  outputs <- list(performances, scores)
  return(outputs)
}

makeBenchmark4 <- function(trainSet, testSet, distribution, meanNoADs, name){
  real <- makeListOfADsPerPatientFunction(unique(testSet$user_id))
  prediction <- as.data.frame(real)
  prediction$prediction <- NA
  for(i in 1:nrow(prediction)){
    for(i in 1:nrow(prediction)){
      (choose <- sample(1:2, size = 1))
      if(choose == 1){
        prediction$prediction[i] <- list(as.vector(unique(sample(0:16, size = floor(meanNoADs), replace = TRUE, prob = distribution$prob))))
      }else 
        if(choose == 2){
          prediction$prediction[i] <- list(as.vector(unique(sample(0:16, size = ceiling(meanNoADs), replace = TRUE, prob = distribution$prob))))
        }
      }
  }
  prediction <- prediction[,c(1,3,2)]
  colnames(prediction) <- c("user_id", "prediction", "real")
  
  prediction <- scoreCalculationFunction(prediction)
  scores <- returnScoreFunction(prediction, name)
  performances <- calculatePerformanceFunction(prediction, name, "no")
  
  outputs <- list(performances, scores)
  return(outputs)
}

usersWithAD$noADs <- (apply(usersWithAD, 1, function(x) sum(!is.na(x))))
usersWithAD$noADs <- usersWithAD$noADs-1

# 5-fold CV
scores <- performances <- c()
set.seed(5)
flds <- createFolds(unique(trainMultiLabel$user_id), k = 5, list = TRUE, returnTrain = FALSE)
for(j in 1:5){
  print(j)
  sampleTestIDs <- train[flds[[j]],]$user_id
  OtherTrainIDs <- setdiff(unique(train$user_id), sampleTestIDs)
  tempTestset <- dplyr::filter(train, user_id %in% sampleTestIDs)
  tempTrainset <- dplyr::filter(train, user_id %in% OtherTrainIDs)
  outputs <- makeBenchmark2(tempTrainset, tempTestset, "Benchmark2")
  performances <- rbind(performances, as.data.frame(outputs[1]))
  scores <- rbind(scores, as.data.frame(outputs[2]))
  tempTrainDist <- makeProbADs(tempTrainset)
  tempMeanNoADs <- calculateMeanNoADs(tempTrainset)
  for(i in 1:100){
    outputs <- makeBenchmark1(tempTrainset, tempTestset, "Benchmark1")
    performances <- rbind(performances, as.data.frame(outputs[1]))
    scores <- rbind(scores, as.data.frame(outputs[2]))
    outputs<- makeBenchmark3(tempTrainset, tempTestset, tempTrainDist, "Benchmark3")
    performances <- rbind(performances, as.data.frame(outputs[1]))
    scores <- rbind(scores, as.data.frame(outputs[2]))
    outputs<- makeBenchmark4(tempTrainset, tempTestset, tempTrainDist, tempMeanNoADs, "Benchmark4")
    performances <- rbind(performances, as.data.frame(outputs[1]))
    scores <- rbind(scores, as.data.frame(outputs[2]))
  }
}


saveRDS(scores, paste0(dir, "/Data/Benchmark/Endscores-Performances/scores.5fold.rda"))
saveRDS(performances, paste0(dir, "/Data/Benchmark/Endscores-Performances/performances.5fold.rda"))


endScores <- endPerformances <- c()
endScores <- rbind(endScores, c("Benchmark1", colMeans(dplyr::filter(scores, Name == "Benchmark1")[,-1])))
endScores <- rbind(endScores, c("Benchmark2", colMeans(dplyr::filter(scores, Name == "Benchmark2")[,-1])))
endScores <- rbind(endScores, c("Benchmark3", colMeans(dplyr::filter(scores, Name == "Benchmark3")[,-1])))
endScores <- rbind(endScores, c("Benchmark4", colMeans(dplyr::filter(scores, Name == "Benchmark4")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark1", colMeans(dplyr::filter(performances, Name == "Benchmark1")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark2", colMeans(dplyr::filter(performances, Name == "Benchmark2")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark3", colMeans(dplyr::filter(performances, Name == "Benchmark3")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark4", colMeans(dplyr::filter(performances, Name == "Benchmark4")[,-1])))

saveRDS(endScores, paste0(dir, "/Data/Benchmark/Endscores-Performances/endScores.5fold.rda"))
saveRDS(endPerformances, paste0(dir, "/Data/Benchmark/Endscores-Performances/endPerformances.5fold.rda"))


write.csv(endScores, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_5Fold.txt"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_5Fold.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_5Fold.txt"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_5Fold.txt"), sep=" & ", row.names=FALSE)

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/Benchmark/EndScores-Performances/end_5Fold.txt"), sep=" & ", row.names=FALSE)

#### VALIDATION

tempTrainDist <- makeProbADs(train)
tempMeanNoADs <- calculateMeanNoADs(train)

scores <- performances <- c()

outputs <- makeBenchmark2(train, validation, "Benchmark2")
performances <- rbind(performances, as.data.frame(outputs[1]))
scores <- rbind(scores, as.data.frame(outputs[2]))
for(i in 1:100){
  outputs <- makeBenchmark1(train, validation, "Benchmark1")
  performances <- rbind(performances, as.data.frame(outputs[1]))
  scores <- rbind(scores, as.data.frame(outputs[2]))
  outputs<- makeBenchmark3(train, validation, tempTrainDist, "Benchmark3")
  performances <- rbind(performances, as.data.frame(outputs[1]))
  scores <- rbind(scores, as.data.frame(outputs[2]))
  outputs<- makeBenchmark4(train, validation, tempTrainDist,tempMeanNoADs, "Benchmark4")
  performances <- rbind(performances, as.data.frame(outputs[1]))
  scores <- rbind(scores, as.data.frame(outputs[2]))
}

saveRDS(scores, paste0(dir, "/Data/Benchmark/Endscores-Performances/scores.Validation.rda"))
saveRDS(performances, paste0(dir, "/Data/Benchmark/Endscores-Performances/performances.Validation.rda"))


endScores <- endPerformances <- c()
endScores <- rbind(endScores, c("Benchmark1", colMeans(dplyr::filter(scores, Name == "Benchmark1")[,-1])))
endScores <- rbind(endScores, c("Benchmark2", colMeans(dplyr::filter(scores, Name == "Benchmark2")[,-1])))
endScores <- rbind(endScores, c("Benchmark3", colMeans(dplyr::filter(scores, Name == "Benchmark3")[,-1])))
endScores <- rbind(endScores, c("Benchmark4", colMeans(dplyr::filter(scores, Name == "Benchmark4")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark1", colMeans(dplyr::filter(performances, Name == "Benchmark1")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark2", colMeans(dplyr::filter(performances, Name == "Benchmark2")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark3", colMeans(dplyr::filter(performances, Name == "Benchmark3")[,-1])))
endPerformances <-rbind(endPerformances, c("Benchmark4", colMeans(dplyr::filter(performances, Name == "Benchmark4")[,-1])))

saveRDS(endScores, paste0(dir, "/Data/Benchmark/Endscores-Performances/endScores.Validation.rda"))
saveRDS(endPerformances, paste0(dir, "/Data/Benchmark/Endscores-Performances/endPerformances.Validation.rda"))


write.csv(endScores, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_Validation.txt"), row.names=FALSE)
write.table(endScores, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_Validation.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_Validation.txt"), row.names=FALSE)
write.table(endPerformances, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_Validation.txt"), sep=" & ", row.names=FALSE)

end <- cbind(endPerformances, endScores)
end <- end[,-5]
write.table(end, paste0(dir, "/Data/Benchmark/EndScores-Performances/end_Validation.txt"), sep=" & ", row.names=FALSE)


### Test
trainValidation <- rbind(trainMultiLabel, validationMultiLabel)
trainValidation <- trainValidation[order(trainValidation$user_id),]

tempTrainDist <- makeProbADs(trainValidation)
tempMeanNoADs <- calculateMeanNoADs(trainValidation)

scores2 <- performances2 <- c()

outputs <- makeBenchmark2(trainValidation, test, "Benchmark2")
performances2 <- rbind(performances2, as.data.frame(outputs[1]))
scores2 <- rbind(scores2, as.data.frame(outputs[2]))
for(i in 1:100){
  outputs <- makeBenchmark1(trainValidation, test, "Benchmark1")
  performances2 <- rbind(performances2, as.data.frame(outputs[1]))
  scores2 <- rbind(scores2, as.data.frame(outputs[2]))
  outputs<- makeBenchmark3(trainValidation, test, tempTrainDist, "Benchmark3")
  performances2 <- rbind(performances2, as.data.frame(outputs[1]))
  scores2 <- rbind(scores2, as.data.frame(outputs[2]))
  outputs<- makeBenchmark4(trainValidation, test, tempTrainDist,tempMeanNoADs, "Benchmark4")
  performances2 <- rbind(performances2, as.data.frame(outputs[1]))
  scores2 <- rbind(scores2, as.data.frame(outputs[2]))
}

saveRDS(scores2, paste0(dir, "/Data/Benchmark/Endscores-Performances/scores.Test.rda"))
saveRDS(performances2, paste0(dir, "/Data/Benchmark/Endscores-Performances/performances.Test.rda"))

endScores2 <- endPerformances2 <- c()
endScores2 <- rbind(endScores2, c("Benchmark1", colMeans(dplyr::filter(scores2, Name == "Benchmark1")[,-1])))
endScores2 <- rbind(endScores2, c("Benchmark2", colMeans(dplyr::filter(scores2, Name == "Benchmark2")[,-1])))
endScores2 <- rbind(endScores2, c("Benchmark3", colMeans(dplyr::filter(scores2, Name == "Benchmark3")[,-1])))
endScores2 <- rbind(endScores2, c("Benchmark4", colMeans(dplyr::filter(scores2, Name == "Benchmark4")[,-1])))
endPerformances2 <-rbind(endPerformances2, c("Benchmark1", colMeans(dplyr::filter(performances2, Name == "Benchmark1")[,-1])))
endPerformances2 <-rbind(endPerformances2, c("Benchmark2", colMeans(dplyr::filter(performances2, Name == "Benchmark2")[,-1])))
endPerformances2 <-rbind(endPerformances2, c("Benchmark3", colMeans(dplyr::filter(performances2, Name == "Benchmark3")[,-1])))
endPerformances2 <-rbind(endPerformances2, c("Benchmark4", colMeans(dplyr::filter(performances2, Name == "Benchmark4")[,-1])))

saveRDS(endScores2, paste0(dir, "/Data/Benchmark/Endscores-Performances/endScores_Test.rda"))
saveRDS(endPerformances2, paste0(dir, "/Data/Benchmark/Endscores-Performances/endPerformances_Test.rda"))


write.csv(endScores2, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_Test.txt"), row.names=FALSE)
write.table(endScores2, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_Test.txt"), sep=" & ", row.names=FALSE)
write.csv(endPerformances2, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_Test.txt"), row.names=FALSE)
write.table(endPerformances2, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_Test.txt"), sep=" & ", row.names=FALSE)

endScores3 <- cbind(endScores, endScores2)
endPerformances3 <- cbind(endPerformances, endPerformances2)
endScores3 <- endScores3[,-6]
endPerformances3 <- endPerformances3[,-5]
write.table(endScores3, paste0(dir, "/Data/Benchmark/EndScores-Performances/endScores_ValTest.txt"), sep=" & ", row.names=FALSE)
write.table(endPerformances3, paste0(dir, "/Data/Benchmark/EndScores-Performances/endPerformances_ValTest.txt"), sep=" & ", row.names=FALSE)


