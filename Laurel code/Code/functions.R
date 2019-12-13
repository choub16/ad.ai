######################################################
######################################################
## Laurel Joannes - laureljoannes@hotmail.com       ##
## November 2018                                    ##
## This script contains all code for libraries      ##
##  functions that are used all over                ##
######################################################
######################################################

packageList <- c("SnowballC",
                 "readr",
                 "plyr",
                 "dplyr",
                 "wordcloud",
                 "tm",
                 "ggplot2",
                 "ggpubr",
                 "scales",
                 "readxl",
                 "cluster",
                 "philentropy",
                 "caret",
                 "randomForest",
                 "rlang",
                 "Rtsne",
                 "e1071",
                 "class",
                 "RColorBrewer",
                 "gridExtra",
                 "mlr",
                 "utiml",
                 "mldr",
                 "klaR", 
                 "icesTAF"
)

lapply(packageList, require, character.only = TRUE)




makeListOfADsPerPatientFunction <- function(dataIds) {
  dataset <- dplyr::filter(usersWithAD, user_id %in% dataIds)
  dataset$AD <- NA
  for(k in 1:nrow(dataset)){
    dataset$AD[k] <- list(na.omit(as.vector(as.numeric(dataset[k,-1]))))
  }
  users <- dataset[,c(1,10)]
  return(users)
}
scoreCalculationFunction <- function(prediction){ #gives for every row a score, every patient
  for(i in 1:nrow(prediction)){
    prediction$percentage[i] <- max(0,(ifelse(((sum(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]) * 100))>100), 100,(sum(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]) * 100)))))
    prediction$score_without[i] <- max(0,(sum(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]) * 100) /
                                            length(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]))))
    prediction$score_penalty[i] <- max(0,(sum(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]) * 100) / 
                                            length(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]))) -
                                         sum(!unlist(prediction$prediction[i]) %in% unlist(prediction$real[i]))*10)
    prediction$score_scaled[i] <- max(0,(sum(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]) * 100) / 
                                           length(unlist(prediction$real[i]) %in% unlist(prediction$prediction[i]))) -
                                        sum((!unlist(prediction$prediction[i]) %in% unlist(prediction$real[i]))*20) /
                                        length(unlist(prediction$real[i])))
  }
  return(prediction) 
}
returnScoreFunction <- function(prediction, Name){ # gives an overall score
  colNames <- c("Name", "Percentage","Score_without", "Score_penalty", "Score_scaled")
  scores <- as.data.frame(setNames(replicate(5,numeric(0), simplify = T), colNames))
  scores[1,] <- NA
  scores$Name[1] <- Name
  scores$Percentage[1] <- sum(prediction$percentage) / nrow(prediction)
  scores$Score_without[1] <- sum(prediction$score_without) /nrow(prediction)
  scores$Score_penalty[1] <- sum(prediction$score_penalty) /nrow(prediction)
  scores$Score_scaled[1] <- sum(prediction$score_scaled) /nrow(prediction)
  return(scores)
}
makeAdsPerPatientFunction <- function(dataIds){
  dataset <- dplyr::filter(usersWithAD, user_id %in% dataIds)
  users <- dataset[,1:2]
  names(users) <- c("user_id", "AD")
  for(i in 3:ncol(dataset)){
    userTemp <- dataset[,c(1,i)]
    userTemp <- userTemp[complete.cases(userTemp),]
    names(userTemp) <- c("user_id", "AD")
    users <- rbind(users, userTemp)
  }
  users <- users[order(users$user_id),]
  return(users)
}
calculatePerformanceFunction <- function(prediction, name, bench){
  total <- merge(prediction[,c('user_id', 'prediction')],makeAdsPerPatientFunction(unique(prediction$user_id)),by="user_id")
  if(bench != "Yes"){
    colNames <- c("user_id", "prediction","AD")
    total2 <- as.data.frame(setNames(replicate(3,numeric(0), simplify = T), colNames))
    total2[1,] <- NA
    for(k in 1:nrow(total)){
      for(j in 1:length(unlist(total$prediction[k][[1]]))){
        total2 <- rbind(total2, c(total$user_id[k], total$prediction[k][[1]][j], total$AD[k]))
      }
    }
    total2 <- total2[-1,]
    total <- total2
  }
  partly <- table(total$prediction, total$AD)
  if(!is_empty(setdiff(total$prediction, total$AD))){
    for(j in 1:length(setdiff(total$prediction, total$AD))){
      partly <- cbind(partly, 0)
      colnames(partly) <- c(colnames(partly)[-length(colnames(partly))], as.character(setdiff(total$prediction, total$AD)[j]))
    }
  }
  if(!is_empty(setdiff(total$AD, total$prediction))){
    for(j in 1:length(setdiff(total$AD, total$prediction))){
      partly <- rbind(partly, 0)
      rownames(partly) <- c(rownames(partly)[-length(rownames(partly))], as.character(setdiff(total$AD, total$prediction)[j]))
    }
  }
  partly <- partly[, as.character(sort(as.numeric(colnames(partly))))]
  partly <- partly[as.character(sort(as.numeric(rownames(partly)))),]
  partly_conf <- confusionMatrix(as.table(partly))
  partly_confClass <- data.frame(partly_conf$byClass)
  partly_confClass[is.na(partly_confClass)] <- 0
  
  colNames <- c("Name", "Macro_Average_Precision","Macro_Average_Recall", "Macro_Average_F-measure")
  performances <-  as.data.frame(setNames(replicate(4,numeric(0), simplify = T), colNames))
  performances[1,] <- NA
  performances$Name[1] <- name
  performances$Macro_Average_Precision[1] <- sum(partly_confClass['Pos.Pred.Value'])/nrow(partly_confClass['Pos.Pred.Value'])
  performances$Macro_Average_Recall[1] <- sum(partly_confClass['Sensitivity'])/nrow(partly_confClass['Sensitivity'])
  performances$Macro_Average_F.measure[1] <- 2*performances$Macro_Average_Precision[1]*performances$Macro_Average_Recall[1]/(performances$Macro_Average_Precision[1]+performances$Macro_Average_Recall[1])
  return(performances)
}



