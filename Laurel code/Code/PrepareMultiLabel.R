dir <- getwd()

train <- read_csv(paste0(dir, "/Data/train.csv"))
test <- read_csv(paste0(dir, "/Data/test.csv"))
validation <- read_csv(paste0(dir, "/Data/validation.csv"))
usersWithAD <- read_csv(paste0(dir, "/Data/usersWithAD.csv"))

allData <- rbind(train, validation, test)
allData <- allData[order(allData$user_id),]

allDataDiseases <- as.matrix(dplyr::filter(usersWithAD, user_id %in% allData$user_id))

userDiseases <- setNames(data.frame(matrix(ncol = 18, nrow= nrow(allData))), c("user_id",paste0("Disease.", seq(1:16)), paste0("Disease.", 0)))
userDiseases$user_id <- allData$user_id
for(i in 1:nrow(allDataDiseases)){
  temp <- as.vector(na.omit(as.numeric(unique(allDataDiseases[i,-1]))))
  if(0 %in% temp){
    userDiseases[i,18] <- TRUE
    temp <- setdiff(temp, 0)
  }
  if(length(temp) > 0){
    for(j in 1:length(temp)){
      userDiseases[i,temp[j]+1] <- TRUE
    }
  }
}

userDiseases[is.na(userDiseases)] <- FALSE
allData <- merge(allData, userDiseases, on = 'user_id')

train2 <- dplyr::filter(allData, user_id %in% train$user_id)
validation2 <- dplyr::filter(allData, user_id %in% validation$user_id)
test2 <- dplyr::filter(allData, user_id %in% test$user_id)

write.csv(train2, paste0(dir, "/trainMultiLabel.csv"), row.names = FALSE)
write.csv(validation2, paste0(dir, "/validationMultiLabel.csv"), row.names = FALSE)
write.csv(test2, paste0(dir, "/testMultiLabel.csv"), row.names = FALSE)
