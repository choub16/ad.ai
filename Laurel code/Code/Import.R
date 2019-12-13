######################################################
######################################################
## Laurel Joannes - laureljoannes@hotmail.com       ##
## November 2018                                    ##
## This script contains all code to read data       ##
######################################################
######################################################

train <- read_csv("/Data/train.csv")
test <- read_csv("/Data/test.csv")
validation <- read_csv("/Data/validation.csv")
usersWithAD <- read_csv("/Data/usersWithAD.csv")
ADList <- read_csv("/Data/ADList.csv")

trainMultiLabel <- read_csv("/Data/trainMultiLabel.csv")
testMultiLabel <- read_csv("/Data/testMultiLabel.csv")
validationMultiLabel <- read_csv("/Data/validationMultiLabel.csv")