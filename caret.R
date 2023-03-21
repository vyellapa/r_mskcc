library("caret")
library("kernlab")
library("ISLR")
library("Hmisc")
data(Wage)
data(spam)
inTrain = createDataPartition(y=spam$type, p=0.6,list=FALSE)
train = spam[inTrain,]
test = spam[-inTrain,]

modelFit = train(type ~ ., data = training, method = "glm")

folds = createFolds(y=spam$type, k = 10,list=TRUE, returnTrain = TRUE)

library("ISLR")
data(Wage)
inTrain = createDataPartition(y=Wage$wage, p=0.6,list=FALSE)

train = Wage[inTrain,]
test = Wage[-inTrain,]

#Cut based on quantile groups
cutWage = cut2(train$wage, g = 5)
train$cutWage = cutWage
table(cutWage)