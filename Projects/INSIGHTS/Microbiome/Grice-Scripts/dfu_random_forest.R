# Random Forest on DFU
library(loescher)
library(dplyr)
library(caret)
library(randomForest)
library(ROCR)

# Read in file
setwd('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/')
load('./blast.Rdata')

# Normalize taxa table
tax <- apply(taxa$Species, 2, function(x) x/sum(x)) %>% t

# Remove rare taxa to make computationally feasible on laptop
tax <- tax[,order(colSums(tax), decreasing = T)[1:50]]

# Make taxa table in same order as dfu table, remove other samples
tax <- tax[dfu$SampleID[dfu$has_sample],]

# OTU: Baseline - Complications ---------------------------------------------------------------

# Subset to samples at baseline
sub <- factor(dfu$comp[dfu$has_sample], labels = c('NoComp','Comp'))
sub <- cbind(outcome = sub, as.data.frame(tax))
sub <- sub[dfu$Visit[dfu$has_sample] == 0,]

# Make training and test sets
set.seed(1234)
inTrain <-  createDataPartition(y = sub$outcome, p = .75, list = F)
training <- sub[inTrain,]
testing <- sub[-inTrain,]

# Define training parameters
ctrl <- trainControl(method = 'repeatedcv', 
                     repeats = 3, 
                     classProbs = T, 
                     summaryFunction = twoClassSummary)

# Make model on training set
set.seed(1234)
(fit <- train(outcome ~ ., data = training, 
              method = 'rf', 
              tuneLength = 7, 
              preProcess = c('center','scale'), 
              trControl = ctrl, 
              metric = 'ROC',
              ntree = 500))

# Predict data
new.pred <- predict(fit, newdata = testing) # gives prediction
new.prob <- predict(fit, newdata = testing, type = "prob") # gives probability of both
confusionMatrix(data = new.pred, testing$outcome, positive = 'Comp')

# Visualize results
varImp(rfFit)
pred <- prediction(rfProbs$comp, testing$overall_comp)
perf <- performance(pred, 'tpr','fpr')
plot(perf,col='red'); abline(0,1)
performance(pred, 'auc')@y.values
plot(performance(pred, 'acc'))


# OTU: Longitudinal - Complications -----------------------------------------------------------

# Add outcome to tax table
sub <- factor(dfu$v.comp[dfu$has_sample], labels = c('NoComp','Comp'))
sub <- cbind(outcome = sub, as.data.frame(tax))

# Make training and test sets
set.seed(1234)
inTrain <-  createDataPartition(y = sub$outcome, p = .75, list = F)
training <- sub[inTrain,]
testing <- sub[-inTrain,]

# Define training parameters
ctrl <- trainControl(method = 'repeatedcv', 
                     repeats = 3, 
                     classProbs = T, 
                     summaryFunction = twoClassSummary)

# Make model on training set
set.seed(1234)
(fit <- train(outcome ~ ., data = training, 
              method = 'rf', 
              #tuneLength = 7, 
              tuneGrid = data.frame(mtry = c(5,10,15,25,35,50)),
              preProcess = c('center','scale'), 
              trControl = ctrl, 
              metric = 'ROC',
              ntree = 500))

# Predict data
new.pred <- predict(fit, newdata = testing) # gives prediction
new.prob <- predict(fit, newdata = testing, type = "prob") # gives probability of both
confusionMatrix(data = new.pred, testing$outcome, positive = 'Comp')

# Visualize results
varImp(rfFit)
pred <- prediction(rfProbs$comp, testing$overall_comp)
perf <- performance(pred, 'tpr','fpr')
plot(perf,col='red'); abline(0,1)
performance(pred, 'auc')@y.values
plot(performance(pred, 'acc'))



# Genera: Baseline - Complications ------------------------------------------------------------

# Subset to samples at baseline
sub <- factor(dfu$comp[dfu$has_sample], labels = c('NoComp','Comp'))
sub <- cbind(outcome = sub, as.data.frame(taxa))
sub <- sub[dfu$visit[dfu$has_sample] == 0,]

# Make training and test sets
set.seed(1234)
inTrain <-  createDataPartition(y = sub$outcome, p = .75, list = F)
training <- sub[inTrain,]
testing <- sub[-inTrain,]

# Define training parameters
ctrl <- trainControl(method = 'repeatedcv', 
                     repeats = 3, 
                     classProbs = T, 
                     summaryFunction = twoClassSummary)

# Make model on training set
set.seed(1234)
(fit <- train(outcome ~ ., data = training, 
              method = 'rf', 
              tuneLength = 7, 
              preProcess = c('center','scale'), 
              trControl = ctrl, 
              metric = 'ROC',
              ntree = 500))

# Predict data
new.pred <- predict(fit, newdata = testing) # gives prediction
new.prob <- predict(fit, newdata = testing, type = "prob") # gives probability of both
confusionMatrix(data = new.pred, testing$outcome, positive = 'Comp')

# Visualize results
varImp(rfFit)
pred <- prediction(rfProbs$comp, testing$overall_comp)
perf <- performance(pred, 'tpr','fpr')
plot(perf,col='red'); abline(0,1)
performance(pred, 'auc')@y.values
plot(performance(pred, 'acc'))


# OTU: Longitudinal - Complications -----------------------------------------------------------

# Add outcome to tax table
sub <- factor(dfu$v.wd[dfu$has_sample], labels = c('NoComp','Comp'))
sub <- cbind(outcome = sub, as.data.frame(taxa))

# Make training and test sets
set.seed(1234)
inTrain <-  createDataPartition(y = sub$outcome, p = .75, list = F)
training <- sub[inTrain,]
testing <- sub[-inTrain,]

# Define training parameters
ctrl <- trainControl(method = 'repeatedcv', 
                     repeats = 3, 
                     classProbs = T, 
                     summaryFunction = twoClassSummary)

# Make model on training set
set.seed(1234)
(fit <- train(outcome ~ ., data = training, 
              method = 'rf', 
              #tuneLength = 7, 
              tuneGrid = data.frame(mtry = c(5,15,25,50,75,100)),
              preProcess = c('center','scale'), 
              trControl = ctrl, 
              metric = 'ROC',
              ntree = 500))

# Predict data
new.pred <- predict(fit, newdata = testing) # gives prediction
new.prob <- predict(fit, newdata = testing, type = "prob") # gives probability of both
confusionMatrix(data = new.pred, testing$outcome, positive = 'Comp')

# Visualize results
varImp(rfFit)
pred <- prediction(rfProbs$comp, testing$overall_comp)
perf <- performance(pred, 'tpr','fpr')
plot(perf,col='red'); abline(0,1)
performance(pred, 'auc')@y.values
plot(performance(pred, 'acc'))


##### Longitudinal Split by SubjectID #####

## Make training and test sets
# Make a random list of samples
set.seed(1234)
pt_id <- sample(levels(dfu.condensed$patient_id), replace = F)
# Add subject samples to a training set until it has 284 samples
i <- 1
train2.pt <- dfu.condensed[0,]
while(nrow(train2.pt) < 284) {
  train2.pt <- rbind(train2.pt, dfu.condensed[dfu.condensed$patient_id == pt_id[i],])
  i <- i + 1
}
# Make test set the remaining patients
test2.pt <- dfu.condensed[dfu.condensed$patient_id %in% pt_id[i:length(pt_id)],]

# Make model on training set
set.seed(1234)
fit.visit_b4_healing.pt <- randomForest(y = as.factor(train2.pt$is_visit_b4_healing), x = train2.pt[,-c(1:8)], ntree = 2000, importance = T)
fit.visit_b4_comp.pt <- randomForest(y = as.factor(train2.pt$visit_b4_complication), x = train2.pt[,-c(1:8)], ntree = 2000, importance = T)

# Visualize results
plot(fit.visit_b4_healing.pt)
varImpPlot(fit.visit_b4_healing.pt)
plot(fit.visit_b4_comp.pt)
varImpPlot(fit.visit_b4_comp.pt)

# Predict data
predict.visit_b4_healing.pt <- predict(fit.visit_b4_healing.pt, test2.pt[,-c(1:8)])
sum(predict.visit_b4_healing.pt == as.factor(test2.pt$is_visit_b4_healing))/nrow(test2.pt)
predict.visit_b4_comp.pt <- predict(fit.visit_b4_comp.pt, test2.pt[,-c(1:8)])
sum(predict.visit_b4_comp.pt == as.factor(test2.pt$visit_b4_complication))/nrow(test2.pt)



# Random GLM baseline basic analysis --------------------------------------------

# Make training and test sets
tmp <- t(apply(taxa.strep.0, 1, function(x) x/sum(x)))
set.seed(1234)
inTrain <- createDataPartition(y = dfu.0.condensed$overall_comp, p = .75, list = F)
training <- tmp[inTrain,]
testing <- tmp[-inTrain,]

# Make randomGLM
set.seed(1234)
rglm <- randomGLM(x = training, y = dfu.0.condensed$overall_comp[inTrain], testing, classify=TRUE, keepModels=TRUE, nThreads=2, randomSeed = 1234)

# Look at out of bag predictions
table(dfu.0.condensed$overall_comp[inTrain], rglm$predictedOOB)
# Look at test set
table(dfu.0.condensed$overall_comp[-inTrain], rglm$predictedTest)



