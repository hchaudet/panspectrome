library(this.path)
Dir.source <- this.dir()
Dir.base <- paste0(rev(rev(unlist(strsplit(Dir.source,.Platform$file.sep)))[-1]),collapse=.Platform$file.sep)
Dir.data <- file.path(Dir.source,'data')

setwd(Dir.source)

library(dynamicTreeCut)
library(dendextend)
library(moduleColor)

# Construction de la matrice des distances = freq relative des co-occurences
adjacencyToDist <- function(adjacency,nspectra) {
  dist <- adjacency
  dist <- dist/nspectra
  dist <- 1-dist
  diag(dist) <- 0
  dist <- as.dist(dist)
  return(dist)
}

hclustplotn <- function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...) 
  #hier1 <- hc
  #Color <- DetectedColors
{
  options(stringsAsFactors=FALSE)
  
  #No.Sets = dim(Color)[[2]];
  C <- Color[hier1$order]
  step <- 1/length(Color)
  #ystep = 1/No.Sets;
  barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
  ind <- seq_along(C)
  xl <- (ind-1) * step
  xr <- ind * step 
  yb <- rep(0, length(C))
  yt <- rep(1, length(C))
  rect(xl, yb, xr, yt, col = as.character(C), border = as.character(C))
  #text('Clusters', pos=2, x=0, y=-0.5, cex=cex.RowLabels, xpd = TRUE)
  lines(x=c(0,1), y=c(1,1))
}

dynamicCutPlot <- function(dynamicCut,hc) {
  DetectedColors <- labels2colors(dynamicCut$labels)
  layout(matrix(c(1,2), nrow = 2, ncol = 1), widths = 1, heights = c(0.975, 0.025));
  par(cex = 0.3)
  par(mar=c(0, 6.5, 2, 0.2))
  plot(hc)
  par(mar=c(0.2, 6.5, 0, 0.2))
  
  hclustplotn(hc,
              DetectedColors,
              main="");
  
}

#labels2colors(1:15)


# load workspace 'analyse_clusters.RData'

peak.freq <- colSums(xfreq)
peak.stat <- data.frame(n=peak.freq,f=peak.freq/nrow(xfreq))
peak.stat <- cbind(peak=rownames(peak.stat),peak.stat)


dist <- adjacencyToDist(adjacency,nrow(xfreq))
hc <- hclust(dist,method="mcquitty")


## 4.Visualisation de la clusterisation

plot(hc,cex=0.3,axes=FALSE,ylab='Relative frequency',xlab='Peaks')
ticks <- seq(0,1,.1)
axis(2,at=ticks,labels=1-ticks)


###
dynamicCut <- cutreeHybrid(hc, minClusterSize=5, distM=as.matrix(dist), cutHeight = .99, maxCoreScatter = .99, minGap = .007, verbose=0)
clusters <- as.data.frame(cbind(peak=as.numeric(hc$labels),labels=dynamicCut$labels,cores=dynamicCut$cores,smallLabels=dynamicCut$smallLabels,onBranch=dynamicCut$onBranch))
clusters <- clusters[clusters[,2]>0,]
clusters <- merge(clusters,peak.stat,by='peak',all.y=FALSE)
clusters <- clusters[order(clusters[,2],clusters[,3]),]

dynamicCutPlot(dynamicCut,hc)


####################################
####################################
####################################

# Etude des carbapénémases
########


intens.carba.pos <- intensities[data$spectrum_id[which(data$carba=='POS')],]
xfreq.carba.pos <- xfreq[data$spectrum_id[which(data$carba=='POS')],]
intens.carba.neg <- intensities[data$spectrum_id[which(data$carba=='NEG')],]
xfreq.carba.neg <- xfreq[data$spectrum_id[which(data$carba=='NEG')],]

peak.freq.carba.pos <- colSums(xfreq.carba.pos)
peak.stat.carba.pos <- data.frame(n=peak.freq.carba.pos,f=peak.freq.carba.pos/nrow(xfreq.carba.pos))
peak.stat.carba.pos <- cbind(peak=rownames(peak.stat.carba.pos),peak.stat.carba.pos)
peak.freq.carba.neg <- colSums(xfreq.carba.neg)
peak.stat.carba.neg <- data.frame(n=peak.freq.carba.neg,f=peak.freq.carba.neg/nrow(xfreq.carba.neg))
peak.stat.carba.neg <- cbind(peak=rownames(peak.stat.carba.neg),peak.stat.carba.neg)


#sac_curve(xfreq.carba.neg)

# Construction de la matrice d'adjacence
adjacency.carba.pos <- t(xfreq.carba.pos)%*%xfreq.carba.pos
diag(adjacency.carba.pos) <- 0

dist.carba.pos <- adjacencyToDist(adjacency.carba.pos,nrow(xfreq.carba.pos))
hc.carba.pos <- hclust(dist.carba.pos,method="mcquitty")

adjacency.carba.neg <- t(xfreq.carba.neg)%*%xfreq.carba.neg
diag(adjacency.carba.neg) <- 0

dist.carba.neg <- adjacencyToDist(adjacency.carba.neg,nrow(xfreq.carba.neg))
hc.carba.neg <- hclust(dist.carba.neg,method="mcquitty")



## 4.Visualisation de la clusterisation
plot(hc.carba.pos,cex=0.3,axes=FALSE,main='Carba pos', ylab='Relative frequency',xlab='Peaks')
ticks <- seq(0,1,.1)
axis(2,at=ticks,labels=1-ticks)

plot(hc.carba.neg,cex=0.3,axes=FALSE,main='Carba neg',ylab='Relative frequency',xlab='Peaks')
ticks <- seq(0,1,.1)
axis(2,at=ticks,labels=1-ticks)





dynamicCut <- cutreeHybrid(hc.carba.pos, minClusterSize=5, distM=as.matrix(dist.carba.pos), cutHeight = .99, maxCoreScatter = .99, minGap = .007, verbose=0)
#dynamicCut <- cutreeDynamicTree(hc.carba.pos, minClusterSize=5, cutHeight = .99, maxCoreScatter = .99, minGap = .007, verbose=0)
clusters.carba.pos <- as.data.frame(cbind(peak=as.numeric(hc.carba.pos$labels),
                                          labels=dynamicCut$labels,
                                          cores=dynamicCut$cores,
                                          smallLabels=dynamicCut$smallLabels,
                                          onBranch=dynamicCut$onBranch))
clusters.carba.pos <- clusters.carba.pos[clusters.carba.pos[,2]>0,]
clusters.carba.pos <- merge(clusters.carba.pos,peak.stat.carba.pos,by='peak',all.y=FALSE)
clusters.carba.pos <- clusters.carba.pos[order(clusters.carba.pos[,2],clusters.carba.pos[,3]),]

dynamicCutPlot(dynamicCut,hc.carba.pos)




dynamicCut <- cutreeHybrid(hc.carba.neg, minClusterSize=5, distM=as.matrix(dist.carba.neg), cutHeight = .99, maxCoreScatter = .99, minGap = .007, verbose=0)
clusters.carba.neg <- as.data.frame(cbind(peak=as.numeric(hc.carba.pos$labels),
                                          labels=dynamicCut$labels,
                                          cores=dynamicCut$cores,
                                          smallLabels=dynamicCut$smallLabels,
                                          onBranch=dynamicCut$onBranch))
clusters.carba.neg <- clusters.carba.neg[clusters.carba.neg[,2]>0,]
clusters.carba.neg <- merge(clusters.carba.neg,peak.stat.carba.neg,by='peak',all.y=FALSE)
clusters.carba.neg <- clusters.carba.neg[order(clusters.carba.neg[,2],clusters.carba.neg[,3]),]

dynamicCutPlot(dynamicCut,hc.carba.neg)




## Listes des disjonctions
clusters.carba.pos[which(clusters.carba.pos$peak %in% setdiff(clusters.carba.pos$peak,clusters.carba.neg$peak)),]
clusters.carba.neg[which(clusters.carba.neg$peak %in% setdiff(clusters.carba.neg$peak,clusters.carba.pos$peak)),]

clusters.carba <- merge(clusters.carba.pos[,-(4:5)],clusters.carba.neg[,-(4:5)],by='peak',all=TRUE)
colnames(clusters.carba) <- c('peak','pos.cluster','pos.core','pos.n','pos.f','neg.cluster','neg.core','neg.n','neg.f')


library(circlize)
clusters.carba.flow <- data.frame(clusters.carba$pos.cluster,clusters.carba$neg.cluster)
clusters.carba.flow[is.na(clusters.carba.flow)] <- 0
colnames(clusters.carba.flow) <- c('pos','neg')
d <- matrix(table(clusters.carba.flow),nrow=length(unique(clusters.carba.flow$pos)),ncol=length(unique(clusters.carba.flow$neg)))
rownames(d) <- paste0("pos-", seq(0,nrow(d)-1))
colnames(d) <- paste0("neg-", seq(0,ncol(d)-1))
par(cex = 1, mar = c(0, 0, 0, 0))
chordDiagram(d, transparency = 0.5)


####################################
####################################
####################################

## Random Forest
library(randomForest)
set.seed(123)

## 1. Tous pics

intens.work <- intensities[data$spectrum_id[which(data$carba %in% c('POS','NEG'))],]
intens.work <- intens.work[,which(colSums(intens.work)>0)]
xfreq.work <- intens.work
xfreq.work[xfreq.work>0] <- 1

work <- data.frame(carba=data$carba[which(data$carba %in% c('POS','NEG'))],xfreq.work)
work$carba <- as.factor(work$carba)

fit <- randomForest(carba ~ ., data = work)
print(fit)
varImpPlot(fit)

work <- data.frame(carba=data$carba[which(data$carba %in% c('POS','NEG'))],intens.work)
work$carba <- as.factor(work$carba)

fit <- randomForest(carba ~ ., data = work)
print(fit)
varImpPlot(fit)




## 2. Pics de clusters

intens.work <- intensities[data$spectrum_id[which(data$carba %in% c('POS','NEG'))],]
peak.list <- as.character(sort(unique(c(clusters.carba.pos$peak,clusters.carba.neg$peak))))
intens.work <- intens.work[,which(colnames(intens.work) %in% peak.list)]
xfreq.work <- intens.work
xfreq.work[xfreq.work>0] <- 1

work <- data.frame(carba=data$carba[which(data$carba %in% c('POS','NEG'))],xfreq.work)
work$carba <- as.factor(work$carba)

fit <- randomForest(carba ~ ., data = work)
print(fit)
varImpPlot(fit)

work <- data.frame(carba=data$carba[which(data$carba %in% c('POS','NEG'))],intens.work)
work$carba <- as.factor(work$carba)

fit <- randomForest(carba ~ ., data = work)
print(fit)
varImpPlot(fit)


#plot(fit$err.rate[, 1], type = "l", xlab = "nombre d'arbres", ylab = "erreur OOB")

#fit <- randomForest(carba ~ ., data = work,ntree=2000)



library(C50)

labels <- data$carba[which(data$carba %in% c('POS','NEG'))]
labels <- as.factor(labels)

treeModel <- C5.0(x = intens.work, y = labels)
treeModel
summary(treeModel)



library(xgboost)

labels <- data$carba[which(data$carba %in% c('POS','NEG'))]
labels <- ifelse(labels=="POS",1,0)

xgb_train = xgb.DMatrix(data = intens.work, label = labels)
xgb_test = xgb.DMatrix(data = intens.work, label = labels)
watchlist = list(train=xgb_train, test=xgb_test)

model = xgboost(data = xgb_train, max.depth = 3, nrounds = 145, objective = "binary:logistic")
pred <- predict(model, xgb_test)
val.pred <- as.numeric(pred > 0.5)
err <- mean(val.pred != labels)
print(paste("test-error=", err))

caret::confusionMatrix(factor(labels), factor(val.pred))



library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

params_booster <- list(booster = 'gbtree', 
                       #eta = .3, 
                       #gamma = 0, 
                       max.depth = 3, 
                       #subsample = 1, 
                       #colsample_bytree = 1, 
                       #min_child_weight = 1, 
                       objective = "binary:logistic")

bst.cv <- xgb.cv(data = xgb_train, 
                 label = labels, 
                 params = params_booster,
                 nrounds = 200, 
                 nfold = 5,
                 print_every_n = 20,
                 verbose = 2)

data.frame(TRAINING_ERROR = bst.cv$evaluation_log$train_logloss_mean, 
           VALIDATION_ERROR = bst.cv$evaluation_log$test_logloss_mean,
           ITERATION = bst.cv$evaluation_log$iter)


res_df <- data.frame(TRAINING_ERROR = bst.cv$evaluation_log$train_logloss_mean, 
                     VALIDATION_ERROR = bst.cv$evaluation_log$test_logloss_mean,
                     ITERATION = bst.cv$evaluation_log$iter) %>%
  mutate(MIN = VALIDATION_ERROR == min(VALIDATION_ERROR))

best_nrounds <- res_df %>%
  filter(MIN) %>%
  pull(ITERATION)

res_df_longer <- pivot_longer(data = res_df, 
                              cols = c(TRAINING_ERROR, VALIDATION_ERROR), 
                              names_to = "ERROR_TYPE",
                              values_to = "ERROR")

ggplot(res_df_longer, aes(x = ITERATION)) +        # Look @ it overfit.
  geom_line(aes(y = ERROR, group = ERROR_TYPE, colour = ERROR_TYPE)) +
  geom_vline(xintercept = best_nrounds, colour = "green") +
  geom_label(aes(label = str_interp("${best_nrounds} iterations gives minimum validation error"), y = 0.2, x = best_nrounds, hjust = 0.1)) +
  labs(
    x = "nrounds",
    y = "Error",
    title = "Test & Train Errors",
    subtitle = str_interp("Note how the training error keeps decreasing after ${best_nrounds} iterations, but the validation error starts \ncreeping up. This is a sign of overfitting.")
  ) +
  scale_colour_discrete("Error Type: ")

model = xgboost(data = xgb_train, max.depth = 3, nrounds = 145, objective = "binary:logistic")
pred <- predict(model, xgb_test)
val.pred <- as.numeric(pred > 0.5)
err <- mean(val.pred != labels)
print(paste("test-error=", err))

caret::confusionMatrix(factor(labels), factor(val.pred))



## XGBoost

library(tidyverse)
library(xgboost)
library(caret)
library(readxl)

set.seed(1337)  # Pour la 'reproductibilité'

labels <- data$carba[which(data$carba %in% c('POS','NEG'))]
labels <- ifelse(labels=="POS",0,1)
inTrain <- createDataPartition(y = labels, p = 0.85, list = FALSE)  # 85% des données dans le train, et le rest dans le test 

### Tous
intens.work <- intensities[data$spectrum_id[which(data$carba %in% c('POS','NEG'))],]
xfreq.work <- intens.work
xfreq.work[xfreq.work>0] <- 1

dim(intens.work)



training <- intens.work[inTrain, ]
testing <- intens.work[-inTrain, ]

training <- xfreq.work[inTrain, ]
testing <- xfreq.work[-inTrain, ]

y_train <- labels[inTrain]
X_train = as.matrix(training)
y_test <- labels[-inTrain]
X_test = as.matrix(testing)

xgb_trcontrol <- trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                              verboseIter = FALSE, returnData = FALSE)
xgbGrid <- expand.grid(nrounds = c(50,100,150,200,250),  
                       max_depth = c(3, 5, 10, 15, 20),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       ## valeurs par défaut : 
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1
)

xgb_model <- train(X_train, y_train, trControl = xgb_trcontrol, tuneGrid = xgbGrid, 
                   method = "xgbTree", objective = "binary:logistic", eval_metric = "auc")

xgb_model$bestTune

#xgb_model$bestTune
#nrounds max_depth eta gamma colsample_bytree min_child_weight subsample
#108     150        20 0.1     0              0.6                1         1

pred <- predict(xgb_model, X_test)
val.pred <- as.numeric(pred > 0.5)
err <- mean(val.pred != y_test)
print(paste("test-error=", err))

caret::confusionMatrix(factor(y_test), factor(val.pred))


xgb.plot.importance(importance_matrix = xgb.importance(colnames(X_train), xgb_model$finalModel))



### Clusters
intens.work <- intensities[data$spectrum_id[which(data$carba %in% c('POS','NEG'))],]
peak.list <- as.character(sort(unique(c(clusters.carba.pos$peak,clusters.carba.neg$peak))))
intens.work <- intens.work[,which(colnames(intens.work) %in% peak.list)]
xfreq.work <- intens.work
xfreq.work[xfreq.work>0] <- 1

labels <- data$carba[which(data$carba %in% c('POS','NEG'))]
labels <- ifelse(labels=="POS",0,1)

dim(intens.work)


training <- intens.work[inTrain, ]
testing <- intens.work[-inTrain, ]

training <- xfreq.work[inTrain, ]
testing <- xfreq.work[-inTrain, ]

y_train <- labels[inTrain]
X_train = as.matrix(training)
y_test <- labels[-inTrain]
X_test = as.matrix(testing)

xgb_trcontrol <- trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                              verboseIter = FALSE, returnData = FALSE)
xgbGrid <- expand.grid(nrounds = c(50,100,150,200,250),  
                       max_depth = c(3, 5, 10, 15, 20),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       ## valeurs par défaut : 
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1
)

xgb_model <- train(X_train, y_train, trControl = xgb_trcontrol, tuneGrid = xgbGrid, 
                   method = "xgbTree", objective = "binary:logistic", eval_metric = "auc")

xgb_model$bestTune


pred <- predict(xgb_model, X_test)
val.pred <- as.numeric(pred > 0.5)
err <- mean(val.pred != y_test)
print(paste("test-error=", err))

caret::confusionMatrix(factor(y_test), factor(val.pred))


xgb.plot.importance(importance_matrix = xgb.importance(colnames(X_train), xgb_model$finalModel))


### Hypothèse : virer le core cluster principal
intens.work <- intensities[data$spectrum_id[which(data$carba %in% c('POS','NEG'))],]
peak.list <- as.character(sort(unique(c(clusters.carba.pos$peak[clusters.carba.pos$cores!=1],
                                        clusters.carba.neg$peak[clusters.carba.neg$cores!=1]))))
intens.work <- intens.work[,which(colnames(intens.work) %in% peak.list)]
xfreq.work <- intens.work
xfreq.work[xfreq.work>0] <- 1

labels <- data$carba[which(data$carba %in% c('POS','NEG'))]
labels <- ifelse(labels=="POS",0,1)

dim(intens.work)

set.seed(1337)  # Pour la 'reproductibilité'

inTrain <- createDataPartition(y = labels, p = 0.85, list = FALSE)  # 85% des données dans le train, et le rest dans le test 

training <- intens.work[inTrain, ]
testing <- intens.work[-inTrain, ]

training <- xfreq.work[inTrain, ]
testing <- xfreq.work[-inTrain, ]

y_train <- labels[inTrain]
X_train = as.matrix(training)
y_test <- labels[-inTrain]
X_test = as.matrix(testing)

xgb_trcontrol <- trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                              verboseIter = FALSE, returnData = FALSE)
xgbGrid <- expand.grid(nrounds = c(50,100,150,200,250),  
                       max_depth = c(3, 5, 10, 15, 20),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       ## valeurs par défaut : 
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1
)

xgb_model <- train(X_train, y_train, trControl = xgb_trcontrol, tuneGrid = xgbGrid, 
                   method = "xgbTree", objective = "binary:logistic", eval_metric = "auc")

xgb_model$bestTune


pred <- predict(xgb_model, X_test)
val.pred <- as.numeric(pred > 0.5)
err <- mean(val.pred != y_test)
print(paste("test-error=", err))

caret::confusionMatrix(factor(y_test), factor(val.pred))


xgb.plot.importance(importance_matrix = xgb.importance(colnames(X_train), xgb_model$finalModel))

##########################
# Equilibrage des effectifs

x.neg <- sample(1:nrow(intens.carba.neg),nrow(intens.carba.pos))
intens.work <- rbind(intens.carba.pos,intens.carba.neg[x.neg,])
peak.list <- as.character(sort(unique(c(clusters.carba.pos$peak,clusters.carba.neg$peak))))
intens.work <- intens.work[,which(colnames(intens.work) %in% peak.list)]
xfreq.work <- intens.work
xfreq.work[xfreq.work>0] <- 1
dim(intens.work)

labels <- c(rep('POS',nrow(intens.carba.pos)),rep('NEG',nrow(intens.carba.pos)))
labels <- ifelse(labels=="POS",0,1)

inTrain <- createDataPartition(y = labels, p = 0.85, list = FALSE)  # 80% des données dans le train, et le rest dans le test 


training <- intens.work[inTrain, ]
testing <- intens.work[-inTrain, ]

training <- xfreq.work[inTrain, ]
testing <- xfreq.work[-inTrain, ]

y_train <- labels[inTrain]
X_train = as.matrix(training)
y_test <- labels[-inTrain]
X_test = as.matrix(testing)

xgb_trcontrol <- trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                              verboseIter = FALSE, returnData = FALSE)
xgbGrid <- expand.grid(nrounds = c(50,100,150,200,250),  
                       max_depth = c(3, 5, 10, 15, 20),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       ## valeurs par défaut : 
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1
)

xgb_model <- train(X_train, y_train, trControl = xgb_trcontrol, tuneGrid = xgbGrid, 
                   method = "xgbTree", objective = "binary:logistic", eval_metric = "auc")

xgb_model$bestTune


pred <- predict(xgb_model, X_test)
val.pred <- as.numeric(pred > 0.5)
err <- mean(val.pred != y_test)
print(paste("test-error=", err))

caret::confusionMatrix(factor(y_test), factor(val.pred))


xgb.plot.importance(importance_matrix = xgb.importance(colnames(X_train), xgb_model$finalModel))


