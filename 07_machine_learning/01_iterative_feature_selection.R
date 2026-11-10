
#h2omoedl_save=h2omoedl
ssave=h2omoedl_save[ooutre_sorted$`rownames(ooutre)`,]
rownames(ssave)=1:nrow(ssave)
h2omoedl_1=model[-which(model$variable==h2omoedl_save$variable),]
h2omoedl=rbind(ssave,h2omoedl_1)
h2omoedl$variable=gsub('-','.',h2omoedl$variable)

inner=exp
colnames(inner)=gsub('-','.',colnames(inner))

inner=cbind(rownames(exp),inner)


inner=merge(nmfg,inner,by = 'V1')

inner=inner[!duplicated(inner$V1),]

rownames(inner)=inner[,1]

inner=inner[,-1]

colnames(inner)[1]='status'

inner$status=as.factor(inner$status)

inner[,-1] <- apply(inner[,-1], 2, as.numeric)

set.seed(123)

ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))


# Assume model AUC vectors auc1 and auc2 are already available
# Initialize index to 1
index <- 1
# Start counting at 1
i <- 1

for (i in 1) {
  hmodel=h2omoedl[i,1]
  
  immer=inner[,c('status',hmodel)]
  immer[,1] <- factor(immer[,1])
  train <- immer[ind==1, ]
  test <- immer[ind==2, ]
  train1=train
  names(train1) <- make.names(names(train1))
  ##rf
  train <- immer[ind==1, ]
  test <- immer[ind==2, ]
  train1=train
  names(train1) <- make.names(names(train1))
  ran <- randomForest(status~.,data = train1,ntree=500,importance=TRUE,mtry=3,proximity=TRUE)
  test1=test
  names(test1) <- make.names(names(test1))
  pre_ran <- predict(ran,newdata=test1)
  pre_ran1=predict(ran,newdata=train1)
  obs_ran <- data.frame(prob=pre_ran,obs=test1$status)
  table(test1$status,pre_ran,dnn = c("Observed","Predicted"))
  ran_roc <- multiclass.roc(test1$status,as.numeric(pre_ran))
  ran_roc_tr=multiclass.roc(train1$status,as.numeric(pre_ran1))
  
  ran_roc
  ran_roc_tr
  
  ####boruta
  gl_m <- glm(status ~.,data = train1, family = binomial(link = "logit"))
  bb=predict.glm(gl_m,newdata = train1)
  pp=predict.glm(gl_m,newdata = test1)
  
  roc_test=multiclass.roc(test1$status,pp)
  roc_train=multiclass.roc(train1$status,bb)
  roc_train$auc
  roc_test$auc
  
  ###xgboost
  
  #Convert predictors to a matrix
  traindata1 <- data.matrix(train[,c(2:ncol(train))])
  
  #Convert to a sparse matrix using Matrix with sparse = TRUE
  traindata2 <- Matrix(traindata1,sparse = T)
  #Convert the response to numeric values
  traindata3 <- as.numeric(gsub('NMF_','',train[,1]))
  traindata3=traindata3-1
  #Combine predictors and response in a list
  traindata4 <- list(data=traindata2,label=traindata3)
  #Construct the xgb.DMatrix from the sparse matrix
  dtrain <- xgb.DMatrix(data=traindata4$data,label=traindata4$label)
  
  testdata1 <- data.matrix(test[,c(2:ncol(test))])
  
  testdata2 <- Matrix(testdata1,sparse = T)
  #Convert the response to numeric values
  testdata3 <- as.numeric(gsub('NMF_','',test[,1]))
  testdata3=testdata3-1
  #Combine predictors and response in a list
  testdata4 <- list(data=testdata2,label=testdata3)
  #Construct the xgb.DMatrix from the sparse matrix
  dtest <- xgb.DMatrix(data=testdata4$data,label=testdata4$label)
  #Confusion matrix
  
  xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=3, nround=25)
  
  pre_tr = round(predict(xgb,newdata = dtrain))
  pre_xgb = round(predict(xgb,newdata = dtest))
  
  table(testdata3,pre_xgb,dnn=c("true","pre"))
  
  #ROC curve
  xgboost_tr <- multiclass.roc(traindata3,as.numeric(pre_tr))
  xgboost_roc <- multiclass.roc(testdata3,as.numeric(pre_xgb))
  
  xgboost_tr$auc
  xgboost_roc$auc
  mean=(ran_roc$auc+roc_test$auc+xgboost_roc$auc)/3 
}

mean
# Iteration condition; adjust as needed
ppppp=matrix(nrow = nrow(h2omoedl),ncol = 4)

#i=1185
i=2
while (i <=  nrow(h2omoedl)) {

  mmudex=index
  hmodel=h2omoedl[c(mmudex,i),1]
  
  immer=inner[,c('status',hmodel)]
  immer[,1] <- factor(immer[,1])
  train <- immer[ind==1, ]
  test <- immer[ind==2, ]
  train1=train
  names(train1) <- make.names(names(train1))
  
  train <- immer[ind==1, ]
  test <- immer[ind==2, ]
  train1=train
  
  ran <- randomForest(status~.,data = train1,ntree=500,importance=TRUE,mtry=3,proximity=TRUE)
  test1=test
  names(test1) <- make.names(names(test1))
  pre_ran <- predict(ran,newdata=test1)
  pre_ran1=predict(ran,newdata=train1)
  obs_ran <- data.frame(prob=pre_ran,obs=test1$status)
  table(test1$status,pre_ran,dnn = c("Observed","Predicted"))
  ran_roc_tr=multiclass.roc(as.numeric(train1$status),as.numeric(pre_ran1))
  ran_roc <- multiclass.roc(as.numeric(test1$status),as.numeric(pre_ran))
  ran_roc
  ran_roc_tr
  
  ####boruta
  gl_m <- glm(status ~.,data = train1, family = binomial(link = "logit"))
  bb=predict.glm(gl_m,newdata = train1)
  pp=predict.glm(gl_m,newdata = test1)
  
  roc_test=multiclass.roc(test1$status,pp)
  roc_train=multiclass.roc(train1$status,bb)
  roc_train$auc
  roc_test$auc
  
  ###xgboost
  
  #Convert predictors to a matrix
  traindata1 <- data.matrix(train[,c(2:ncol(train))])
  
  #Convert to a sparse matrix using Matrix with sparse = TRUE
  traindata2 <- Matrix(traindata1,sparse = T)
  #Convert the response to numeric values
  traindata3 <- as.numeric(gsub('NMF_','',train[,1]))
  traindata3=traindata3-1
  #Combine predictors and response in a list
  traindata4 <- list(data=traindata2,label=traindata3)
  #Construct the xgb.DMatrix from the sparse matrix
  dtrain <- xgb.DMatrix(data=traindata4$data,label=traindata4$label)
  
  testdata1 <- data.matrix(test[,c(2:ncol(test))])
  
  testdata2 <- Matrix(testdata1,sparse = T)
  #Convert the response to numeric values
  testdata3 <- as.numeric(gsub('NMF_','',test[,1]))
  testdata3=testdata3-1
  #Combine predictors and response in a list
  testdata4 <- list(data=testdata2,label=testdata3)
  #Construct the xgb.DMatrix from the sparse matrix
  dtest <- xgb.DMatrix(data=testdata4$data,label=testdata4$label)
  #Confusion matrix
  
  xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=3, nround=25)
  
  pre_tr = round(predict(xgb,newdata = dtrain))
  pre_xgb = round(predict(xgb,newdata = dtest))
  
  table(testdata3,pre_xgb,dnn=c("true","pre"))
  
  #ROC curve
  xgboost_tr <- multiclass.roc(traindata3,as.numeric(pre_tr))
  xgboost_roc <- multiclass.roc(testdata3,as.numeric(pre_xgb))
  
  xgboost_tr$auc
  xgboost_roc$auc

  # Compute the mean AUC across the three models
  mean_auc=(ran_roc$auc+roc_test$auc+xgboost_roc$auc)/3
  # Retain candidate i if mean AUC improves and each model AUC exceeds 0.5
  if (mean_auc > mean && ran_roc$auc > 0.5 && roc_test$auc > 0.5 && xgboost_roc$auc > 0.5) {
    index <- c(index, i)
    mean=mean_auc
  }
  # Otherwise retain the current index
  else {
    index <- index
  }
  # Record the current iteration, selected indices, and AUC values
  ppppp[i,]=c(i,paste(index, collapse = " "),mean,mean_auc)
  # Increment the counter
  i <- i + 1
  
}

View(p)
finalmodel=h2omoedl[index,1]
finalmodel
{
  
  hmodel=finalmodel
  
  immer=inner[,c('status',hmodel)]
  immer[,1] <- factor(immer[,1])
  train <- immer[ind==1, ]
  test <- immer[ind==2, ]
  train1=train
  names(train1) <- make.names(names(train1))
  
  train <- immer[ind==1, ]
  test <- immer[ind==2, ]
  train1=train
  
  ran <- randomForest(status~.,data = train1,ntree=500,importance=TRUE,mtry=3,proximity=TRUE)
  ran1 <- ran$importance
  ran1 <- as.data.frame(ran1)
  ran2 <- ran1[order(ran1$MeanDecreaseGini,decreasing = T),]
  ran3 <- subset(ran2,ran2$MeanDecreaseGini>0.1)
  test1=test
  names(test1) <- make.names(names(test1))
  pre_ran <- predict(ran,newdata=test1)
  pre_ran1=predict(ran,newdata=train1)
  obs_ran <- data.frame(prob=pre_ran,obs=test1$status)
  table(test1$status,pre_ran,dnn = c("Observed","Predicted"))
  ran_roc <- roc(test1$status,as.numeric(pre_ran))
  ran_roc_tr=roc(train1$status,as.numeric(pre_ran1))
  
  ran_roc
  ran_roc_tr
  
  ####boruta
  gl_m <- glm(status ~.,data = train1, family = binomial(link = "logit"))
  bb=predict.glm(gl_m,newdata = train1)
  pp=predict.glm(gl_m,newdata = test1)
  
  roc_test=roc(test1$status,pp)
  roc_train=roc(train1$status,bb)
  roc_train$auc
  roc_test$auc
  
  ###xgboost
  
  #Convert predictors to a matrix
  traindata1 <- data.matrix(train[,c(2:ncol(train))])
  
  #Convert to a sparse matrix using Matrix with sparse = TRUE
  traindata2 <- Matrix(traindata1,sparse = T)
  #Convert the response to numeric values
  traindata3 <- as.numeric(gsub('NMF_','',train[,1]))
  traindata3=traindata3-1
  #Combine predictors and response in a list
  traindata4 <- list(data=traindata2,label=traindata3)
  #Construct the xgb.DMatrix from the sparse matrix
  dtrain <- xgb.DMatrix(data=traindata4$data,label=traindata4$label)
  
  testdata1 <- data.matrix(test[,c(2:ncol(test))])
  
  testdata2 <- Matrix(testdata1,sparse = T)
  #Convert the response to numeric values
  testdata3 <- as.numeric(gsub('NMF_','',test[,1]))
  testdata3=testdata3-1
  #Combine predictors and response in a list
  testdata4 <- list(data=testdata2,label=testdata3)
  #Construct the xgb.DMatrix from the sparse matrix
  dtest <- xgb.DMatrix(data=testdata4$data,label=testdata4$label)
  #Confusion matrix
  
  xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=3, nround=25)
  
  pre_tr = round(predict(xgb,newdata = dtrain))
  pre_xgb = round(predict(xgb,newdata = dtest))
  
  table(testdata3,pre_xgb,dnn=c("true","pre"))
  
  #ROC curve
  xgboost_tr <- roc(traindata3,as.numeric(pre_tr))
  xgboost_roc <- roc(testdata3,as.numeric(pre_xgb))
  
  xgboost_tr$auc
  xgboost_roc$auc
  
  # Compute the mean AUC across the three models
  mean_auc=(ran_roc$auc+roc_test$auc+xgboost_roc$auc)/3
  
  
  }



save.image(file = 'chris18625model.rds')
