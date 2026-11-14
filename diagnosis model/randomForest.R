library(GSVA)
library(GSEABase)
library(ggplot2)
library(simplevis)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(ggprism)
gsva_matrix=gsva(as.matrix(express),
                   list, 
                   method='ssgsea',
                   ssgsea.norm = TRUE,
                   kcdf='Gaussian',  #Character string denoting the kernel to use during the non-parametric estimation of the cumulative distribution function of expression levels across samples when method="gsva". By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson".
                   abs.ranking=TRUE) #When abs.ranking=TRUE the original Kuiper statistic that sums the largest positive and negative random walk deviations, is used. In this latter case, gene sets with genes enriched on either extreme (high or low) will be regarded as 'highly' activated.
input=data.frame(flag=NA,gsva_matrix)
input$flag=sample_type[,2]
input$flag=as.factor(input$flag)
set.seed(234)
ind=sample(1:10,nrow(input),replace=TRUE)   
out=matrix(0,10,14)
colnames(out)=c("Accuracy","Kappa","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision","Recall","F1","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy","AUC")
for(i in 1:10){
  rf_try=randomForest(flag~.,data=input[ind!=i,],mtry=mtry,ntree=ntree,nodesize=nodesize,maxnodes=maxnodes,importance=TRUE)
  fold_pred=predict(rf_try,input[ind==i,])
  result=confusionMatrix(factor(as.vector(fold_pred)),as.factor(input[ind==i,1]),mode="everything",positive='normal')
  result
  out[i,1]=result$overall["Accuracy"]
  out[i,2]=result$overall["Kappa"]
  out[i,3:13]=as.numeric(result$byClass)
  predict_pred=predict(rf_try,input[ind==i,],type='prob')
  pred=prediction(predict_pred[,'tumor'],input[ind==i,1])
  auc_try=performance(pred,'auc')
  out[i,14]=auc_try@y.values[[1]][1]             
}
