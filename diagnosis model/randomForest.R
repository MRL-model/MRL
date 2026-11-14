input=data.frame(flag=NA,gsva_matrix)
input$flag=sample_type[,2]
input$flag=as.factor(input$flag)
set.seed(234)
ind=sample(1:10,nrow(input),replace=TRUE)   
out=matrix(0,10,14)#########
colnames(out)=c("Accuracy","Kappa","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision","Recall","F1","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy","AUC")
for(i in 1:10){
  rf_try=randomForest(flag~.,data=input[ind!=i,],mtry=2,ntree=2000,nodesize=1,maxnodes=500,importance=TRUE)
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
