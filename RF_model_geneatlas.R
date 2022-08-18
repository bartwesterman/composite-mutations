###########################################
### PREDICTION MODEL RECURRENT MUTATIONS###
### random forest model for the top 4   ###
### genes that show the most recurrent  ###
### mutations. Both for chromosomal loc ###
### and for the mutation (0 - 1 ) of co ###
### -mutations. Combine models with     ###
### bagging                             ###
###########################################

##installing packages
# install.packages("caret")
# install.packages("randomForestSRC")

##load needed libraries
library(caret) #Splitting data into training and test set
library(randomForestSRC) #Package for Random forest
library(ggplot2) #Package for plots
library(ROCR) #package for ROC curve

##load input data from input_RF_model_geneatlas.R
##input chr data for the RF model

# in_BRAFc 
# in_EGFRc 
# in_KRASc 
# in_PIK3CAc 

##input mut data for the RF model
# in_BRAFm
# in_EGFRm 
# in_KRASm
# in_PIK3CAm

##function to run the RF model as an input you need the inputfile in_XXX and the interested gene name
##as an output you get two lists
##1. pred.conv the predicted values from the model
##2. testData which is the true labels of the testData
RF_model = function(input, gene){
  
  dataset = input[which(input$object %in% c("0", gene)),]
  dataset$object = as.factor(dataset$object)
  ### Downsample
  set.seed(85)
  total.conv = sum(table(dataset$object[dataset$object != 0])) #all cases with recurrent mutations
  nonconv = dataset[dataset$object == 0,] #all cases with no recurrent mutations
  set.seed(85)
  down_samplec = nonconv[sample(nrow(nonconv), total.conv),]
  convergentc = dataset[dataset$object != 0,]
  datasetT = rbind(down_samplec, convergentc)
  
  ## Split data
  set.seed(85)
  inTrainc = createDataPartition(datasetT$object, p = 2/3, list = FALSE)#split data according to label, so each label is been split into 2/3 training 
  trainData=datasetT[inTrainc,]
  testData=datasetT[-inTrainc,]
  ###Run RF model
  set.seed(85)
  model <- rfsrc(object ~ .,  data = trainData[,-which(names(trainData) %in% "Sample")],
                 #case.wt = randomForestSRC:::make.wt(trainData$conv),
                 #sampsize = randomForestSRC:::make.size(y))
                 var.used="all.trees", ntree=100,importance="none")

  
  #Testset
  pred.conv = predict(model, testData[,-which(names(testData) %in% "Sample")])
  output<-list(pred.conv,testData)
  return(output)
}

#Save output of model for co-mutation and chromosmal
out_B  = RF_model(in_BRAFm, "BRAF")
out_Bc = RF_model(in_BRAFc, "BRAF")

out_P  = RF_model(in_PIK3CAm, "PIK3CA")
out_Pc = RF_model(in_PIK3CAc, "PIK3CA")

out_E  = RF_model(in_EGFRm, "EGFR")
out_Ec = RF_model(in_EGFRc, "EGFR")

out_K  = RF_model(in_KRASm, "KRAS")
out_Kc = RF_model(in_KRASc, "KRAS")


###ENSEMBLE method Random forests
### https://www.analyticsvidhya.com/blog/2017/02/introduction-to-ensembling-along-with-implementation-in-r/
##training performance

##function to combine the co-mutation and chr model together. As an input you need
## 1. result of chr model
## 2. result of co-mut model
## 3. the testData 
## all derived from function RF_model
## as an OUTPUT you will get a list with the following results
## 1. cmfinal the confusion matrix of the ensemble model
## 2. roc.t the result for the ROC curve

ensemble_model = function(input_co,input_chr, testData, gene){
  chr_prob = as.data.frame(input_chr$predicted)
  co_prob = as.data.frame(input_co$predicted)
  truelabel = as.data.frame(testData$object)
  ##WEIGHTED VOTE
  #Taking weighted average of predictions
  prob.pred = ((co_prob*0.75)+(chr_prob*0.25))/2
  prob.pred$label = colnames(prob.pred)[max.col(prob.pred,ties.method="first")]
  prob.pred$true = testData$object
  
  #confusion matrix
  cm_weight = confusionMatrix(factor(prob.pred$label), truelabel$`testData$object`)
  cmT = cm_weight$table
  
  cmfinal = as.data.frame(cm_weight$table)
  cmfinal$Freq[cmfinal$Freq == 0] = 0.0001

  TN = sum(cmT) - ((rowSums(cmT) + colSums(cmT)) - diag(cmT))
  FN = rowSums(cmT) - diag(cmT)
  FP = colSums(cmT) - diag(cmT)
  TP = diag(cmT)

  accuracy <- sum(diag(cmT)) / sum(cmT)
  Baccuracy = (TP / (TP+FN) + TN/ (TN+FP)) /2
  precision = diag(cmT) / colSums(cmT)
  recall = diag(cmT) / rowSums(cmT)

  print("accuracy")
  print(accuracy)
  print("")
  print("balanced accuracy")
  print(Baccuracy)
  print("")
  print("precision")
  print(precision)
  print("")
  print("recall")
  print(recall)
  print("")
  
  #ROC curve
  gene.prob = prob.pred[,c(gene, "label", "true")]
  gene.prob$p_pred = as.factor(ifelse(gene.prob$label == gene, 1, 0))
  gene.prob$actual = as.factor(ifelse(gene.prob$true == gene, 1, 0))
  print(gene.prob)

  pred.gene = prediction(as.vector(as.numeric(gene.prob[,1])), as.vector(as.numeric(gene.prob$actual)))
  rfperf.gene = performance(pred.gene, "tpr", "fpr")

  auc_ROCR = performance(pred.gene, "auc")
  print("auc")
  print(auc_ROCR@y.values[[1]])

  roc.x = as.data.frame(unlist(rfperf.gene@x.values))
  roc.y = as.data.frame(unlist(rfperf.gene@y.values))
  roc.t = cbind(roc.x, roc.y)
  colnames(roc.t) = c("xval", "yval")

  output2 = list(cmfinal,roc.t)
  return(output2)
}

#save the output of the ensemble model
ens_B = ensemble_model(out_B[[1]], out_Bc[[1]], out_B[[2]], "BRAF")
ens_P = ensemble_model(out_P[[1]], out_Pc[[1]], out_P[[2]], "PIK3CA")
ens_E = ensemble_model(out_E[[1]], out_Ec[[1]], out_E[[2]], "EGFR")
ens_K = ensemble_model(out_K[[1]], out_Kc[[1]], out_K[[2]], "KRAS")

as.data.frame(out_B[[2]]$object)
as.data.frame(out_B[[2]]$Sample)
out_P[[2]]$object
as.data.frame(out_K[[2]]$Sample)
out_E[[2]]$object

out_K[[2]]$object


##PLOTS
##1. ROC curve
##2. Confusion matrix

roc.tB = ens_B[[2]]
roc.tP = ens_P[[2]]
roc.tE = ens_E[[2]]
roc.tK = ens_K[[2]]

## annotations
roc.tP$gene = "PIK3CA 0.81"
roc.tP$colour = "coral"
roc.tE$gene = "EGFR 0.62"
roc.tE$colour = "coral3"
roc.tB$gene = "BRAF 0.91"
roc.tB$colour = "coral4"
roc.tK$gene = "KRAS 0.82"
roc.tK$colour = "chocolate4"

total = rbind(roc.tB, roc.tP, roc.tE, roc.tK)

location = '~/Desktop/geneatlas/rf_model/results/'
png(paste(location, filename = "ROC.png", sep = ""), width = 10 , height = 5, units = 'in', res = 300)
ggplot(total, aes(x = xval, y= yval, group = gene, colour = gene)) +
  geom_line() +
  scale_color_brewer(palette = "PuOr") +
  #scale_colour_manual(values=c("coral", "coral3", "coral4", "chocolate4")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_text(size = 12, hjust = 0),
        axis.text.y=element_text(size = 12 , hjust = 0),
        legend.text = element_text(size = 18 , hjust = 0),
        text = element_text(size=12),
        legend.justification = c(1.1, -0.01), legend.position = c(1, 0),
        panel.border = element_rect(colour = "gray", fill=NA, size = 1),
        legend.title=element_blank(),
        legend.key = element_blank()) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  labs(title= "ROC curve", 
       x = "FPR", 
       y = "TPR")

dev.off()


cm.tB = ens_B[[1]]
cm.tP = ens_P[[1]]
cm.tE = ens_E[[1]]
cm.tK = ens_K[[1]]

#save confustion matrixes
location = '~/Desktop/geneatlas/rf_model/results/'
b = c(0,10, 50 ,2000, 2500)

cmfinal = cm.tE
gene = "KRAS"

png(paste(location, filename = paste("performance_", gene, ".png", sep = "")), width = 10 , height = 5, units = 'in', res = 300)
ggplot(data = cmfinal,
       mapping = aes(x = Prediction,
                      y = Reference)) +
    geom_tile(aes(fill = Freq)) + #, colour = "gray"
    #geom_abline(intercept = 0, slope = 2.5, color = "gray") +
    geom_text(data = subset(cmfinal, Freq != 0.0001), aes(label = sprintf("%1.0f", Freq)), vjust = 0.5, size = 5) +
    #labs(x = "Prediction", y = "Reference") + 
    theme(axis.text.x=element_text(angle = -45, hjust = 0)) +
    scale_fill_gradientn(colours = c("white", "#FF0000"),
                         trans = "log", # if your results aren't quite as clear as the above example
                         breaks=b, labels=format(b))
dev.off()





#function to test the individual models
conf_matrix = function(pred,testData ){
  #Confusion matrix
  cm = confusionMatrix(pred$class, testData$object)
  cmT = cm$table
  cmfinal = as.data.frame(cm$table)
  cmfinal$Freq[cmfinal$Freq == 0] = 0.0001
  print(cmT)
  
  TN = sum(cmT) - ((rowSums(cmT) + colSums(cmT)) - diag(cmT))
  FN = rowSums(cmT) - diag(cmT)
  FP = colSums(cmT) - diag(cmT)
  TP = diag(cmT)

  accuracy <- sum(diag(cmT)) / sum(cmT)
  Baccuracy = (TP / (TP+FN) + TN/ (TN+FP)) /2
  precision = diag(cmT) / colSums(cmT)
  recall = diag(cmT) / rowSums(cmT)

  print("accuracy")
  print(accuracy)
  print("")
  print("balanced accuracy")
  print(Baccuracy)
  print("")
  print("precision")
  print(precision)
  print("")
  print("recall")
  print(recall)
  print("")

  return(cmfinal)
}
conf_matrix(out_B[[1]], out_B[[2]])

prep_ROC = function(pred,testData, gene){
  prob.pred = as.data.frame(pred$predicted)
  prob.pred$labels = pred$class
  prob.pred$true = testData$object

  gene.prob = prob.pred[,c(gene, "labels", "true")]
  gene.prob$p_pred = as.factor(ifelse(gene.prob$labels == gene, 1, 0))
  gene.prob$actual = as.factor(ifelse(gene.prob$true == gene, 1, 0))

  pred.gene = prediction(as.vector(as.numeric(gene.prob[,1])), as.vector(as.numeric(gene.prob$actual)))
  rfperf.gene = performance(pred.gene, "tpr", "fpr")
  
  auc_ROCR = performance(pred.gene, "auc")
  print("auc")
  print(auc_ROCR@y.values[[1]])

  roc.x = as.data.frame(unlist(rfperf.gene@x.values))
  roc.y = as.data.frame(unlist(rfperf.gene@y.values))
  roc.t = cbind(roc.x, roc.y)
  colnames(roc.t) = c("xval", "yval")

  return(roc.t)
}
prep_ROC(out_B[[1]], out_B[[2]], "BRAF")


