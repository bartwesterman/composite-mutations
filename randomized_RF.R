###################################################################
### RANDOMIZED MODEL RECURRENT MUTATIONS                        ###
### random forest model for the top 4 genes that show the genes ###
### that show the most recurrent mutations. Both for chromosomal###
### loc and for the mutation (0 - 1 ) of co-mutations. Combine  ###
### models with bagging. Now we randomly shuffle the data, with ###
### increasing amount of 10# each time up till 100# suffled data###
### http://www.datasciencemadesimple.com/sample-function-in-r/  ###
###################################################################
library(reshape2)
#function to shuffle the columns
shuffle_c= function(df,i){
  data = df[ , -which(names(df) %in% c("Sample","object"))]
  part_data = sample(ncol(data) * i)
  new_data <- as.data.frame(apply(data[part_data], 2, sample)) #column
  rem_data = data[,-c(part_data)] #column
  total_data = cbind(new_data, rem_data) #column
  total_data$object = df$object
  return(total_data)
}


#input data after downsample
set.seed(100)
down_sample = function(input, gene){
  
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
  datasetT$object = as.factor(datasetT$object)
  return(datasetT)
}

set.seed(100)
test1 = down_sample(in_BRAFm, "BRAF")
test2 = down_sample(in_BRAFc, "BRAF")

set.seed(100)
test1 = down_sample(in_EGFRm, "EGFR")
test2 = down_sample(in_EGFRc, "EGFR")

test1 = down_sample(in_PIK3CAm, "PIK3CA")
test2 = down_sample(in_PIK3CAc, "PIK3CA")

test1 = down_sample(in_KRASm, "KRAS")
test2 = down_sample(in_KRASc, "KRAS")

#create list with shuffle-ing data in increasing amount from 10% till 100%
RF_random = function(test1){
  set.seed(100)
  d = list()
  counter = 1
  set.seed(100)
  for(i in seq(0.1, 1.0, 0.1)){
    total_data = shuffle_c(test1,i)
    d[[counter]] = total_data
    counter = counter +1
  }
  
  set.seed(100)
  train_test <- lapply(d, function(x) { 
    return(createDataPartition(as.character(x$object), p = 2/3, list = FALSE))
  })
  
  set.seed(100)
  model_list = purrr::map2(d, train_test, 
                           function(df, train_index)  {
                             rfsrc(object ~ .,  data = df[train_index,],
                                   var.used="all.trees", ntree=100,importance="none")
                           })
  
  testData = purrr::map2(d, train_test,
                         function(df, train_index) {
                           df[-train_index,]
                         })
  
  pred.conv = purrr::map2(model_list, testData, 
                          function(x, y)  {
                            predict(x , y)
                          })
  out = list(testData, pred.conv)
  return(out)
  
}

brafm = RF_random(test1)
brafc = RF_random(test2)



input_chr = brafc[[2]]
input_co = brafm[[2]]
testData = brafm[[1]]

prob.predB = list()
for (i in 1:10){
  
  chr_prob = as.data.frame(input_chr[[i]][["predicted"]])
  co_prob = as.data.frame(input_co[[i]][["predicted"]])
  truelabel = as.data.frame(testData[[i]][["object"]])
  colnames(truelabel) = "true"
  ##WEIGHTED VOTE
  #Taking weighted average of predictions
  prob.pred = ((co_prob*0.75)+(chr_prob*0.25))/2
  prob.pred$label = colnames(prob.pred)[max.col(prob.pred,ties.method="first")]
  prob.pred$true = truelabel[1]
  prob.predB[[i]] = prob.pred
}


cm = purrr::map2(prob.predB, testData, 
                 function(x, y)  {
                   #print(x)
                   #print(y[["conv"]])
                   confusionMatrix(factor(x$label) , y$object)
                 })

cm_random = function(cm, gene){
  cmT = lapply(cm, function(x) x["table"])
  
  
  TN = lapply(cmT, function(j) sum(j[[1]]) - ((rowSums(j[[1]]) + colSums(j[[1]])) - diag(j[[1]])))
  FN = lapply(cmT, function(k) colSums(k[[1]]) - diag(k[[1]]))
  FP = lapply(cmT, function(l) rowSums(l[[1]]) - diag(l[[1]]))
  TP = lapply(cmT, function(m) diag(m[[1]]))
  
  
  Baccuracy = as.data.frame(mapply(function(tn, fn, fp, tp) (tp/ (tp + fn) + tn/ (tn + fp)) /2 , TN, FN, FP, TP))
  precision = as.data.frame(lapply(cmT, function(i) diag(i[[1]] / rowSums(i[[1]]))))
  recall = as.data.frame(lapply(cmT, function(i) diag(i[[1]] / colSums(i[[1]]))))
  acc = as.data.frame(lapply(cmT, function(i) diag(i[[1]] / sum(i[[1]]))))
  
  colnames(Baccuracy) = c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  Baccuracy$gene = rownames(Baccuracy)
  return(Baccuracy)
}
braf_random = cm_random(cm)
egfr_random = cm_random(cm)
pik3ca_random = cm_random(cm)
kras_random = cm_random(cm)

total_cm = rbind(braf_random,egfr_random,pik3ca_random,kras_random)
s_baccuracy = total_cm[total_cm$gene %in% c("BRAF","EGFR", "PIK3CA", "KRAS"),] 
m_B_acc = melt(s_baccuracy)

location = '~/Desktop/geneatlas/rf_model/results/'
png(paste(location, filename = "balanced_accuracyrandom.png", sep = ""), width = 10 , height = 5, units = 'in', res = 300)
ggplot(m_B_acc, aes(x = variable, y = value, group = gene, color = gene))+
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("PIK3CA" = "#7030A0", "KRAS" = "#0070C0", "EGFR" = "#385723", "BRAF" =  "#FF7979")) +
  labs(title= "Random data - col", 
       x = "Randomization %", 
       y = "Balanced accuracy") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_text(size = 12, hjust = 0),
        axis.text.y=element_text(size = 12 , hjust = 0),
        #legend.text = element_text(size = 18 , hjust = 0),
        text = element_text(size=12),
        #legend.justification = c(1.1, -0.01), legend.position = c(1, 0),
        panel.border = element_rect(colour = "gray", fill=NA, size = 1),
        legend.title=element_blank(),
        legend.key = element_blank())
dev.off()

#ROC curve

gene = "KRAS"
roc_gene = list()
for(r in 1:10){
  gene.prob = prob.predB[[r]][,c(gene, "label", "true")]
  gene.prob$p_pred = as.factor(ifelse(gene.prob$label == gene, 1, 0))
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
  roc.t$random = 10 * r
  roc_gene[[r]] = roc.t
}

#create one dataframe of mutation data
names(roc_gene) <- NULL
roc_kras = do.call("rbind", roc_gene) 
roc_kras$random = as.character(roc_kras$random)

location = '~/Desktop/geneatlas/rf_model/results/'
png(paste(location, filename = "roc_randomKRAS", sep = ""), width = 10 , height = 5, units = 'in', res = 300)
ggplot(roc_kras, aes(x = xval, y= yval, group = random, colour = random)) +
  geom_line() +
  scale_color_manual(values = c("10" = "#4575b4", "20" = "#74add1", "30" = "#abd9e9", "40" =  "#e0f3f8", "50" = "black",
                                "60" = "#fee090", "70" = "#fdae61", "80" = "#f46d43", "90" =  "#d73027", "100" = "#b2182b")) +
  #scale_colour_manual(values=c("coral", "coral3", "coral4", "chocolate4")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_text(size = 12, hjust = 0),
        axis.text.y=element_text(size = 12 , hjust = 0),
        #legend.text = element_text(size = 18 , hjust = 0),
        text = element_text(size=12),
        #legend.justification = c(1.1, -0.01), legend.position = c(1, 0),
        panel.border = element_rect(colour = "gray", fill=NA, size = 1),
        legend.title=element_blank(),
        legend.key = element_blank()) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  labs(title= "ROC curve", 
       x = "FPR", 
       y = "TPR")

dev.off()
