###########################################
### preparation input data for RF model ###
###########################################

## load data 
## 1. co-mutation data from fishers
setwd("") #drive for the input data
filelist1 = list.files(pattern = ".*comp.txt")
datalist1 = lapply(filelist1, function(x)read.table(x, header=T)) 
dataset1 = do.call("rbind", datalist1) 
head(dataset1)
table(dataset1$tumor)


dataset.S = dataset1
dataset.S$tumor = 
  ifelse(dataset1$tumor == "blca_tcga", "bladder",
         ifelse(dataset1$tumor == "coadread_tcga", "bowel",
                ifelse(dataset1$tumor == "gbm_tcga_mutations", "brain",
                       ifelse(dataset1$tumor == "brca_tcga", "breast",
                              ifelse(dataset1$tumor == "cesc_tcga", "cervix",
                                     ifelse(dataset1$tumor == "stes_tcga_pub", "eso",
                                            ifelse(dataset1$tumor == "hnsc_tcga", "headneck",
                                                   ifelse(dataset1$tumor == "kirc_tcga", "kidney",
                                                          ifelse(dataset1$tumor == "lihc_tcga", "liver",
                                                                 ifelse(dataset1$tumor == "nsclc_tcga_broad_2016", "lung",
                                                                        ifelse(dataset1$tumor == "paad_tcga", "pancreas", 
                                                                               ifelse(dataset1$tumor == "prad_p1000", "prostate",
                                                                                      ifelse(dataset1$tumor == "skcm_tcga", "skin",
                                                                                             ifelse(dataset1$tumor == "ucec_tcga", "uterus", "0"))))))))))))))
table(dataset.S$tumor)
table(dataset1$tumor)

## 2. chr location data for all cases
path.input = "" #path to input_data
chr_braf = read.table(paste(path.input,  "BRAF_inputCHR.txt",sep=""), head = TRUE, sep = "\t", row.names = 1) #co-mutation genelist
chr_egfr = read.table(paste(path.input,  "EGFR_inputCHR.txt",sep=""), head = TRUE, sep = "\t", row.names = 1) #co-mutation genelist
chr_kras = read.table(paste(path.input,  "KRAS_inputCHR.txt",sep=""), head = TRUE, sep = "\t", row.names = 1) #co-mutation genelist
chr_pik3ca = read.table(paste(path.input,  "PIK3CA_inputCHR.txt",sep=""), head = TRUE, sep = "\t", row.names = 1) #co-mutation genelist
table(chr_pik3ca$tumor)
#chr_data = rbind(chr_braf, chr_egfr,chr_kras,chr_pik3ca)
# chr_data[is.na(chr_data)] = 0

## 3. mutation data for all cases
mut_data = read.table(paste(path.input,  "masterfile_mutcnv.txt",sep=""), head = TRUE, sep = "\t") #co-mutation genelist
mut_data[is.na(mut_data)] = 0

### convergent gene vs 0 : choose only tumor type with convergent gene!!!
# function to select interested convergent gene and split the data
sel_gene = function(input, gene){
  comut_set = dataset.S[which(dataset.S$conv == gene),]
  input_sample = input[,colnames(input) %in% c(as.character(comut_set$gene), "Sample", "object")]
  input_sample$object = droplevels(input_sample)$object
  input_sample[is.na(input_sample)] = 0
  return(input_sample)
}

##input chr data for the RF model
in_BRAFc = sel_gene(chr_braf, "BRAF")
in_EGFRc = sel_gene(chr_egfr, "EGFR")
in_KRASc = sel_gene(chr_kras, "KRAS")
in_PIK3CAc = sel_gene(chr_pik3ca, "PIK3CA")

##input mut data for the RF model

### convergent gene vs 0 : choose only tumor type with convergent gene!!!
# function to select interested convergent gene and split the data
sel_geneM = function(input, sample, gene){
  input_sample = input[which(input$Sample %in% sample$Sample),]
  input_sample1 = input_sample[,colnames(input_sample) != gene] #remove convergent gene feature because this will be the best feature to seperate the two groups
  input_sample1 = input_sample1[!duplicated(input_sample1$Sample),]
  
  comut_set = dataset.S[which(dataset.S$conv == gene),]
  input_sample1 = input_sample1[,colnames(input_sample1) %in% c(as.character(comut_set$gene), "Sample", "object")]
  input_sample1$object = droplevels(input_sample1)$object
  input_sample1[is.na(input_sample1)] = 0
  input_sample1 = input_sample1[which(input_sample1$object %in% c(gene, "0")),]
  return(input_sample1)
}


in_BRAFm = sel_geneM(mut_data, chr_braf, "BRAF")
in_EGFRm = sel_geneM(mut_data, chr_egfr, "EGFR")
in_KRASm = sel_geneM(mut_data, chr_kras, "KRAS")
in_PIK3CAm = sel_geneM(mut_data, chr_pik3ca, "PIK3CA")


in_BRAFc = in_BRAFc[which(in_BRAFc$Sample %in% in_BRAFm$Sample),]
in_EGFRc = in_EGFRc[which(in_EGFRc$Sample %in% in_EGFRm$Sample),]
in_KRASc = in_KRASc[which(in_KRASc$Sample %in% in_KRASm$Sample),]
in_PIK3CAc = in_PIK3CAc[which(in_PIK3CAc$Sample %in% in_PIK3CAm$Sample),]

