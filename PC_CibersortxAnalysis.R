mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


#############
path = "/home/magda/Desktop/cibersortX/PC/CV_filtered"
setwd(path)
files <- list.files(path=path)

#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_PC_",".txt_GEPs_CVFilteredAll.txt"),c("",""),i)
  namesGEP = c(namesGEP, nam)
}

#put all files in the list
filteredGEP = list()
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  filteredGEP[[i]] = r
}
names(filteredGEP) = namesGEP

gene_count = data.frame()
for(i in namesGEP){
  gene_count = rbind(gene_count, apply(filteredGEP[[i]], 2, function(x){sum(!is.na(x))}))
}

rownames(gene_count) = namesGEP
colnames(gene_count) =colnames(filteredGEP[[1]])

#extract all gene names from the files and put them in the list to find intersect 
#of genes acros all the histological subtypes
results = list()
for (i in namesGEP){
  results[[i]] = as.character(filteredGEP[[i]][,1])
}

universe = Reduce(intersect, results)

#extract genes in monocyte signature common for all samples 
Mono_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Mono = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,6)] #select gene_names and Monocytes
  Mono_ex = cbind(Mono_ex,Mono)
}

all(as.character(Mono_ex[,2])==as.character(Mono_ex[,4]))
rownames(Mono_ex) = as.character(Mono_ex[,2])
Mono_ex = Mono_ex[, seq(3, length(Mono_ex), 2)]
colnames(Mono_ex) = namesGEP

Mono_ex[Mono_ex<=1] <- NA
predicted_genes_Monocytes = apply(Mono_ex, 2, function(x){sum(!is.na(x))})

#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 2,]
al3_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 3,]

#extract genes in NK cells signature common for all samples##################################################################################################
NK_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  NK = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,5)] #select gene_names and NK
  NK_ex = cbind(NK_ex,NK)
}

all(as.character(NK_ex[,2])==as.character(NK_ex[,4]))
rownames(NK_ex) = as.character(NK_ex[,2])
NK_ex = NK_ex[, seq(3, length(NK_ex), 2)]
colnames(NK_ex) = namesGEP

NK_ex[NK_ex<=1] <- NA
predicted_genes_NK = apply(NK_ex, 2, function(x){sum(!is.na(x))})
predicted_genes_NK
#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_NK = NK_ex[rowSums(!is.na(NK_ex)) >= 2,]

###############################################################################################################################################################

NT_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  NT = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,10)] #select gene_names and NT
  NT_ex = cbind(NT_ex,NT)
}

all(as.character(NT_ex[,2])==as.character(NT_ex[,4]))
rownames(NT_ex) = as.character(NT_ex[,2])
NT_ex = NT_ex[, seq(3, length(NT_ex), 2)]
colnames(NT_ex) = namesGEP

NT_ex[NT_ex<=1] <- NA
predicted_genes_NT = apply(NT_ex, 2, function(x){sum(!is.na(x))})
predicted_genes_NT
#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_NT = NT_ex[rowSums(!is.na(NT_ex)) >= 2,]

################################################################################################################################################################
Tc_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Tc = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,4)] #select gene_names and Tcells
  Tc_ex = cbind(Tc_ex,Tc)
}

all(as.character(Tc_ex[,2])==as.character(Tc_ex[,4]))
rownames(Tc_ex) = as.character(Tc_ex[,2])
Tc_ex = Tc_ex[, seq(3, length(Tc_ex), 2)]
colnames(Tc_ex) = namesGEP

Tc_ex[Tc_ex<=1] <- NA
predicted_genes_Tc = apply(Tc_ex, 2, function(x){sum(!is.na(x))})
predicted_genes_Tc
#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_Tc = Tc_ex[rowSums(!is.na(Tc_ex)) >= 2,]


##################################################################################################################################################################

B_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  B = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,2)] #select gene_names and NT
  B_ex = cbind(B_ex,B)
}

all(as.character(B_ex[,2])==as.character(B_ex[,4]))
rownames(B_ex) = as.character(B_ex[,2])
B_ex = B_ex[, seq(3, length(B_ex), 2)]
colnames(B_ex) = namesGEP

B_ex[B_ex<=1] <- NA
predicted_genes_B = apply(B_ex, 2, function(x){sum(!is.na(x))})
predicted_genes_B
#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_B = B_ex[rowSums(!is.na(B_ex)) >= 2,]

# library(heatmaply)
# heatmaply_na(al3_mono, showticklabels = c(T,F))

##########
#plot immune related gene values in different stages 
###########
library(plotly)
library(reshape2)
library(dplyr)
library(ggplot2)
to_reshape = data.frame(gene.names = rownames(al2_mono), al2_mono)
to_reshape = data.frame(gene.names = rownames(al2_NK), al2_NK)
to_reshape = data.frame(gene.names = rownames(al2_NT), al2_NT)
to_reshape = data.frame(gene.names = rownames(al2_Tc), al2_Tc)
to_reshape = data.frame(gene.names = rownames(al2_B), al2_B)
colnames(to_reshape)


#compute log2 Fold Changes
change = mutate(to_reshape, normal_BPH = log2(to_reshape$BPH/to_reshape$normal), 
                BPH_LGPIN = log2(to_reshape$LGPIN/to_reshape$BPH),
                LPGIN_HGPIN = log2(to_reshape$HGPIN/to_reshape$LGPIN),
                BPH_LGHGPIN = log2(to_reshape$HG_LG_PIN/to_reshape$BPH),
                HPGIN_T3 = log2(to_reshape$T3/to_reshape$HGPIN),
                LGHGPIN_T3T34 = log2(to_reshape$HG_LG_PIN/to_reshape$T3_34),
                T3_T34 = log2(to_reshape$T34/to_reshape$T3),
                T34_T4 = log2(to_reshape$T4/to_reshape$T34),
                T4_T45 = log2(to_reshape$T45/to_reshape$T4),
                T45_T5 = log2(to_reshape$T5/to_reshape$T45),
                T3T34_T4T45T5 = log2(to_reshape$T4_45_5/to_reshape$T3_34),
                T3T34_T45T5 = log2(to_reshape$T45_T5/to_reshape$T3_34))

to_GO = list()
for(nm in tail(names(change),11)){
  to_GO[[paste0("DE_",nm,"_up")]] = change[which(change[[nm]] >= 0.5), 1]
  to_GO[[paste0("DE_",nm,"_down")]] = change[which(change[[nm]] <= -0.5), 1]
}

lengths(to_GO)

library("Rgraphviz")
library("DOSE")
library("clusterProfiler")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)


to_GO_EI = list()
for(element in names(to_GO)){
  to_GO_EI[[element]] = bitr(to_GO[[element]], 
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
}

GO_res = list()
for(element in names(to_GO_EI)){
  GO_res[[element]] = enrichGO(gene         = to_GO_EI[[element]]$ENTREZID,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable = TRUE)
}

for(element in names(GO_res)){
  plot(dotplot(GO_res[[element]],showCategory=50))
}

###Tcell results
GO_res$DE_LGHGPIN_T3T34_up@result[,c("ID","Description","geneID")]# Tcells ERBB signaling GO:0038127 GO:1901184 GO:1901185 
GO_res$DE_T3T34_T4T45T5_up@result[,c("ID","Description","geneID")]#neutrophils GO:0043312 GO:0042119 GO:0036230 GO:0002446 
                                                                  #ER UPR GO:0030968  GO:0034620 GO:0034976 GO:0036498
GO_res$DE_T3T34_T45T5_up@result[,c("ID","Description","geneID")]#TGFb responce GO:0071560 GO:0007179 #estradiol responce GO:0032355 #platelet activation GO:0030168
                                                                # ER UPR GO:0030968 GO:0006986 GO:0034620
