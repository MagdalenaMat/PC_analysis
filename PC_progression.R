##PC_progression

design = read.table("/home/magda/Desktop/prostate_cancer_progression/design.file.prostate.225.b1-8.revised.csv", sep = ",", header = T) 
data = read.table("/home/magda/Desktop/prostate_cancer_progression/prostate.b1_b8.225.raw.nodup.csv", sep = ",", header = T, row.names = 1)


design$Id = sub("-",".", design$Id)
design$Id = sub("-",".", design$Id)
design = design[1:225,]
design$Id[1:116] = apply(design[1:116,], 1, function(x){paste0("X",x[1])})
design$Id == colnames(data)

design$Dx = sub("NN","normal", design$Dx)
design$Dx = sub("^N$","normal", design$Dx)
design$Dx[is.na(design$Dx)] = "normal"
design$Id == colnames(data)

# #normalize with rlog from DESeq2
# library(DESeq2)
# data = data[rowSums(data) > 0,]
# dds <- DESeq(DESeqDataSetFromMatrix(data, design, ~ Dx), parallel = T, fitType = "local")
# rld <- assay(rlog(dds, fitType = "local"))

#exclude rows with mostly 0s
#data = data[rowSums(data) > 200,]
#TPM
#tpm = 1E6 * sweep(data, 2, colSums(data), "/")


Cinput = list()
for(name in c("normal","BPH", "HGPIN","LGPIN","^T3$","T34","^T4$","T45","T5")){
  Cinput[[name]] = tpm[,grepl(name, design$Dx)]
}

names(Cinput) = c("normal", "BPH","HGPIN","LGPIN","T3","T34","T4","T45","T5")




Cinput[["HG_LG_PIN"]] = cbind(Cinput$HGPIN, Cinput$LGPIN)
Cinput[["T3_34"]] = cbind(Cinput$T3, Cinput$T34)
Cinput[["T45_T5"]] = cbind(Cinput$T45, Cinput$T5)
Cinput[["T4_45_5"]] = cbind(Cinput$T4, Cinput$T45, Cinput$T5)
lengths(Cinput)

for(name in names(Cinput)){
  Cinput[[name]] = cbind(rownames(Cinput[[name]]), Cinput[[name]])
  Cinput[[name]][,1] = gsub("\\..*","",Cinput[[name]][,1])
  Cinput[[name]] = Cinput[[name]][!duplicated(Cinput[[name]][,1]),]
}

setwd("/home/magda/Desktop/prostate_cancer_progression/data")
#write to the dics
for(name in names(Cinput)){
  write.table(Cinput[[name]], paste0("PC_",name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

write.table(Cinput$T45_T5, "PC_T45_T5.txt", sep = "\t", quote = FALSE, row.names = FALSE)
