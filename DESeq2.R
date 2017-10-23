library(DESeq2)
library(plotly)
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

data = data[rowSums(data) > 225,]

Cinput = list()
for(name in c("normal","BPH","LGPIN","HGPIN","^T3$","T34","^T4$","T45","T5")){
  Cinput[[name]] = data[,grepl(name, design$Dx)]
}

names(Cinput) = c("normal","BPH","LGPIN","HGPIN","T3","T34","T4","T45","T5")

data.ordered = data.frame(matrix(NA, nrow = nrow(data), ncol = 1))
for(name in names(Cinput)){
  data.ordered = cbind(data.ordered, Cinput[[name]])
}

data.ordered = data.ordered[,-1]

colnames(data.ordered)

design$Id[is.na(pmatch(design$Id,colnames(data.ordered)))]

design.v = as.character()
for(name in names(Cinput)){
  design.v = c(design.v, rep(name, length(Cinput[[name]])))
}

design.df = data.frame(name = colnames(data.ordered), subtype = design.v)
design.df$CaseID = design$CaseId[match(colnames(data.ordered), design$Id)]
design.df$batch_nam = design.df$CaseID %in% c(9:14)
design.df$batch_nam = sub("TRUE","b2",design.df$batch)
design.df$batch_nam = sub("FALSE","b1",design.df$batch)
design.df$batch_nam = as.factor(design.df$batch_nam)

design.df$batch = design.df$CaseID %in% c(6,9:14)   #PCA analysis after DESeq vst shows there are 2 batches like that
design.df$batch = sub("TRUE","b2",design.df$batch)
design.df$batch = sub("FALSE","b1",design.df$batch)
design.df$batch[which(design.df$name == "P107.T34")] = "b2"
design.df$batch = as.factor(design.df$batch)

rownames(design.df) = design.df$name

# batch1 = read.table("/home/magda/Desktop/prostate_cancer_progression/data/170711_prostate_batch_1.csv", sep = ",", header = T) 
# batch2 = read.table("/home/magda/Desktop/prostate_cancer_progression/data/07272017_prostate_batch_2.csv", sep = ",", header = T)
# batch3 = read.table("/home/magda/Desktop/prostate_cancer_progression/data/2017_07_13_breast_cancer_bulk_and_prostate_batch3.csv", sep = ",", header = T)
# batch3 = batch3[17:32,]
# batch4 = read.table("/home/magda/Desktop/prostate_cancer_progression/data/08012017_Prostate_cancer_batch4.csv", sep = ",", header = T)
# batch2let = read.table("/home/magda/Desktop/prostate_cancer_progression/data/2016_11_10_prostate_cancer_batch2.csv", sep = ",", header = T)
# batch4let = read.table("/home/magda/Desktop/prostate_cancer_progression/data/2016_11_23_prostate_cancer_batch4.csv", sep = ",", header = T)
# 
# batches = list(batch1, batch2, batch3, batch4, batch2let, batch4let)
# names(batches) = c("batch1", "batch2", "batch3", "batch4", "batch2let", "batch4let")
# batches.df = data.frame(ID = 1, bh = "batch0")
# for(batch in names(batches)){
#   temp.df = data.frame(ID = batches[[batch]][,1], bh = batch)
#   batches.df = rbind(batches.df, temp.df)
# }
# batches.df = batches.df[2:173,]
# str(batches.df)
# batches.df$ID = sub("_",".", batches.df$ID)
# batches.df$ID = sub("_",".", batches.df$ID)
# batches.df$ID = gsub("\\_*","",batches.df$ID)
# merge(design.df, batches.df, by.x = )
# 
# 
# write.table(design.df, "design.df.txt", sep = "\t", row.names = F)
# write.table(data.ordered, "PCdataForDE.txt", sep = "\t", row.names = T)
# design.df = read.table(file.choose(), row.names = 1, sep = "\t")

TotRead=colSums(data.ordered)
barplot(TotRead, las=2, main="Total number of reads")
boxplot(log2(data.ordered+1), col=design.df$subtype, las=2)

pca1=prcomp(t(log2(data.ordered+1)))
plot(pca1$x[ ,1:2], col=design.df$batch, pch=as.numeric(design.df$subtype), main="orig data")

##################
log.data.ordered = log2(data.ordered+1)
centered = log.data.ordered - rowMeans(log.data.ordered)
class(centered)
centered.m = as.matrix(centered)
se <- SummarizedExperiment(centered.m,
                           colData = design.df)
PCA = plotPCA( DESeqTransform( se ), "batch" )
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))


##batch correction with edgeR
library(edgeR)
y=DGEList(counts=data.ordered)
y2=calcNormFactors(y)
logCPM=cpm(y2, log=T, prior.counts=1)
mod=model.matrix(~subtype, data=design.df)
logCPM=removeBatchEffect(logCPM, batch=design.df$batch, design=mod)

boxplot(logCPM, col=design.df$subtype, las=2)

pca2=prcomp(t(logCPM), scale=F)
plot(pca2$x[ ,1:2], col=design.df$batch, pch=as.numeric(design.df$subtype), main= "edgeR_corr")
plot(pca2$x[ ,1:2], col=design.df$subtype, pch=as.numeric(design.df$batch), main= "edgeR_corr")

SampC = design.df$subtype
str(SampC)
SampC=factor(SampC, levels=c("normal", "BPH", "LGPIN", "HGPIN", "T3", "T34","T4","T45","T5"))
library(rgl)
plot3d(pca2$x[ ,1], pca2$x[ ,2], pca2$x[ ,3],radius=15, type="s", col=as.numeric(SampC),
       xlab = "PC1", ylab = "PC2", zlab="PC3")

library(DESeq2)
library(plotly)
str(design.df$subtype)
design.df$CaseID = as.factor(design.df$CaseID)
str(design.df$CaseID)

dds <- DESeqDataSetFromMatrix(countData = data.ordered,
                              colData = design.df,
                              design= ~ subtype)
levels(dds$subtype)
dds$subtype = relevel(dds$subtype,"normal")

vsd <- vst(dds, blind=FALSE)

PCA = plotPCA(vsd, "batch_nam")
dataPCA = PCA$data
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))
p

PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))
p

#correct for batch
dds <- DESeqDataSetFromMatrix(countData = data.ordered,
                              colData = design.df,
                              design= ~ batch + subtype)

dds$sutype = relevel(dds$subtype,"normal")
vsd <- vst(dds, blind=FALSE)
vsdCounts = assay(vsd)
boxplot(vsdCounts, main="Vst transformed data")
PCA = plotPCA(vsd, "batch")


# perform correction using limma?
library(limma)
mod=model.matrix(~subtype, data = design.df)
newvstCounts=removeBatchEffect(vsdCounts, batch = design.df$batch, design=mod)
boxplot(newvstCounts, main="Vst transformed batch corrected")
se <- SummarizedExperiment(newvstCounts,
                           colData = design.df)
PCA = plotPCA( DESeqTransform( se ), "batch" )
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))

pca3=prcomp(t(newvstCounts), scale=F)
plot(pca3$x[ ,1:2], col=design.df$batch, pch=as.numeric(design.df$subtype), main="coloured by batch")
plot(pca3$x[ ,1:2], col=design.df$subtype, pch=as.numeric(design.df$batch), main="coloured by subtype")

#DE
library("BiocParallel")
register(MulticoreParam(4))
all(rownames(design.df) == colnames(data.ordered))
dds = DESeqDataSetFromMatrix(countData = data.ordered,
                             colData = design.df,
                             design = ~ batch + subtype)

dds$subtype = relevel(dds$subtype, "normal")

dds = DESeq(dds, parallel = TRUE)
res = results(dds, contrast=c("subtype", "T45", "normal"))
resultsNames(dds)
resLFC = lfcShrink(dds, contrast=c("subtype", "T45", "normal"), res = res)
DESeq2::plotMA(dds, ylim=c(-2,2))
DESeq2::plotMA(resLFC, ylim=c(-2,2))



#DEresN_BPH = results(dds, contrast=c("subtype", "BPH", "normal"))
DEresLG_HGPIN = results(dds, contrast=c("subtype", "HGPIN", "LGPIN"), parallel = TRUE)                          
DEresBPH_LGPIN = results(dds, contrast=c("subtype", "LGPIN", "BPH"))
DEresHGPIN_T3 = results(dds, contrast=c("subtype", "T3", "HGPIN"))
DEresT3_T34 = results(dds, contrast=c("subtype", "T34", "T3")) 
DEresT34_T4 = results(dds, contrast=c("subtype", "T4", "T34"))
DEresT4_T45 = results(dds, contrast=c("subtype", "T45", "T4"))
DEresT45_T5 = results(dds, contrast=c("subtype", "T5", "T45"))


sum(res.sort$padj < 0.05, na.rm=TRUE)
