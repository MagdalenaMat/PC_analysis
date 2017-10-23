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


########### gene centered and log2
log.data.ordered = log2(data.ordered+1)
centered = log.data.ordered - rowMeans(log.data.ordered)
class(centered)
centered.m = as.matrix(centered)
se <- SummarizedExperiment(centered.m,
                           colData = design.df)
PCA = plotPCA( DESeqTransform( se ), "batch" )
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))%>%
  layout(title = "log2 transform and centered original data")

##batch correction with edgeR
library(edgeR)
y=DGEList(counts=data.ordered)
y2=calcNormFactors(y)
logCPM=cpm(y2, log=T, prior.counts=1)
mod=model.matrix(~subtype, data=design.df)
logCPM=removeBatchEffect(logCPM, batch=design.df$batch, design=mod)

se <- SummarizedExperiment(logCPM,
                           colData = design.df)
PCA = plotPCA( DESeqTransform( se ), "batch" )
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))%>%
  layout(title = "edgeR batch corrected")

###########vst
design.df$subtype = factor(design.df$subtype, levels = c("normal", "BPH", "LGPIN", "HGPIN", "T3","T34", "T4", "T45", "T5"))
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
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))%>%
  layout(title = "DESeq2::vst transformed original data")
p

PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))%>%
  layout(title = "DESeq2::vst transformed original data")
p

#correct for batch
class(design.df$subtype)
dds <- DESeqDataSetFromMatrix(countData = data.ordered,
                              colData = design.df,
                              design= ~ batch + subtype)

dds$sutype = relevel(dds$subtype,"normal")
vsd <- vst(dds, blind=FALSE)
vsdCounts = assay(vsd)
#boxplot(vsdCounts, main="Vst transformed data")
PCA = plotPCA(vsd, "batch")


# perform correction using limma?
library(limma)
mod=model.matrix(~subtype, data = design.df)
newvstCounts=removeBatchEffect(vsdCounts, batch = design.df$batch, design=mod)
#boxplot(newvstCounts, main="Vst transformed batch corrected")
se <- SummarizedExperiment(newvstCounts,
                           colData = design.df)
PCA = plotPCA( DESeqTransform( se ), "batch" )
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"))%>%
  layout(title = "DESeq2::vst transformed limma batch corrected")


PCA = plotPCA( DESeqTransform( se ), "CaseID" )
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = c("red","blue"), marker = list(size = 20))%>%
  layout(title = "DESeq2::vst transformed limma batch corrected, color ~ stage")

plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set1", marker = list(size = 20))%>%
  layout(title = "DESeq2::vst transformed limma batch corrected, color ~ CaseID")
