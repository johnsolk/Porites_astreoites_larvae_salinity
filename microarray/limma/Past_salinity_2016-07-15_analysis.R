library(limma)
#library("arrayQualityMetrics") # need to install this 07/03/2016
#library(preprocessCore) # need to install this 07/03/2016
library("Biobase")
library("RColorBrewer")
library("gplots")
library(ggplot2)
library(lattice)
# work laptop path:
setwd("~/Documents/HBOI/Past_salinity/salinity_microarray/")
# home Ubuntu laptop path:
# setwd("~/Documents/salinity_microarray")

# import raw data, which were extracted by hand
# (in the future, don't do this. Use automated function to extract from raw GenePix files, which we no longer have for this data set...)
salinity_data<-read.csv("Past_salinity_microarray_data.csv")
# convert data structure into data frame type
salinity_dataframe<-data.frame(salinity_data)
# sort by RefNum
salinity_dataframe_sort<-salinity_dataframe[order(salinity_dataframe$RefNum),]
head(salinity_dataframe_sort)
# set rownames to RefNum
rownames(salinity_dataframe_sort)<-salinity_dataframe_sort$RefNumber
head(salinity_dataframe_sort)
colnames(salinity_dataframe_sort)
head(salinity_dataframe_sort)
# dimensions of the data frame should be
# [1] 15744    16
dim(salinity_dataframe_sort)
# get rid of data that were flagged "BAD" = -100
salinity.filtered.flags <- salinity_dataframe_sort[salinity_dataframe_sort$Flags!=-100,]
head(salinity.data.filtered.flags)
dim(salinity.data.filtered.flags)
salinity<-salinity.data.filtered.flags[,c(8:15)]
dim(salinity)
head(salinity)
colnames(salinity)

# normalize
normDat<-normalizeBetweenArrays(salinity.data.matrix,method="quantile")
# or Loess
#normDat<-normalizeBetweenArrays(salinity.data.matrix,method="cyclicloess")
head(normDat)
###


# import microarray annotation file
annotation<-read.csv("FeatAnnotFile_SEE_10-18-12.csv")
dim(annotation)
# assume salinity.data.filtered.flags is in same order as annotation file
# this filters only annotations in data file
# cuts out some rows
# [1] 15732    23
annotation.filtered<-annotation[,][match(rownames(normDat),annotation$RefNumber),]
dim(annotation.filtered)
head(annotation.filtered)
class(annotation.filtered)

dim(normDat)
colnames(normDat)
#plotMD(tmp_data)

# set all negative values to 0
salinity.data[salinity.data<0]=0
head(salinity.data)

# transform all data + 1
salinity.data<-normDat+1
head(salinity.data)
# log base2 transform all data
salinity.data.log<-log2(salinity.data)
head(salinity.data.log)
class(salinity.data.log)

###
# merge annotation with normalized data
# so that data frame can be collapsed by probe
RefNumber<-rownames(salinity.data.log)
salinity.norm.data.log.all_merge<-as.data.frame(cbind(salinity.data.log,RefNumber))
salinity.norm.data.all.ann<-merge(salinity.norm.data.log.all_merge,annotation.filtered,by="RefNumber")

# remove controls
annot.keep<-salinity.norm.data.all.ann[salinity.norm.data.all.ann$ControlType=="FALSE",]
colnames(annot.keep)
rownames(annot.keep)<-annot.keep$RefNumber
# displays number of annotations that are not controls, should be 12,345
dim(annot.keep)
head(annot.keep)
annot.keep.data<-annot.keep[,c(2:15)]
# displays number of unique ID (probes), should be 4118
length(unique(annot.keep.data$ID))

colnames(annot.keep)
head(annot.keep)
class(annot.keep)
salinity.data.probe.rep<-aggregate(annot.keep.data[,c(1:9)], by=annot.keep.data['ID'], FUN = function(x) mean(as.numeric(as.character(x))))
#salinity.data.probe.rep<-aggregate(annot.keep.data[,c(1:9)], by=annot.keep.data['Name'], FUN = function(x) mean(as.numeric(as.character(x))))
dim(salinity.data.probe.rep)
head(salinity.data.probe.rep)

# PCA
#x<-normDat
x<-salinity.data.probe.rep
colnames(x)
colnames(x)<-c("sal25ppt_A","sal25ppt_B","sal25ppt_C","sal30ppt_A","sal30ppt_B","sal30ppt_C","sal35ppt_A","sal35ppt_B")
colnames(x)
pca = prcomp(t(x))
names = colnames(x)
#fac = factor(sapply(names,function(x){strsplit(x,'Rep')[[1]][1]}))
fac= factor(c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt"))
colours = c("red","blue","green")
xyplot(
  PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
  panel=function(x, y, ...) {
    panel.xyplot(x, y, ...);
    ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=1)
  },
  aspect = "fill", col=colours
  #main = draw.key(key = list(rect = list(col = list(col=colours), text = list(levels(fac)), rep = FALSE)))
)

#

# put ID as rownames
rownames(salinity.data.probe.rep)<-salinity.data.probe.rep$ID
colnames(salinity.data.probe.rep)
head(salinity.data.probe.rep)
# only data:
salinity.data<-salinity.data.probe.rep[,c(2:9)]
colnames(salinity.data)
head(salinity.data)
dim(salinity.data)

# remove low expressing probes
dim(salinity.data[rowMeans(salinity.data)>1,])
salinity.norm.filtered<-salinity.data[rowMeans(salinity.data)>1,]
head(salinity.norm.filtered)
dim(salinity.norm.filtered)
colnames(salinity.norm.filtered)

# PCA
x<-salinity.norm.filtered
#x<-norm.data.filtered
colnames(x)
colnames(x)<-c("sal25ppt_A","sal25ppt_B","sal25ppt_C","sal30ppt_A","sal30ppt_B","sal30ppt_C","sal35ppt_A","sal35ppt_B")
colnames(x)
pca = prcomp(t(x))
names = colnames(x)
#fac = factor(sapply(names,function(x){strsplit(x,'Rep')[[1]][1]}))
fac= factor(c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt"))
colours = c("red","blue","green")
xyplot(
  PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
  panel=function(x, y, ...) {
    panel.xyplot(x, y, ...);
    ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=1)
  },
  aspect = "fill", col=colours
  #main = draw.key(key = list(rect = list(col = list(col=colours), text = list(levels(fac)), rep = FALSE)))
)

head(annotation.filtered)
head(salinity.norm.filtered)
annot.keep<-subset(annotation.filtered,annotation.filtered$ID %in% rownames(salinity.data))
dim(annot.keep)
colnames(annot.keep)
annot.keep<-annot.keep[!duplicated(annot.keep$ID),][,-c(22,23)]
annot.keep<-annot.keep[order(annot.keep$ID),]
rownames(annot.keep)<-annot.keep$ID
salinity.norm.filtered<-salinity.norm.filtered[!duplicated(rownames(salinity.norm.filtered)),]
dim(annot.keep)
dim(salinity.norm.filtered)
an<-new("AnnotatedDataFrame",data=annot.keep)
experimentData<-new("MIAME",name="Sara Edge, Lisa Cohen, Ana María González Angel",
                    lab="Marine Genomics, HBOI",
                    contact="ljcohen@ucdavis.edu",
                    title="Porites astreoides larvae actuve salinity exposure: 25 ppt, 30 ppt, 35 ppt",
                    abstract="ExpressionSet for P. astreoides Microarray Data",
                    other=list(notes="Extracted from raw data with GenePixPro software"))
# no outliers removed:
#sample_type<-c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt","sal35ppt")
# remove outlier
sample_type<-c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt")
# remove: ,"c253295110054_ar2_jul1112_C_C"
sample_name<-c("c253295110048_ar2_jun1912_A_SL","c253295110048_ar1_jun1912_B_SL","c253295110048_ar4_jun1912_C_SL",
               "c253295110048_ar3_jun1912_A_L","c253295110048_ar6_jun1912_B_L","c253295110048_ar5_jun1912_C_L",
               "c253295110048_ar8_jun1912_A_C","c253295110048_ar7_jun1912_B_C")
sample_info<-cbind(sample_type,sample_name)
sample_info<-as.data.frame(sample_info)
rownames(sample_info)<-sample_name
sample_info
pd<-new("AnnotatedDataFrame",data=sample_info)

# unfiltered, unnormalized
#salinity_ExpressionSet_unfiltered<-new("ExpressionSet",exprs=salinity.norm,phenoData=pd,experimentData=experimentData,featureData=an)
# filtered, normalized
salinity.norm.filtered_matrix<-as.matrix(salinity.data)
dim(salinity.norm.filtered_matrix)
dim(annot.keep)
rownames(annot.keep)
colnames(salinity.norm.filtered_matrix)

salinity_ExpressionSet<-new("ExpressionSet",exprs=salinity.norm.filtered_matrix,phenoData=pd,experimentData=experimentData,featureData=an)
#arrayQualityMetrics(expressionset = salinity_ExpressionSet_filtered,outdir = "salinity_ArrayQualityMetrics_log_collapsed_by_transcriptName_loessnorm_filtered_one35pptoutlierremoved",force = TRUE,intgroup = c("sample_type"),do.logtransform=FALSE)
###
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3)))
colnames(design) <- c("sal35ppt", "sal30ppt", "sal25ppt")
design
fit <- lmFit(salinity_ExpressionSet_filtered,design)
contrast.matrix <- makeContrasts(sal25ppt-sal30ppt, sal25ppt-sal35ppt, sal30ppt-sal35ppt, levels= design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2
tt_25ppt_v_30ppt<-topTable(fit2,coef=1,adjust = "BH",n=5000)
tt_25ppt_v_35ppt<-topTable(fit2,coef=2,adjust = "BH",n=5000)
tt_30ppt_v_35ppt<-topTable(fit2,coef=3,adjust = "BH",n=5000)
tt_25ppt_v_30ppt_sig<-tt_25ppt_v_30ppt[order(-tt_25ppt_v_30ppt$adj.P.Val),]
tt_25ppt_v_35ppt_sig<-tt_25ppt_v_35ppt[order(-tt_25ppt_v_35ppt$adj.P.Val),]
tt_30ppt_v_35ppt_sig<-tt_30ppt_v_35ppt[order(-tt_30ppt_v_35ppt$adj.P.Val),]
write.csv(tt_25ppt_v_30ppt_sig,"P_ast_microarray_25pptvs30ppt.csv")
write.csv(tt_25ppt_v_35ppt_sig,"P_ast_microarray_25pptvs35ppt.csv")
write.csv(tt_30ppt_v_35ppt_sig,"P_ast_microarray_30pptvs35ppt.csv")



colnames(tt_25ppt_v_30ppt)
norm.data.reps<-cbind(norm.data.filtered,rownames(norm.data.filtered))
dim(norm.data.reps)
head(norm.data.reps)
colnames(norm.data.reps)<-c("c253295110048_ar2_jun1912_A_SL","c253295110048_ar1_jun1912_B_SL","c253295110048_ar4_jun1912_C_SL",
                            "c253295110048_ar3_jun1912_A_L","c253295110048_ar6_jun1912_B_L","c253295110048_ar5_jun1912_C_L",
                            "c253295110048_ar8_jun1912_A_C","c253295110048_ar7_jun1912_B_C","ID")
colnames(norm.data.reps)
head(norm.data.reps)


tt_25ppt_v_30ppt_reps<-merge(norm.data.reps,tt_25ppt_v_30ppt,by="ID")
tt_25ppt_v_35ppt_reps<-merge(norm.data.reps,tt_25ppt_v_35ppt,by="ID")
tt_30ppt_v_35ppt_reps<-merge(norm.data.reps,tt_30ppt_v_35ppt,by="ID")
dim(tt_25ppt_v_30ppt_reps)
dim(tt_25ppt_v_35ppt_reps)
dim(tt_30ppt_v_35ppt_reps)
#write.csv(tt_25ppt_v_30ppt_reps,"P_ast_microarray_25pptvs30ppt_reps_byprobeID_quantile.csv")
#write.csv(tt_25ppt_v_35ppt_reps,"P_ast_microarray_25pptvs35ppt_reps_byprobeID_quantile.csv")
#write.csv(tt_30ppt_v_35ppt_reps,"P_ast_microarray_30pptvs35ppt_reps_byprobeID_quantile.csv")
# filter
tt_25ppt_v_30ppt_sig<-subset(tt_25ppt_v_30ppt_reps,tt_25ppt_v_30ppt_reps$adj.P.Val<0.1)
tt_25ppt_v_35ppt_sig<-subset(tt_25ppt_v_35ppt_reps,tt_25ppt_v_35ppt_reps$adj.P.Val<0.1)
tt_30ppt_v_35ppt_sig<-subset(tt_30ppt_v_35ppt_reps,tt_30ppt_v_35ppt_reps$padj.P.Val<0.1)
tt_25ppt_v_30ppt_sig<-tt_25ppt_v_30ppt_sig[order(-tt_25ppt_v_30ppt_sig$adj.P.Val),]
tt_25ppt_v_35ppt_sig<-tt_25ppt_v_35ppt_sig[order(-tt_25ppt_v_35ppt_sig$adj.P.Val),]
tt_30ppt_v_35ppt_sig<-tt_30ppt_v_35ppt_sig[order(-tt_30ppt_v_35ppt_sig$adj.P.Val),]
dim(tt_25ppt_v_30ppt_sig)
dim(tt_25ppt_v_35ppt_sig)
dim(tt_30ppt_v_35ppt_sig)
colnames(tt_25ppt_v_35ppt_sig)
dim(subset(tt_25ppt_v_35ppt_sig,logFC< -1 | logFC>1))
logFC_1<-subset(tt_25ppt_v_35ppt_sig,logFC< -1 | logFC>1)
head(logFC_1)
dim(logFC_1)
#write.csv(logFC_1,"P_ast_salinity_67sig_logFC1.csv")
#write.csv(logFC_1,"P_ast_salinity_54sig_logFC1.csv")
source('~/Documents/scripts/overLapper.R')
setlist <- list(sal25ppt_v_sal30ppt_padj0.25=as.vector(tt_25pptv_30ppt_sig$Name),sal25ppt_v_sal35ppt_padj0.25=as.vector(tt_25pptv_35ppt_sig$Name),sal25ppt_v_sal35ppt_padj0.25=as.vector(tt_25pptv_30ppt_sig$Name))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
sal25ppt_35ppt_overlap<-intersect(tt_25pptv_30ppt_sig$Name,tt_25pptv_35ppt_sig$Name)
length(sal25ppt_35ppt_overlap)
#####
# make heatmap of 25ppt vs. 35 ppt (145 genes)
head(logFC_1)
dim(logFC_1)
d <- as.matrix(logFC_1[,c(2:9)])
class(d)<-"numeric"
rownames(d) <- logFC_1[,18]
head(d)
colnames(d)<-c("25ppt_A","25ppt_B","25ppt_C","30ppt_A","30ppt_B","30ppt_C","35ppt_A","35ppt_B")
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
#png("Salinity_heatmap_1.png", width = 7*300,height = 7*300,res = 1200,pointsize = 2)
heatmap.2(d, main="25ppt vs. 35 ppt (padj<0.25,log2FC+-1)", 
          Rowv=as.dendrogram(hr),
          cexRow=0.5,cexCol=1,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
#dev.off()
sessionInfo()
