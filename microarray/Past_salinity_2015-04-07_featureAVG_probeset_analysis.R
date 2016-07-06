library(limma)
library("arrayQualityMetrics")
library(preprocessCore)
library("Biobase")
library("RColorBrewer")
library("gplots")
# work laptop path:
setwd("~/Documents/HBOI/Past_salinity/salinity_microarray/")
# home Ubuntu laptop path:
# setwd("~/Documents/salinity_microarray")
# data extracted by hand
salinity_data<-read.csv("Past_salinity_microarray_data.csv")
salinity_dataframe<-data.frame(salinity_data)
salinity_dataframe_sort<-salinity_dataframe[order(salinity_dataframe$RefNum),]
head(salinity_dataframe_sort)
rownames(salinity_dataframe_sort)<-salinity_dataframe_sort$RefNumber
head(salinity_dataframe_sort)
colnames(salinity_dataframe_sort)
head(salinity_dataframe_sort)
# get rid of data that were flagged "BAD" = -100
salinity.data.filtered.flags <- salinity_dataframe_sort[salinity_dataframe_sort$Flags!=-100,]
head(salinity.data.filtered.flags)


# annotation
annotation<-read.csv("FeatAnnotFile_SEE_10-18-12.csv")
dim(annotation)
# assume salinity.data.filtered.flags is in same order as annotation file
# this filters only annotations in data file
# cuts out some rows
annotation.filtered<-annotation[,][match(rownames(salinity.data.filtered.flags),annotation$RefNumber),]
dim(annotation.filtered)
# remove controls
annot.keep<-annotation[annotation$ControlType=="FALSE",]
colnames(annot.keep)
# displays number of annotations that are not controls, should be 12,354
dim(annot.keep)
# displays number of unique gene Names, should be 2206
length(unique(annot.keep$Name))

# Remove data rows corresponding to control probes
# Do we want to look at control probe intensities? Yes. How?
salinity.data.filtered<-salinity.data.filtered.flags[,][match(annot.keep$RefNumber,rownames(salinity.data.filtered.flags)),]

# Collapse to both probe replicate and probeset levels:
# 1. collapse all probe reps
# probe rep = ID col
# group by ID (should be 3 reps of each ID, maybe less filtered)
# This is the number of unique ID:
# should be 4119
u<-unique(salinity.data.filtered$ID)
length(u)
# for each item in u
# salinity.data.rep<-subset(salinity.data.filtered, mean at each c(8:16) 
# mean of all rep intensities
salinity.data.probe.rep<-aggregate(salinity.data.filtered[,c(8:16)], by=salinity.data.filtered['ID'], mean)
class(salinity.data.probe.rep)
dim(salinity.data.probe.rep)

# 2. Take ID column in annotation, split by "_"
#salinity.data.probeset.ID <- strsplit(as.character(salinity.data.filtered$ID),'_',fixed=TRUE)
# Example ID syntax:
# CUST_3_PI426266566
# Get unique 3rd position in ID
#v<-unique(sapply(strsplit(as.character(salinity.data.filtered$ID),'_',fixed=TRUE),"[",3))
#length(v)
# Example Name syntax:
# 75858820
# Looks like there are 9 replicates of each name, probes designed at 3 positions on each transcript, replicated 3 times
# This is the number of unique Names:
#x<-unique(salinity.data.filtered$Name)
#length(x)
# http://stackoverflow.com/questions/1355355/how-to-avoid-a-loop-in-r-selecting-items-from-a-list
#t <- c("bob_smith","mary_jane","jose_chung","michael_marx","charlie_ivan")
#pieces <- strsplit(t,"_")
#sapply(pieces, "[", 3)
salinity.probeset.group<-aggregate(salinity.data.filtered[,c(8:16)], by=list(sapply(strsplit(as.character(salinity.data.filtered$ID),'_',fixed=TRUE),"[",3)), mean)
dim(salinity.probeset.group)
salinity.probeset.name<-aggregate(salinity.data.filtered[,c(8:16)], by=salinity.data.filtered['Name'], mean)
dim(salinity.probeset.name)
# probeset is 3rd position
# group by third position
# mean of intensities for all matching third position
# this is "probeset" = same trancript/gene?
# Analyze data collapsed by name = salinity.probeset.name
# To-Do for another time:
# it would be nice to see variation between individual probes for each Name 

head(salinity.probeset.name)

# put Name as rownames
rownames(salinity.probeset.name)<-gsub("[[:space:]]", "",salinity.probeset.name$Name)
colnames(salinity.probeset.name)
# only data:
salinity.data<-salinity.probeset.name[,c(2:10)]
colnames(salinity.data)
# remove outliers:
#Column 3 is main outlier to remove
salinity.data<-salinity.data[,-9]
#salinity.data<-salinity.data[,-3]
head(salinity.data)
# set all negative values to 0

salinity.data[salinity.data<0]=0
head(salinity.data)

# transform all data + 1
salinity.data<-salinity.data+1
head(salinity.data)
# log base2 transform all data
salinity.data.log<-log2(salinity.data)
head(salinity.data.log)


salinity.data.matrix<-as.matrix(salinity.data.log)
salinity.data.matrix<-na.omit(salinity.data.matrix)
head(salinity.data.matrix)
salinity.data.matrix<-salinity.data.matrix[order(rownames(salinity.data.matrix)),]


head(annot.keep)
annot.keep<-subset(annot.keep,annot.keep$Name %in% rownames(salinity.data.matrix))
annot.keep<-annot.keep[!duplicated(annot.keep$Name),]
dim(annot.keep)
rownames(annot.keep)<-gsub("[[:space:]]", "", annot.keep$Name)
annot.keep<-annot.keep[order(rownames(annot.keep)),]
head(annot.keep)
# normalize
# quantile looks slightly better in arrayQualityMetrics
normDat<-normalizeBetweenArrays(salinity.data.matrix,method="quantile")
# Loess not so bad either
#normDat<-normalizeBetweenArrays(salinity.data.matrix,method="cyclicloess")
head(normDat)
salinity.norm<-normDat
#salinity.norm<-normalize.quantiles(salinity.data.matrix)
#salinity.norm<-normalizeCyclicLoess(salinity.data.matrix)
head(salinity.norm)
#salinity.norm <- backgroundCorrect(salinity.data.matrix, method="normexp", normexp.method="saddle",offset=5)
dim(salinity.norm)

#coln<-colnames(salinity.data.matrix)
#rown<-rownames(salinity.data.matrix)
#colnames(salinity.norm)<-coln
#rownames(salinity.norm)<-rown


#salinity.norm<-na.omit(salinity.norm)
#class(salinity.norm)
# filter after normalization
dim(salinity.norm)
# average signal probes
dim(salinity.norm[rowMeans(salinity.norm)>1,])
salinity.norm.filtered<-salinity.norm[rowMeans(salinity.norm)>1,]
norm.data.filtered<-salinity.norm.filtered
head(norm.data.filtered)

annot.keep<-subset(annot.keep,rownames(annot.keep) %in% rownames(norm.data.filtered))
# unfiltered:
#annot.keep<-subset(annot.keep,rownames(annot.keep) %in% rownames(salinity.norm))
an<-new("AnnotatedDataFrame",data=annot.keep)
experimentData<-new("MIAME",name="Sara Edge, Lisa Cohen, Ana María González Angel",
                    lab="Marine Genomics, HBOI",
                    contact="lisa.johnson.cohen@gmail.com",
                    title="Porites astreoides larvae actuve salinity exposure: 25 ppt, 30 ppt, 35 ppt",
                    abstract="ExpressionSet for P. astreoides Microarray Data",
                    other=list(notes="Extracted from raw data with GenePixPro software"))
# no outliers removed:
sample_type<-c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt","sal35ppt")
# remove outliers
#sample_type<-c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt")
# remove: ,,,"c253295110048_ar4_jun1912_C_SL",,"c253295110054_ar2_jul1112_C_C"
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
salinity_ExpressionSet_filtered<-new("ExpressionSet",exprs=norm.data.filtered,phenoData=pd,experimentData=experimentData,featureData=an)
#arrayQualityMetrics(expressionset = salinity_ExpressionSet_filtered,outdir = "salinity_ArrayQualityMetrics_log_collapsed_by_transcriptName_loessnorm_filtered_one35pptoutlierremoved",force = TRUE,intgroup = c("sample_type"),do.logtransform=FALSE)
###
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3)))
colnames(design) <- c("sal35ppt", "sal30ppt", "sal25ppt")
design
fit <- lmFit(salinity_ExpressionSet_filtered,design)
contrast.matrix <- makeContrasts(sal25ppt-sal30ppt, sal25ppt-sal35ppt, levels=design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
head(fit2)
tt_25ppt_v_30ppt<-topTable(fit2,coef="sal25ppt - sal30ppt", n=1000,adjust="BH")
tt_25ppt_v_35ppt<-topTable(fit2,coef="sal25ppt - sal35ppt",n=1000,adjust="BH")
tt_25ppt_v_30ppt<-topTable(fit2,coef="sal25ppt - sal30ppt",n=1000,adjust="BH")
colnames(tt_25ppt_v_30ppt)
norm.data.reps<-cbind(norm.data.filtered,rownames(norm.data.filtered))
colnames(norm.data.reps)<-c("c253295110048_ar2_jun1912_A_SL","c253295110048_ar1_jun1912_B_SL","c253295110048_ar4_jun1912_C_SL",
                            "c253295110048_ar3_jun1912_A_L","c253295110048_ar6_jun1912_B_L","c253295110048_ar5_jun1912_C_L",
                            "c253295110048_ar8_jun1912_A_C","c253295110048_ar7_jun1912_B_C","Name")
colnames(norm.data.reps)
tt_25ppt_v_30ppt_reps<-merge(norm.data.reps,tt_25ppt_v_30ppt,by="Name")
tt_25ppt_v_35ppt_reps<-merge(norm.data.reps,tt_25ppt_v_35ppt,by="Name")
#write.csv(tt_25ppt_v_30ppt_reps,"P_ast_microarray_25pptvs30ppt_reps_byName_quantile.csv")
#write.csv(tt_25ppt_v_35ppt_reps,"P_ast_microarray_25pptvs35ppt_reps_byName_quantile.csv")
# filter by padj<0.25
#overlap
tt_25pptv_30ppt_sig<-subset(tt_25ppt_v_30ppt_reps,tt_25ppt_v_30ppt_reps$adj.P.Val<0.25)
tt_25pptv_35ppt_sig<-subset(tt_25ppt_v_35ppt_reps,tt_25ppt_v_35ppt_reps$adj.P.Val<0.25)
tt_25pptv_30ppt_sig<-tt_25pptv_30ppt_sig[order(-tt_25pptv_30ppt_sig$adj.P.Val),]
tt_25pptv_35ppt_sig<-tt_25pptv_35ppt_sig[order(-tt_25pptv_35ppt_sig$adj.P.Val),]
dim(tt_25pptv_30ppt_sig)
dim(tt_25pptv_35ppt_sig)
colnames(tt_25pptv_35ppt_sig)
dim(subset(tt_25pptv_35ppt_sig,logFC< -1 | logFC>1))
logFC_1<-subset(tt_25pptv_35ppt_sig,logFC< -1 | logFC>1)
head(logFC_1)
#write.csv(logFC_1,"P_ast_salinity_67sig_logFC1.csv")
source('~/Documents/scripts/overLapper.R')
setlist <- list(sal25ppt_v_sal30ppt_padj0.25=as.vector(tt_25pptv_30ppt_sig$Name),sal25ppt_v_sal35ppt_padj0.25=as.vector(tt_25pptv_35ppt_sig$Name))
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