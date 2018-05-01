#### The relationship between gene expression and GC content
TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", header = T, row.names = 1, stringsAsFactors = F)
CellSum <- apply(TPM,1,function(x){sum(x>0)})
GEmean <- apply(TPM,1,function(x){mean(x[x>0])})
sum(CellSum==0)
[1] 10893
GEmean[is.na(GEmean)] <- 0
sum(CellSum>0 & GEmean>0)
[1] 47158
TPM <- TPM[CellSum>0 & GEmean>0,]
CellSum <- apply(TPM,1,function(x){sum(x>0)})
GEmean <- apply(TPM,1,function(x){mean(x[x>0])})

levels <- factor(cut(CellSum,unique(quantile(CellSum, probs=seq(0,1,.01))),include.lowest=T),labels=1:78)
geneLevel <- data.frame(gene = names(CellSum), level = as.character(levels), stringsAsFactors = F)


gtf <- read.table("/mnt/data5/BGI/UCB/tangchao/Homo_sapiens.GRCh38.87.gtf", header = FALSE, sep = "\t", stringsAsFactors = F)
gtf <- gtf[gtf$V3=="gene",]
test <- as.data.frame(do.call(rbind, strsplit(as.character(gtf$V9), split = "[; ]")))

valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),  
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})), stringsAsFactors=F)

gene_gtf <- cbind(gtf[,1:8],gtf_infor)

GLI <- merge(x = gene_gtf, y = geneLevel, by.x = "gene_id", by.y = "gene", all.y = T)
GLI <- GLI[GLI$V1 %in% c(1:22,"X","Y"),]

library(GenomicRanges)
GR <- GRanges(seqnames = paste("chr", GLI$V1,sep = ""),IRanges(start = GLI$V4, end = GLI$V5), strand = GLI$V7)


library(BSgenome)
available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))
library(BSgenome.Hsapiens.UCSC.hg38)

GR_dna <- getSeq(Hsapiens, GR)
MY_GC <- letterFrequency(GR_dna, "GC", as.prob=TRUE)

GLI$GC <- as.numeric(MY_GC)

GLI$level <- factor(GLI$level, levels = 1:78)
aggregate(GLI$GC,by=list(GLI$level),mean)
aggregate(GLI$GC,by=list(GLI$gene_biotype),mean)


aggregate(GLI$GC,by=list(GLI$level),boxplot)
library(ggplot2)
ggplot(GLI, aes(x=level, y=GC)) + 
    geom_boxplot(outlier.shape=NA)



biomart_GC <- read.table("/mnt/data5/BGI/UCB/tangchao/GEGC/Gene_GC_content.txt", sep = "\t", header = T)
colnames(biomart_GC) <- c("gene","GC_biomart")
GLI <- merge(x = GLI, y = biomart_GC, by.x = "gene_id", by.y = "gene", all.x = T)

aggregate(GLI$GC,by=list(GLI$level),mean)
aggregate(GLI$GC_biomart,by=list(GLI$level),function(x){mean(x, na.rm = T)})







