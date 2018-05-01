

#### splicing QTL 
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
require(data.table)  # v1.6.6
require(gdata) 
FastRemoveMissingValues = function(DT) {
  # either of the following for loops
  
  # by name :
  for (j in names(DT))
    set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

te <- as.data.frame(psi_sj_same_start_table)
FastRemoveMissingValues(te)

FastRemoveMissingValues = function(DT) {
  # either of the following for loops
  
  # by name :
  for (j in names(DT))
    set(DT,which(DT[[j]]>0),j,1)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]>0),j,1)
}
FastRemoveMissingValues(te)

TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/UCB_rsem_TPM_mat.txt", header = T, row.names = 1, stringsAsFactors = F)
TPM <- TPM[,colnames(te)]

f0 <- apply(te,1,function(x) sum(x==0)/length(x))
sum(f0>0.1 & f0<0.9)
[1] 7214
summary(f0[f0>0.1 & f0<0.9])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 0.1004  0.7278  0.8281  0.7634  0.8724  0.8998

te_tu <- te[f0>0.1 & f0<0.9,]
f1 <- apply(te_tu,1,function(x) sum(x==0)/length(x))
summary(f1)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 0.1004  0.7278  0.8281  0.7634  0.8724  0.8998
dim(te_tu)
[1] 7214  3574
dim(TPM)
[1] 31907  3574

library(MatrixEQTL)

snps = SlicedData$new()
snps$CreateFromMatrix(as.matrix(te_tu))


gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(TPM))
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);


me <- Matrix_eQTL_engine(snps = snps,
                         gene = gene,
                         cvrt = SlicedData$new(),
                         output_file_name = " ",
                         pvOutputThreshold = 1e-6,
                         useModel = modelANOVA, 
                         verbose = TRUE,
                         pvalue.hist = "qqplot",
                         min.pv.by.genesnp = FALSE,
                         noFDRsaveMemory = FALSE)

save(me, file = "/mnt/data5/BGI/UCB/tangchao/QTL/me_anova.RData")

me <- Matrix_eQTL_engine(snps = snps,
                         gene = gene,
                         cvrt = SlicedData$new(),
                         output_file_name = " ",
                         pvOutputThreshold = 1e-6,
                         useModel = modelLINEAR, 
                         verbose = TRUE,
                         pvalue.hist = "qqplot",
                         min.pv.by.genesnp = FALSE,
                         noFDRsaveMemory = FALSE)
save(me, file = "/mnt/data5/BGI/UCB/tangchao/QTL/me_linear.RData")


dim(me$all$eqtls)
[1] 17669932        5
my_eqtls <- me$all$eqtls

my_eqtls[which.max(my_eqtls$beta),]

head(my_eqtls)
test_1 <- data.frame(AS = as.numeric(te_tu["7:134442771-134445237",]), GE = as.numeric(TPM["ENSG00000085662", ]))
test_1$AS <- as.factor(test_1$AS)
aggregate(test_1[,2],by=list(test_1$AS),mean)
aggregate(test_1[,2],by=list(test_1$AS),median)

library(dplyr)
group_by(test_1, AS) %>% summarize_each(funs(mean), GE)
group_by(test_1, AS) %>% summarize(mean(GE))


for(i in 1:10){
  test_1 <- data.frame(AS = as.numeric(te_tu[my_eqtls$snps[i],]), GE = as.numeric(TPM[my_eqtls$gene[i],]))
  test_1$AS <- as.factor(test_1$AS)	
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/QTL/Figure/", my_eqtls$snps[i], "_", my_eqtls$gene[i], ".pdf", sep = ""))
  print(ggplot(test_1, aes(x=AS, y=GE)) +
          geom_boxplot(fill="cornflowerblue",color="black",outlier.size = 0, outlier.alpha = 0)+
          geom_point(position="jitter", color="blue", alpha=.5)+
          geom_rug(color="black")+
          xlab(label = my_eqtls$snps[2])+
          ylab(label = my_eqtls$gene[2])+
          ggtitle(paste("FDR =",signif(my_eqtls$FDR[i],3), 
          				"; p0 =", paste(round(as.numeric(f1[my_eqtls$snps[i]]),2)), 
          				"; Median =", paste(round(aggregate(test_1[,2],by=list(test_1$AS),median)[,2],2), collapse="/", sep = ""), 
          				"; Mean =", paste(round(aggregate(test_1[,2],by=list(test_1$AS),mean)[,2],2), collapse="/", sep = ""), sep=" ")))
  dev.off()
}


testyn <- vector()
for(i in 1:1000){
  test_1 <- data.frame(AS = as.numeric(te_tu[my_eqtls$snps[i],]), GE = as.numeric(TPM[my_eqtls$gene[i],]))
  test_1$AS <- as.factor(test_1$AS)
  testyn[i] <- diff(aggregate(test_1[,2],by=list(test_1$AS),median)[,2])>0
  }

sum(testyn)
[1] 127

library(parallel)
mcmapply(function(i){
	test_1 <- data.frame(AS = as.numeric(te_tu[my_eqtls$snps[i],]), GE = as.numeric(TPM[my_eqtls$gene[i],]))
  	test_1$AS <- as.factor(test_1$AS)
  	diff(aggregate(test_1[,2],by=list(test_1$AS),median)[,2])>0
	},1:100) -> testyn

my_eqtls[which(testyn),] -> my_eqtls


for(i in 1:16){
  test_1 <- data.frame(AS = as.numeric(te_tu[my_eqtls$snps[i],]), GE = log2(as.numeric(TPM[my_eqtls$gene[i],])+1))
  test_1$AS <- as.factor(test_1$AS)	
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/QTL/Figure/", my_eqtls$snps[i], "_", my_eqtls$gene[i], "10.pdf", sep = ""))
  print(ggplot(test_1, aes(x=AS, y=GE)) +
          geom_boxplot(fill="cornflowerblue",color="black",outlier.size = 0, outlier.alpha = 0)+
          geom_point(position="jitter", color="blue", alpha=.5)+
          geom_rug(color="black")+
          xlab(label = my_eqtls$snps[2])+
          ylab(label = my_eqtls$gene[2])+
          ggtitle(paste("FDR =",signif(my_eqtls$FDR[i],3), 
          				"; p0 =", paste(round(as.numeric(f1[my_eqtls$snps[i]]),2)), 
          				"; Median =", paste(round(aggregate(test_1[,2],by=list(test_1$AS),median)[,2],2), collapse="/", sep = ""), 
          				"; Mean =", paste(round(aggregate(test_1[,2],by=list(test_1$AS),mean)[,2],2), collapse="/", sep = ""), sep=" ")))
  dev.off()
}



















me <- Matrix_eQTL_engine(snps = snps,
                         gene = gene,
                         cvrt = SlicedData$new(),
                         output_file_name = " ",
                         pvOutputThreshold = 1e-6,
                         useModel = modelLINEAR, 
                         verbose = TRUE,
                         pvalue.hist = "qqplot",
                         min.pv.by.genesnp = FALSE,
                         noFDRsaveMemory = FALSE)
save(me, file = "/mnt/data5/BGI/UCB/tangchao/QTL/me_linear.RData")

my_eqtls <- me$all$eqtls
my_eqtls <- my_eqtls[order(my_eqtls$beta,decreasing=T),]
dim(my_eqtls)
[1] 3541327       5

library(parallel)
p0 <- mcmapply(function(i){
		as.numeric(f1[i])
		}, my_eqtls$snps, mc.cores = 2)
my_eqtls$p0 <- p0


mcmapply(function(i){
	test_1 <- data.frame(AS = as.numeric(te_tu[my_eqtls$snps[i],]), GE = as.numeric(TPM[my_eqtls$gene[i],]))
  	test_1$AS <- as.factor(test_1$AS)
  	diff(aggregate(test_1[,2],by=list(test_1$AS),median)[,2])
	},1:nrow(my_eqtls), mc.cores = 4) -> diff_Median

my_eqtls[which(testyn),] -> my_eqtls






for(i in 1:10){
  test_1 <- data.frame(AS = as.numeric(te_tu[my_eqtls$snps[i],]), GE = as.numeric(TPM[my_eqtls$gene[i],]))
  test_1$AS <- as.factor(test_1$AS)	
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/QTL/Figure/", my_eqtls$snps[i], "_", my_eqtls$gene[i], ".pdf", sep = ""))
  print(ggplot(test_1, aes(x=AS, y=GE)) +
          geom_boxplot(fill="cornflowerblue",color="black",outlier.size = 0, outlier.alpha = 0)+
          geom_point(position="jitter", color="blue", alpha=.5)+
          geom_rug(color="black")+
          xlab(label = my_eqtls$snps[2])+
          ylab(label = my_eqtls$gene[2])+
          ggtitle(paste("FDR =",signif(my_eqtls$FDR[i],3), 
          				"; p0 =", paste(round(as.numeric(f1[my_eqtls$snps[i]]),2)), 
          				"; Median =", paste(round(aggregate(test_1[,2],by=list(test_1$AS),median)[,2],2), collapse="/", sep = ""), 
          				"; Mean =", paste(round(aggregate(test_1[,2],by=list(test_1$AS),mean)[,2],2), collapse="/", sep = ""), sep=" ")))
  dev.off()
}



library(ggplot2)

pdf("/mnt/data5/BGI/UCB/tangchao/QTL/Figure/QQplot.pdf")
plot(me)
dev.off()





