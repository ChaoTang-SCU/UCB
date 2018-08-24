#### CD45 Version7: use exon expression as cutoff and splice junction for validate
library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "median",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_exp/")
system.time(median_merge<-multmerge(path))
median_merge <- as.data.frame(median_merge)
colnames(median_merge) <- substr(colnames(median_merge),1,10)

colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB3", replacement="UCB4")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB1", replacement="UCB3")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB5", replacement="UCB1")

row.names(median_merge) <- paste("E", 3:7, sep = "")
median_merge <- t(median_merge)
median_merge <- data.frame(median_merge)
median_merge$Exon_Sum <- as.numeric(rowSums(median_merge))

hist(log10(median_merge$Exon_Sum), breaks=200)

median_merge$difE3E7 <- median_merge$E3/(median_merge$E3+median_merge$E7)
median_merge[is.na(median_merge)] <- 0

hist(median_merge$difE3E7, breaks=200)

hist(log10(median_merge$Exon_Sum), breaks=200)
abline(v = 1.3, col = 2)
sum(log10(median_merge$Exon_Sum)>=1.3)
[1] 2663

plot(density(median_merge$difE3E7, adjust=.5))
abline(v = .4, col = 2)
abline(v = .75, col = 2)
sum(median_merge$difE3E7>=.4&median_merge$difE3E7<=.75)
[1] 2364

sum(log10(median_merge$Exon_Sum)>=1.3 & median_merge$difE3E7>=.4 & median_merge$difE3E7<=.75)
[1] 2109

exon_exp_tu <- median_merge[log10(median_merge$Exon_Sum)>=1.3 & median_merge$difE3E7>=.4 & median_merge$difE3E7<=.75, ]
hist(log10(exon_exp_tu$Exon_Sum), breaks=200)
plot(density(exon_exp_tu$difE3E7, adjust=.5))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Exon_expression_cutoff.pdf")
hist(log10(median_merge$Exon_Sum), breaks=200, xlab = "log10(Read counts of CD45 exon3 to exon7)", main = "")
abline(v = 1.3, col = 2)
plot(density(median_merge$difE3E7, adjust=.5), xlab = "exon3/(exon3+exon7)", main = "")
abline(v = .4, col = 2)
abline(v = .75, col = 2)
dev.off()


#### CD45 of salmon
library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(paste(filenames,"quant.sf",sep = "/"), function(x){
    tmp<-fread(x, header=TRUE, select = 5, sep = "\t")
    setnames(tmp, "NumReads",strsplit(x,"/")[[1]][9])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/xiaqiuqi/result/3039_cd45_result")

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")
colnames(cov_merge) <- substr(colnames(cov_merge), nchar(colnames(cov_merge))-9, nchar(colnames(cov_merge)))

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_merge[,row.names(tsne)] -> cov_tu

#cov_tu[cov_tu<20] <- 0

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

sum(colnames(cov_tu) %in% row.names(exon_exp_tu))
# [1] 1784
1784/3039
# [1] 0.5870352

#### Cutoff use exon expression
cov_tu[, !colnames(cov_tu) %in% row.names(exon_exp_tu)] <- 0

#### validate use splice junction

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
SJ34 = "1:198692374-198696711"
SJ35 = "1:198692374-198699563"
SJ45 = "1:198696910-198699563"
SJ46 = "1:198696910-198702386"
SJ56 = "1:198699705-198702386"
SJ57 = "1:198699705-198703297"
SJ67 = "1:198702531-198703297"
SJ47 = "1:198696910-198703297"
SJ36 = "1:198692374-198702386"
SJ37 = "1:198692374-198703297"

as.data.frame(te[.(c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37))])[,-1] -> CD45_SJ
colnames(CD45_SJ) <- substr(colnames(CD45_SJ),1,10)
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

CD45_SJ <- CD45_SJ[,colnames(cov_tu)]

cov_tu["RA", apply(CD45_SJ[c("SJ34","SJ47"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAB", apply(CD45_SJ[c("SJ34","SJ45","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RABC", apply(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAC", apply(CD45_SJ[c("SJ34","SJ46","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RB", apply(CD45_SJ[c("SJ35","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RBC", apply(CD45_SJ[c("SJ35","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RC", apply(CD45_SJ[c("SJ36","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RO", apply(CD45_SJ[c("SJ37"),], 2, function(x) sum(is.na(x))!=0)] <- 0

cov_tu[cov_tu<20] <- 0
table(colSums(cov_tu>0))
#    0    1    2    3    4    5
# 1284  795  536  188  225   11
cov_tu[cov_tu<50] <- 0
table(colSums(cov_tu>0))
#    0    1    2    3    4    5
# 1343  940  437  142  173    4

save(cov_tu, file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")

library(ggplot2)
t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_RABC_v7.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()

test <- cov_tu[c(1,4,7),]
test2 <- data.frame(colSums(test))
test3 <- do.call(rbind,apply(test,2,function(x) row.names(test)[which(x>0)]))
test4 <- merge(test2,test3,by=0, all.x=T)
test4[,3] <- as.numeric(test4[,3])
test4[is.na(test4)] <- 0

test4[order(test4[,2]),] -> test4

f9 <- ggplot(data = tsne[test4[,1],], aes(x = tSNE_1, y = tSNE_2, shape = as.factor(test4[,3])))+
  geom_point(aes(colour = log10(as.numeric(test4[,2]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_shape(solid = TRUE) +
  ggtitle("")+
  scale_shape_manual(values = c(19,4,5,6),labels=c("Other","RA","RAC","RC"))+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_RABC_v7_f5.pdf", height = 6, width = 12)
library(cowplot)
plot_grid(f2, f3, f5, f6, f8, f9, nrow = 2)
dev.off()



#### Only plot the major isoform


for(i in 1:ncol(cov_tu)){
  cov_tu[,i][-which.max(cov_tu[,i])] <- 0
}


t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = 1)+
  scale_colour_gradient(high = 'Grey',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_RABC_v7_2.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()



test <- cov_tu[c(1,4,7),]
test2 <- data.frame(colSums(test))
test3 <- do.call(rbind,apply(test,2,function(x) row.names(test)[which(x>0)]))
test4 <- merge(test2,test3,by=0, all.x=T)
test4[,3] <- as.numeric(test4[,3])
test4[is.na(test4)] <- 0

test4[order(test4[,2]),] -> test4

f9 <- ggplot(data = tsne[test4[,1],], aes(x = tSNE_1, y = tSNE_2, shape = as.factor(test4[,3])))+
  geom_point(aes(colour = log10(as.numeric(test4[,2]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_shape(solid = TRUE) +
  ggtitle("")+
  scale_shape_manual(values = c(19,4,6),labels=c("Other","RA","RC"))+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_RABC_v7_2_f5.pdf", height = 6, width = 12)
library(cowplot)
plot_grid(f2, f3, f5, f6, f8, f9, nrow = 2)
dev.off()




cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[, !colnames(cov_tu) %in% row.names(exon_exp_tu)] <- 0
cov_tu["RA", apply(CD45_SJ[c("SJ34","SJ47"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAB", apply(CD45_SJ[c("SJ34","SJ45","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RABC", apply(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAC", apply(CD45_SJ[c("SJ34","SJ46","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RB", apply(CD45_SJ[c("SJ35","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RBC", apply(CD45_SJ[c("SJ35","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RC", apply(CD45_SJ[c("SJ36","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RO", apply(CD45_SJ[c("SJ37"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu[cov_tu<50] <- 0
table(colSums(cov_tu>0))
#    0    1    2    3    4    5
# 1343  940  437  142  173    4



#### Only plot the cells have only one isoform
for(i in 1:ncol(cov_tu)){
  if(sum(cov_tu[,i] > 0 )>1){
    cov_tu[,i] <- 0
  } 
}


t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = 1)+
  scale_colour_gradient(high = 'Grey',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = 1)+
  scale_colour_gradient(high = 'Grey',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_RABC_v7_3.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()




cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[, !colnames(cov_tu) %in% row.names(exon_exp_tu)] <- 0
cov_tu["RA", apply(CD45_SJ[c("SJ34","SJ47"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAB", apply(CD45_SJ[c("SJ34","SJ45","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RABC", apply(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAC", apply(CD45_SJ[c("SJ34","SJ46","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RB", apply(CD45_SJ[c("SJ35","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RBC", apply(CD45_SJ[c("SJ35","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RC", apply(CD45_SJ[c("SJ36","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RO", apply(CD45_SJ[c("SJ37"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu[cov_tu<50] <- 0
table(colSums(cov_tu>0))


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of cell isoform numbers.pdf")
x = barplot(table(apply(cov_tu,2,function(x) sum(x>0))), ylim = c(0,1500),
  ylab = "No. Cells", col = "white", xlab = "No. Isoforms")
n <- table(apply(cov_tu,2,function(x) sum(x>0)))
text(x = x, y = n, labels = n, cex = .8, pos = 3)
dev.off()



cov_tu[,row.names(Cell_type)] -> cov_tu
cov_tu_p <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu_p[is.na(cov_tu_p)] <- 0
major_isoform_p <- apply(cov_tu_p, 2, max)
identical(row.names(Cell_type), names(major_isoform_p))
#[1] TRUE
Cell_type$P_major <- major_isoform_p
Cell_type$P_major2 <- Cell_type$P_major
Cell_type[Cell_type$P_major2==0,"P_major2"] <-NA


ggplot(data = Cell_type, aes(x = major_isoform_p, colour = CellType))+
    geom_density()

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Density plot of major isoform percentage.pdf")
ggplot(data = Cell_type, aes(x = CellType, y = P_major, fill = CellType))+
    geom_violin(adjust = .2, show.legend = FALSE)+
    ylab("Percentage of major isoform")

ggplot(data = Cell_type, aes(x = CellType, y = P_major2, fill = CellType))+
    geom_violin(adjust = 1, show.legend = FALSE)+
    ylab("Percentage of major isoform")
dev.off()




Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
#                   0   1   2   3   4   5
#  B Cells        329 207  49  12  10   0
#  CD4+ T Cells   684 508 275  97 111   3
#  CD8+ T Cells   240 173  93  30  42   1
#  Megakaryocytes  11   8   2   1   2   0
#  Monocytes       33  13   1   2   3   0
#  NK cells        46  31  17   0   5   0

#### Bar chart of cell isoform numbers

Cell_type$Isoforms <- as.factor(Cell_type$Isoforms)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of cell isoform numbers.pdf")
ggplot(Cell_type, aes(Isoforms))+ 
	geom_bar(aes(fill = CellType))+
	geom_text(data = data.frame(table(Cell_type$Isoforms)), aes(x=Var1,y=Freq, label=Freq), nudge_y=20)
dev.off()





#### Bar chart of cell type isoform numbers



pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of cell type isoform numbers2.pdf")
ggplot(Cell_type, aes(CellType))+ 
	geom_bar(aes(fill = Isoforms))+
	geom_text(data = data.frame(table(Cell_type$CellType)), aes(x=Var1,y=Freq, label=Freq), nudge_y=20)

ggplot(Cell_type, aes(CellType))+ 
	geom_bar(aes(fill = Isoforms), position = "fill")+
	scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	ylab("Fraction")
dev.off()



#### Only care about the major isoform


for(i in 1:ncol(cov_tu)){
  cov_tu[,i][-which.max(cov_tu[,i])] <- 0
}
cov_tu <- cov_tu[,colSums(cov_tu)!=0]

to_plot <- data.frame(apply(cov_tu,2,function(x) row.names(cov_tu)[x>0]))
colnames(to_plot) <- "CD45_Type"

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

plot_tab <- cbind(Cell_type[row.names(to_plot),],to_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of cell type major isoform fraction.pdf")
ggplot(plot_tab, aes(CellType))+ 
	geom_bar(aes(fill = CD45_Type), position = "fill")+
	scale_x_discrete(labels=paste(levels(plot_tab$CellType), "(", table(plot_tab$CellType), ")", sep = ""))+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	ylab("Fraction")
dev.off()


#save(cov_tu, file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")


load("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")
cov_tu_2 <- cov_tu[,colSums(cov_tu>0)==2]
apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))
table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")))
# RABC|RBC  RABC|RC  RABC|RO RAB|RABC   RAB|RB   RAB|RO  RAC|RBC   RA|RAB
#      134        1        4       72       36        2        1        1
#   RBC|RO   RB|RBC    RB|RC    RB|RO
#        6      160        1       19
test <- table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Distribution of isoform_2.pdf")
par(mar=c(10.1, 4.1, 4.1, 4.1))
barplot(sort(test, decreasing=T), las=2, ylab = "Fequency")
dev.off()


to_plot <- data.frame(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")))
colnames(to_plot) <- "CD45_Type"

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

plot_tab <- cbind(Cell_type[row.names(to_plot),],to_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of cell type major isoform_2 fraction.pdf")
ggplot(plot_tab, aes(CellType))+ 
	geom_bar(aes(fill = CD45_Type), position = "fill")+
	scale_x_discrete(labels=paste(levels(plot_tab$CellType), "(", table(plot_tab$CellType), ")", sep = ""))+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	ylab("Fraction")
dev.off()




to_plot <- data.frame(colSums(cov_tu>0))
colnames(to_plot) <- "Isoforms"

plot_tab <- cbind(median_merge[row.names(to_plot),], to_plot)
plot_tab$Isoforms <- factor(plot_tab$Isoforms)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Box plot of exon expression of different isoform numbers.pdf")
ggplot(plot_tab, aes(x = Isoforms, y = log10(Exon_Sum+1)))+
    geom_boxplot()
dev.off()



#### bulk CD45 Sashimi plot



path<-file.path("/mnt/data5/BGI/UCB/xiaqiuqi/result/cell_20_samples_results_files")

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(paste(filenames,"quant.sf",sep = "/"), function(x){
    tmp<-fread(x, header=TRUE, select = 5, sep = "\t")
    setnames(tmp, "NumReads",strsplit(x,"/")[[1]][9])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

system.time(bulk_merge<-multmerge(path))

bulk_merge <- as.data.frame(bulk_merge)
row.names(bulk_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

mon_merge <- bulk_merge[,grep("mono", colnames(bulk_merge))]
tcl_merge <- bulk_merge[,grep("mono", colnames(bulk_merge), invert=T)]

tcl_merge[tcl_merge<50] <- 0

rowSums(tcl_merge>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   0   10   10    0    9   10    0    2

mon_merge[mon_merge<50] <- 0
rowSums(mon_merge>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   0   10   10    0   10   10    0   10


tcl_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD4+ T Cells"]]), tcl_merge)
tcl_merge <- apply(tcl_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of cd4 sc and bulk.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(tcl_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,10)))
legend("topright", inset=c(-.2,0), legend = rownames(tcl_merge),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()



mon_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Monocytes"]]), mon_merge)
mon_merge <- apply(mon_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bar chart of Monocytes sc and bulk.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(mon_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,10)))
legend("topright", inset=c(-.2,0), legend = rownames(mon_merge),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()



list.files("/mnt/data5/BGI/UCB/xiaqiuqi/result/cell_20_samples_results_files")
tcel <- list.files("/mnt/data5/BGI/UCB/xiaqiuqi/result/cell_20_samples_results_files")[grep("mono", list.files("/mnt/data5/BGI/UCB/xiaqiuqi/result/cell_20_samples_results_files"), invert=T)]
mono <- list.files("/mnt/data5/BGI/UCB/xiaqiuqi/result/cell_20_samples_results_files")[grep("mono", list.files("/mnt/data5/BGI/UCB/xiaqiuqi/result/cell_20_samples_results_files"), invert=F)]

tcel <- substr(tcel,1,15)
mono <- substr(mono,1,15)

Tcel_tsv <- read.table("/mnt/nfs_nas/Lab/People/tangchao/Project/IR_QTL/data/phenotype/BAM_tsv/Tcells.tsv", header = F, stringsAsFactors = F, sep = "\t")
Mono_tsv <- read.table("/mnt/nfs_nas/Lab/People/tangchao/Project/IR_QTL/data/phenotype/BAM_tsv/Monocytes.tsv", header = F, stringsAsFactors = F, sep = "\t")

mapply(function(x) grep(x,Mono_tsv$V2), mono)
write.table(Mono_tsv[mapply(function(x) grep(x,Mono_tsv$V2), mono),], file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Monocytes_sashimiplot_bam.tsv", sep = "\t", row.names = F, quote = F, col.names = F)

python /mnt/data5/BGI/UCB/tangchao/soft/ggsashimi/sashimi-plot.py \
                -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Monocytes_sashimiplot_bam.tsv \
                -c chr1:198660476-198673501 \
                -g /mnt/nfs_nas/Lab/People/tangchao/Project/IR_QTL/data/software/ggsashimi/gencode.v15.annotation.transcript.exon.gtf \
                -M 11 \
                --alpha 0.3 \
                --base-size=20 \
                --ann-height=12 \
                --height=4 \
                --width=18 \
                -o  /mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bulk_Monocytes.pdf \
                -C 3



list.files("/mnt/nfs_nas/EGAD00001002671/EGAD00001002671_bam_decrypted/", pattern = "bam$")
tcel_bam <- list.files("/mnt/nfs_nas/EGAD00001002671/EGAD00001002671_bam_decrypted", pattern = "bam$", full.names = T)[grep("mRNA", list.files("/mnt/nfs_nas/EGAD00001002671/EGAD00001002671_bam_decrypted/", pattern = "bam$"))[1:10]]
data.frame(V1 = tcel, V2 = tcel_bam, V3 = "T Cells")
write.table(data.frame(V1 = tcel, V2 = tcel_bam, V3 = "T Cells"), file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Tcells_sashimiplot_bam.tsv", sep = "\t", row.names = F, quote = F, col.names = F)

python /mnt/data5/BGI/UCB/tangchao/soft/ggsashimi/sashimi-plot.py \
                -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Tcells_sashimiplot_bam.tsv \
                -c chr1:198660476-198673501 \
                -g /mnt/nfs_nas/Lab/People/tangchao/Project/IR_QTL/data/software/ggsashimi/gencode.v15.annotation.transcript.exon.gtf \
                -M 11 \
                --alpha 0.3 \
                --base-size=20 \
                --ann-height=12 \
                --height=4 \
                --width=18 \
                -o  /mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Bulk_CD4TCellss.pdf \
                -C 3




#### CD45 gene expression bias 

library(data.table)
tpm <- fread("/mnt/data5/BGI/UCB/ExpMat_NewID/UCB_rsem_TPM_mat.txt", header = T, sep = "\t")
setkey(tpm,V1)
CD45_tpm <- data.frame(tpm["ENSG00000081237",][,-1])

table(colSums(cov_tu>0))
#    0    1    2    3    4    5
# 1343  940  437  142  173    4
cell_iso <- as.data.frame(colSums(cov_tu>0))
colnames(cell_iso) <- "Isoforms"
CD45_tpm[,row.names(cell_iso)]
cell_iso <- cbind(cell_iso, t(CD45_tpm[,row.names(cell_iso)]))
colnames(cell_iso)[2] <- "TPM"
cell_iso$Isoforms <- factor(cell_iso$Isoforms)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_gene_expression_of_different_sioforms.pdf")
ggplot(cell_iso, aes(x = Isoforms, y = log10(TPM+1)))+
    geom_boxplot()
dev.off()



#### sequence depth bias 

sta <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", header = T, sep = "\t", row.names = 1)
cell_iso$Input_reads <- sta[row.names(cell_iso), "Input_reads"]
cell_iso$Mapped_reads <- sta[row.names(cell_iso), "Mapped_reads"]

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/sequence_depth_of_different_sioforms.pdf")
ggplot(cell_iso, aes(x = Isoforms, y = log10(Input_reads+1)))+
    geom_boxplot()
ggplot(cell_iso, aes(x = Isoforms, y = log10(Mapped_reads+1)))+
    geom_boxplot()
dev.off()



#### 5' or 3' bias

test <- fread("/mnt/data5/BGI/UCB/tangchao/BCseq/featurecount_ensembl/UCB3.00575.epCount", skip = 1, header = T, sep = "\t")
setkey(test,Geneid)

test2 <- test["ENSG00000081237",]

gtf <- fread("/mnt/data5/BGI/UCB/tangchao/Homo_sapiens.GRCh38.87.transccript.exon.gtf", header = F, sep = "\t")
grep("ENST00000442510", gtf[[9]])

gtf <- gtf[grep("ENST00000442510", gtf[[9]]),]
gtf <- gtf[-1, 1:7]

test2 <- as.data.frame(test2)
colnames(test2)[7] <- strsplit(colnames(test2)[7], "[/_]")[[1]][10]

sum(test2$Start %in% gtf[[4]] & test2$End %in% gtf[[5]])
[1] 33
test2 <- test2[test2$Start %in% gtf[[4]] & test2$End %in% gtf[[5]], ]
test2 <- test2[order(test2$Start, test2$End), ]

files <- list.files("/mnt/data5/BGI/UCB/tangchao/BCseq/featurecount_ensembl", pattern = "epCount$", full.names=T)

mclapply(files, function(x){
  test <- fread(x, skip = 1, header = T, sep = "\t")
  setkey(test,Geneid)
  test2 <- test["ENSG00000081237",]
  test2 <- as.data.frame(test2)
  colnames(test2)[7] <- strsplit(colnames(test2)[7], "[/_]")[[1]][10]
  test2 <- test2[test2$Start %in% gtf[[4]] & test2$End %in% gtf[[5]], ]
  test2 <- test2[order(test2$Start, test2$End), ]
  return(test2)
  }, mc.cores = 10) -> tang


tang <- list()
for(i in 1:length(files)){
  print(i)
  test <- fread(files[i], skip = 1, header = T, sep = "\t")
  setkey(test,Geneid)
  test2 <- test["ENSG00000081237",]
  test2 <- as.data.frame(test2)
  colnames(test2)[7] <- strsplit(colnames(test2)[7], "[/_]")[[1]][10]
  test2 <- test2[test2$Start %in% gtf[[4]] & test2$End %in% gtf[[5]], ]
  test2 <- test2[order(test2$Start, test2$End), ]
  tang[[i]] <- test2
}

save(tang, file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/exon_exp_of_CD45.RData")

test <- do.call(rbind,lapply(tang, function(x) as.data.frame(t(x[,7]))))
row.names(test) <- unlist(lapply(tang, function(x) colnames(x)[7]))
colnames(test) <- paste("exon", 1:ncol(test), sep = "")

#for(i in 1:nrow(test)){
#  print(i)
#  test[i, (test[i,]/sum(test[i,])) < 0.05] <- 0
#}

test[test<10] <- 0



##Mfuzz#########
library(Mfuzz)
library(reshape2)
library(tidyverse)
library(marray)
##### A slight change for plot function#####
cos_fuc <- function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels,
                     time.points, ylim.set = c(0, 0), xlab = "Time", ylab = "Expression changes",
                     x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black",
                     col.lab = "black", col.main = "black", col.sub = "black",
                     col = "black", centre = FALSE, centre.col = "black", centre.lwd = 2,
                     Xwidth = 5, Xheight = 5, single = FALSE,clu = clu, ...)
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(exprs(eset))[[1]])
  if (missing(colo)) {
    colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800",
              "#FF6000", "#FF7800", "#FF8F00", "#FFA700", "#FFBF00",
              "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00",
              "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28",
              "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F",
              "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF",
              "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF",
              "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF",
              "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF",
              "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F",
              "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }
  else {
    if (colo == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))),
                     c(255:150))
      colo <- rgb(b = fancy.blue/255, g = fancy.green/255,
                  r = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  for (j in 1:dim(cl[[1]])[[1]]) {
    if (single)
      j <- single
    tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
    tmpmem <- memship[clusterindex == j, j]
    if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
      if (x11)
        X11(width = Xwidth, height = Xheight)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      if (sum(ylim.set == c(0, 0)) == 2) {
        ylim <- c(ymin, ymax)
      }
      else {
        ylim <- ylim.set
      }
      if (!is.na(sum(mfrow))) {
        par(mfrow = mfrow, bg = bg, col.axis = col.axis,
            col.lab = col.lab, col.main = col.main, col.sub = col.sub,
            col = col)
      }
      else {
        par(bg = bg, col.axis = col.axis, col.lab = col.lab,
            col.main = col.main, col.sub = col.sub, col = col)
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points)))
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim,
                   xlab = xlab, ylab = ylab, main = paste("Cluster",
                                                          j,'(',clu[j,],')'), axes
                   = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]),
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points,
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]], time.labels,
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col,
             ...)
        axis(2, col = ax.col, ...)
      }
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      if (sum(ylim.set == c(0, 0)) == 2) {
        ylim <- c(ymin, ymax)
      }
      else {
        ylim <- ylim.set
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points)))
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim,
                   xlab = xlab, ylab = ylab, main = paste("Cluster",
                                                          j,'(',clu[j,],')'), axes
                   = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]),
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points,
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]], time.labels,
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col,
             ...)
        axis(2, col = ax.col, ...)
      }
    }
    if (length(tmpmem) > 0) {
      for (jj in 1:(length(colorseq) - 1)) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <=
                     colorseq[jj + 1])
        if (sum(tmpcol) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            if (missing(time.points)) {
              lines(tmp[tmpind[k], ], col = colo[jj])
            }
            else lines(time.points, tmp[tmpind[k], ],
                       col = colo[jj])
          }
        }
      }
    }
    if (centre) {
      lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
    }
    if (single)
      return()
  }
}

cos_bar <- function (x, horizontal = TRUE, cex.axis = cex.axis,col = heat.colors(50), scale = 1:length(x), k = 10, ...)
{
  if (is.numeric(x)) {
    x <- x
    colmap <- col
  }
  else {
    colmap <- x
    low <- range(scale)[1]
    high <- range(scale)[2]
    x <- seq(low, high, length = length(x))
  }
  if (length(x) > k)
    x.small <- seq(x[1], x[length(x)], length = k)
  else x.small <- x
  if (horizontal) {
    image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "",
          ylab = "", col = colmap, ...)
    axis(1, at = rev(x.small),cex.axis = cex.axis, labels = signif(rev(x.small),
                                                                   2), srt = 270)
  }
  if (!horizontal) {
    image(1, x, matrix(x, 1, length(x)), axes = FALSE, xlab = "",
          ylab = "", col = colmap, ...)
    par(las = 1)
    axis(4, at = rev(x.small),cex.axis = cex.axis, labels = signif(rev(x.small),
                                                                   2))
    par(las = 0)
  }
  box()
}



library(Biobase)
eset <- new('ExpressionSet', exprs=log10(as.matrix(test)+1))
pdf("~/test.pdf")
eset <- filter.std(eset,min.std = 0.001)
dev.off()
eset.s <- Mfuzz::standardise(eset)
cl <- mfuzz(eset.s,c = 10, m = mestimate(eset.s)[1])
clu <- as.data.frame(table(as.data.frame(cl$cluster)[,1]))
clu$Var1 <- NULL


col <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700",
         "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00",
         "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
         "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40",
         "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7",
         "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
         "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF",
         "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF",
         "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF",
         "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
         "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078",
         "#FF0060", "#FF0048", "#FF0030", "#FF0018")

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Mfuzz_test.pdf")
pdf('/mnt/data2/xiaoxia/APA_PDUI_V2/figure/test3.pdf',25,26)
par(mar = c(8,6,3,0.5))
cos_fuc(eset.s,cex.lab = 2.5,cex.axis = 2,cex.main = 3,cl=cl,clu = clu,xlab = '',mfrow=c(6,6),
  time.labels=1:ncol(test),x11 = F)
par(mar= c(15,2,12,0))
cos_bar(seq(0,1,.001),k=11, col = col, horizontal = TRUE,cex.axis = 2)
dev.off()



eset <- new('ExpressionSet', exprs=log10(as.matrix(test)+1))
pdf("~/test.pdf")
eset <- filter.std(eset,min.std = 0.001)
dev.off()
eset.s <- Mfuzz::standardise(eset)
cl <- mfuzz(eset.s,c = 12, m = mestimate(eset.s)[1])
clu <- as.data.frame(table(as.data.frame(cl$cluster)[,1]))
clu$Var1 <- NULL


col <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700",
         "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00",
         "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
         "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40",
         "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7",
         "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
         "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF",
         "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF",
         "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF",
         "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
         "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078",
         "#FF0060", "#FF0048", "#FF0030", "#FF0018")

#pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/Mfuzz_test.pdf")
pdf('/mnt/data2/xiaoxia/APA_PDUI_V2/figure/test3.pdf',25,26)
par(mar = c(8,6,3,0.5))
cos_fuc(eset.s,cex.lab = 2.5,cex.axis = 2,cex.main = 3,cl=cl,clu = clu,xlab = '',mfrow=c(6,6),
  time.labels=1:ncol(test),x11 = F)
par(mar= c(15,2,12,0))
cos_bar(seq(0,1,.001),k=11, col = col, horizontal = TRUE,cex.axis = 2)
dev.off()


names(cl$cluster[cl$cluster == 6])
names(cl$cluster[cl$cluster == 7])

sashimi(junction = "1:198692374-198703297", cell = sample(names(cl$cluster[cl$cluster == 2]),20))
sashimi(junction = "1:198692374-198703297", cell = sample(names(cl$cluster[cl$cluster == 6]),20))
sashimi(junction = "1:198692374-198703297", cell = sample(names(cl$cluster[cl$cluster == 9]),20))
sashimi(junction = "1:198692374-198703297", cell = sample(names(cl$cluster[cl$cluster == 11]),20))




test <- do.call(rbind,lapply(tang, function(x) as.data.frame(t(x[,7]/x[,6]))))

row.names(test) <- unlist(lapply(tang, function(x) colnames(x)[7]))
colnames(test) <- paste("exon", 1:ncol(test), sep = "")


library(Biobase)
eset <- new('ExpressionSet', exprs=log10(as.matrix(test)+1))
pdf("~/test.pdf")
eset <- filter.std(eset,min.std = 0.001)
dev.off()
eset.s <- Mfuzz::standardise(eset)
cl <- mfuzz(eset.s,c = 10, m = mestimate(eset.s)[1])
clu <- as.data.frame(table(as.data.frame(cl$cluster)[,1]))
clu$Var1 <- NULL

pdf('/mnt/data2/xiaoxia/APA_PDUI_V2/figure/test3.pdf',25,26)
par(mar = c(8,6,3,0.5))
cos_fuc(eset.s,cex.lab = 2.5,cex.axis = 2,cex.main = 3,cl=cl,clu = clu,xlab = '',mfrow=c(6,6),
  time.labels=1:ncol(test),x11 = F)
par(mar= c(15,2,12,0))
cos_bar(seq(0,1,.001),k=11, col = col, horizontal = TRUE,cex.axis = 2)
dev.off()





#### RO vs. RABC ================================================================================================================================



load(file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")
sum(cov_tu[8,]>0)
#[1] 111
cov_tu <- cov_tu[, cov_tu[8,]>0]
cov_tu_p <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu_p[is.na(cov_tu_p)] <- 0

plot(x=cov_tu_p[8,], y=colSums(cov_tu_p[grep("A",row.names(cov_tu)),]))
cov_tu_p -> cov_tu_p_salmon
plot(x=cov_tu_p_salmon[8,], y=cov_tu_p_salmon[grep("ABC",row.names(cov_tu_p_salmon)),])
plot(x=cov_tu_p_salmon[8,], y=colSums(cov_tu_p_salmon[grep("A",row.names(cov_tu_p_salmon)),]), 
     xlab = "Percentage of RO", ylab = "Percentage of RA+", main = "Salmon")


cov_tu <- read.table("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2/CD45_RABC.txt", header = T, row.names = 1)
cov_tu <- cov_tu[, cov_tu[1,]>0]
cov_tu_p <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu_p[is.na(cov_tu_p)] <- 0

plot(x=cov_tu_p[1,], y=colSums(cov_tu_p[grep("A",row.names(cov_tu)),]), 
     xlab = "Percentage of RO", ylab = "Percentage of RA+", main = "Splice Junction")

pdf('/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RO_VS_RA+.pdf',width = 12, height = 6)
par(mfrow = c(1,2))
plot(x=cov_tu_p_salmon[8,], y=colSums(cov_tu_p_salmon[grep("A",row.names(cov_tu_p_salmon)),]), 
     xlab = "Percentage of RO", ylab = "Percentage of RA+", main = "Salmon")

plot(x=cov_tu_p[1,], y=colSums(cov_tu_p[grep("A",row.names(cov_tu)),]), 
     xlab = "Percentage of RO", ylab = "Percentage of RA+", main = "Splice Junction")
dev.off()











#### Bulk Cord blood CD45 Salmon ==========================================

load(file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")

all_salmon_quant <- list.files("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/EGA_CDD45_Salmon/transcripts_quant", recursive=T, pattern="quant.sf", full.names=T)

my_list <- list()
for (i in 1:length(all_salmon_quant)) {
  my_list[[i]] <- read.table(all_salmon_quant[[i]], sep = "\t", header = T)[5]
}
cd45_tab <- do.call(cbind, my_list)

colnames(cd45_tab) <- do.call(rbind, strsplit(all_salmon_quant, "/"))[,11]
row.names(cd45_tab) <- c("RA", "RAB", "RABC", "RAC", "RB", "RBC", "RC", "RO")
cd45_tab <- cd45_tab[, order(substr(colnames(cd45_tab),10,100))]
cd45_tab <- t(cd45_tab)
cd45_tab[cd45_tab<50] <- 0

all_info <- read.table("/mnt/raid61/Science_2014_chen/STAR_output/All_Sample_Information.txt", sep = "\t", header = T, stringsAsFactors=F)
TISSUE_TYPE <- data.frame(TISSUE_TYPE=all_info[,"TISSUE_TYPE"], row.names=paste(all_info$ID, all_info$CELL_TYPE_ABB, sep = "_"), stringsAsFactors = F)
TISSUE_TYPE <- data.frame(TISSUE_TYPE = TISSUE_TYPE[row.names(cd45_tab),],row.names=row.names(cd45_tab))

cd45_tab <- cd45_tab[row.names(TISSUE_TYPE)[TISSUE_TYPE$TISSUE_TYPE=="cord blood"],]


#### B Cells ------------------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)
brewer.pal(8, "Set3") -> my_col
names(my_col) <- c("RA", "RAB", "RABC", "RAC", "RB", "RBC", "RC", "RO")

UCB_B <- data.frame(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="B Cells"]])/sum(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="B Cells"]])))
colnames(UCB_B) <- "Single Cell"
EGA_B <- cd45_tab[grep("B_CELL", row.names(cd45_tab)),]
EGA_B <- data.frame(colSums(EGA_B)/sum(colSums(EGA_B)))
colnames(EGA_B) <- "Bulk"
merge(UCB_B, EGA_B, by = 0) -> BCell
library(reshape2)
reshape2::melt(BCell) -> BCell
colnames(BCell) <- c("Isoform", "Method", "Fraction")
BCell$Fraction <- round(BCell$Fraction*1000)

BCell <- do.call(rbind, lapply(split(BCell, 1:nrow(BCell)), function(x) mefa:::rep.data.frame(x, x[3])))
BCell$Isoform <- as.factor(BCell$Isoform)

library(ggplot2)
ggplot(BCell, aes(Method))+ 
  geom_bar(aes(fill = Isoform), position = "fill")+
  #scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("B Cells")+
  ylab("Fraction")+
  xlab("")+
  scale_fill_manual(values = my_col[levels(BCell$Isoform)]) -> p_bcell






#### CD4 T Cells ------------------------------------------------------------------------------------------------------------------------------------

UCB_CD4 <- data.frame(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD4+ T Cells"]])/sum(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD4+ T Cells"]])))
colnames(UCB_CD4) <- "Single Cell"
EGA_CD4 <- cd45_tab[grep("CD4T", row.names(cd45_tab)),]
EGA_CD4 <- data.frame(colSums(EGA_CD4)/sum(colSums(EGA_CD4)))
colnames(EGA_CD4) <- "Bulk"
merge(UCB_CD4, EGA_CD4, by = 0) -> CD4
library(reshape2)
reshape2::melt(CD4) -> CD4
colnames(CD4) <- c("Isoform", "Method", "Fraction")
CD4$Fraction <- round(CD4$Fraction*1000)

CD4 <- do.call(rbind, lapply(split(CD4, 1:nrow(CD4)), function(x) mefa:::rep.data.frame(x, x[3])))
CD4$Isoform <- as.factor(CD4$Isoform)

library(ggplot2)
ggplot(CD4, aes(Method))+ 
  geom_bar(aes(fill = Isoform), position = "fill")+
  #scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Fraction")+
  xlab("")+
  ggtitle("CD4+ T Cells")+
  scale_fill_manual(values = my_col[levels(CD4$Isoform)]) -> p_cd4




#### CD8 T Cells ------------------------------------------------------------------------------------------------------------------------------------


UCB_CD8 <- data.frame(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD8+ T Cells"]])/sum(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD8+ T Cells"]])))
colnames(UCB_CD8) <- "Single Cell"
EGA_CD8 <- cd45_tab[grep("CD8T", row.names(cd45_tab)),]
EGA_CD8 <- data.frame(colSums(EGA_CD8)/sum(colSums(EGA_CD8)))
colnames(EGA_CD8) <- "Bulk"
merge(UCB_CD8, EGA_CD8, by = 0) -> CD8
library(reshape2)
reshape2::melt(CD8) -> CD8
colnames(CD8) <- c("Isoform", "Method", "Fraction")
CD8$Fraction <- round(CD8$Fraction*1000)

CD8 <- do.call(rbind, lapply(split(CD8, 1:nrow(CD8)), function(x) mefa:::rep.data.frame(x, x[3])))
CD8$Isoform <- as.factor(CD8$Isoform)

library(ggplot2)
ggplot(CD8, aes(Method))+ 
  geom_bar(aes(fill = Isoform), position = "fill")+
  #scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Fraction")+
  xlab("")+
  ggtitle("CD8+ T Cells")+
  scale_fill_manual(values = my_col[levels(CD8$Isoform)]) -> p_cd8




#### MK Cells ------------------------------------------------------------------------------------------------------------------------------------


UCB_MK <- data.frame(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Megakaryocytes"]])/sum(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Megakaryocytes"]])))
colnames(UCB_MK) <- "Single Cell"
EGA_MK <- cd45_tab[grep("MK", row.names(cd45_tab)),]
EGA_MK <- data.frame(colSums(EGA_MK)/sum(colSums(EGA_MK)))
colnames(EGA_MK) <- "Bulk"
merge(UCB_MK, EGA_MK, by = 0) -> MK
library(reshape2)
reshape2::melt(MK) -> MK
colnames(MK) <- c("Isoform", "Method", "Fraction")
MK$Fraction <- round(MK$Fraction*1000)

MK <- do.call(rbind, lapply(split(MK, 1:nrow(MK)), function(x) mefa:::rep.data.frame(x, x[3])))
MK$Isoform <- as.factor(MK$Isoform)

library(ggplot2)
ggplot(MK, aes(Method))+ 
  geom_bar(aes(fill = Isoform), position = "fill")+
  #scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Fraction")+
  xlab("")+
  ggtitle("Megakaryocytes")+
  scale_fill_manual(values = my_col[levels(MK$Isoform)]) -> p_mk



#### Monocytes Cells ------------------------------------------------------------------------------------------------------------------------------------


UCB_Monocytes <- data.frame(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Monocytes"]])/sum(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Monocytes"]])))
colnames(UCB_Monocytes) <- "Single Cell"
EGA_Monocytes <- cd45_tab[grep("monocyte", row.names(cd45_tab)),]
EGA_Monocytes <- data.frame(colSums(EGA_Monocytes)/sum(colSums(EGA_Monocytes)))
colnames(EGA_Monocytes) <- "Bulk"
merge(UCB_Monocytes, EGA_Monocytes, by = 0) -> Monocytes
library(reshape2)
reshape2::melt(Monocytes) -> Monocytes
colnames(Monocytes) <- c("Isoform", "Method", "Fraction")
Monocytes$Fraction <- round(Monocytes$Fraction*1000)

Monocytes <- do.call(rbind, lapply(split(Monocytes, 1:nrow(Monocytes)), function(x) mefa:::rep.data.frame(x, x[3])))
Monocytes$Isoform <- factor(Monocytes$Isoform, levels = c("RA", "RAB", "RABC", "RB", "RBC", "RO"))

library(ggplot2)
ggplot(Monocytes, aes(Method))+ 
  geom_bar(aes(fill = Isoform), position = "fill")+
  #scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Fraction")+
  xlab("")+
  ggtitle("Monocytes")+
  scale_fill_manual(values = my_col[levels(Monocytes$Isoform)]) -> p_monocytes




#### NK Cells ------------------------------------------------------------------------------------------------------------------------------------


UCB_NK <- data.frame(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="NK cells"]])/sum(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="NK cells"]])))
colnames(UCB_NK) <- "Single Cell"
EGA_NK <- cd45_tab[grep("NK", row.names(cd45_tab)),]
EGA_NK <- data.frame(colSums(EGA_NK)/sum(colSums(EGA_NK)))
colnames(EGA_NK) <- "Bulk"
merge(UCB_NK, EGA_NK, by = 0) -> NK
library(reshape2)
reshape2::melt(NK) -> NK
colnames(NK) <- c("Isoform", "Method", "Fraction")
NK$Fraction <- round(NK$Fraction*1000)

NK <- do.call(rbind, lapply(split(NK, 1:nrow(NK)), function(x) mefa:::rep.data.frame(x, x[3])))
NK$Isoform <- as.factor(NK$Isoform)

library(ggplot2)
ggplot(NK, aes(Method))+ 
  geom_bar(aes(fill = Isoform), position = "fill")+
  #scale_x_discrete(labels=paste(levels(Cell_type$CellType), "(", table(Cell_type$CellType), ")", sep = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Fraction")+
  xlab("")+
  ggtitle("NK cells")+
  scale_fill_manual(values = my_col[levels(NK$Isoform)]) -> p_nk


pdf('/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/CD45_Salmon_Single_cell_VS_Bulk.pdf',width = 9, height = 6)
library(cowplot)
plot_grid(p_bcell, p_cd4, p_cd8, p_mk, p_monocytes, p_nk, nrow = 2)
dev.off()


