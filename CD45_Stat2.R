#### CD45 Stat version2
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

cov_merge[cov_merge<50] <- 0

cov_merge[,row.names(tsne)] -> cov_tu

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

apply(cov_tu[,row.names(Cell_type)], 2, function(x))

Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
#                   0   1   2   3   4   5   6   7
#  B Cells        240 211  88  47  20   1   0   0
#  CD4+ T Cells   524 549 326 125 113  31   8   2
#  CD8+ T Cells   175 196 105  44  45  12   2   0
#  Megakaryocytes   9  10   2   1   1   1   0   0
#  Monocytes       26  14   5   4   2   1   0   0
#  NK cells        32  33  22   4   6   1   1   0

cell_iso <- table(Cell_type$CellType, Cell_type$Isoforms)


#### Bar chart of cell isoform numbers


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell isoform numbers.pdf")
par(mar=c(10.1, 4.1, 4.1, 2.1))
x <- barplot(cell_iso, col = RColorBrewer::brewer.pal(n = 6, name = "Pastel2"), ylim = c(0,1200),
			 ylab = "No. Cells", xlab = "No. Isoforms", names.arg = colnames(cell_iso))
legend("topright", legend = rownames(cell_iso),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 6, name = "Pastel2"), title = "Cell Type", xpd=TRUE)

n <- colSums(cell_iso)
p <- round(n/sum(n)*100,2)
text(x = x, y = n, labels = paste(n,"\n",p, sep=""), cex = .8, pos = 3)
dev.off()



#### Bar chart of cell type isoform numbers



cell_iso <- table(Cell_type$Isoforms, Cell_type$CellType)
cell_iso_plot <- apply(cell_iso,2,function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell type isoform numbers.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_iso_plot*100, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = paste(colnames(cell_iso_plot), "(", colSums(cell_iso), ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_iso_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()

cell_iso <- table(Cell_type$Isoforms, Cell_type$CellType)[-1,]
cell_iso_plot <- apply(cell_iso,2,function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell type isoform numbers2.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_iso_plot*100, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3")[-1], 
			 ylab = "%", names.arg = paste(colnames(cell_iso_plot), "(", colSums(cell_iso), ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_iso_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3")[-1], title = "Isoforms", xpd=TRUE)
dev.off()




#### Cells has more than one isoforms



cov_merge[,row.names(Cell_type)] -> cov_tu

Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
Cell_type[Cell_type$Isoforms<2,"Isoforms"] <- 0
Cell_type[Cell_type$Isoforms>0,"Isoforms"] <- 1
table(Cell_type$Isoforms)
#   0    1
# 2542  497



Cell_type[order(Cell_type$Isoforms),] -> test

ggplot(data = tsne[row.names(test),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = test$Isoforms), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())





Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
Cell_type[Cell_type$Isoforms!=4,"Isoforms"] <- 0
table(Cell_type$Isoforms)
#   0    1
# 2542  497



pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Cells have more than one isoforms.pdf")
Cell_type[order(Cell_type$Isoforms),] -> test
ggplot(data = tsne[row.names(test),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = test$Isoforms), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()






#### 



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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cd4 sc and bulk.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(tcl_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,10)))
legend("topright", inset=c(-.2,0), legend = rownames(tcl_merge),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()



mon_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Monocytes"]]), mon_merge)
mon_merge <- apply(mon_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of Monocytes sc and bulk.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(mon_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,10)))
legend("topright", inset=c(-.2,0), legend = rownames(mon_merge),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()






#### Box plot of Exon expression of different CD45 isoforms 




cov_merge[,row.names(tsne)] -> cov_tu


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


sum(colnames(cov_tu) %in% row.names(median_merge))
# [1] 3039

median_merge[colnames(cov_tu)[colSums(cov_tu>0)==0],]

I0 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==0],])
I1 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==1],])
I2 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==2],])
I3 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==3],])
I4 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==4],])
I5 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==5],])
I6 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==6],])
I7 <- rowSums(median_merge[colnames(cov_tu)[colSums(cov_tu>0)==7],])

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Box plot of Exon expression of different CD45 isoforms.pdf")
boxplot(log10(I0+1),log10(I1+1),log10(I2+1),log10(I3+1),log10(I4+1),log10(I5+1),log10(I6+1),log10(I7+1),
	names=0:7, xlab = "No. Isoforms", ylab = "log10(exon expression)")
dev.off()




#### bar chart of differents cell types CD45 transcript



cov_merge[,row.names(tsne)] -> cov_tu

cov_tu_1 <- cov_tu[,colSums(cov_tu>0)==1]
colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]
sum(colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"])
#[1] 211
cov_tu_1_Bcel <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]]
rowSums(cov_tu_1_Bcel>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   3    1  179    0    7   17    3    1


cov_tu_1_cell <- list()
for(i in 1:length(unique(Cell_type$CellType))){
	cov_tu_1_cell[[i]] <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType==as.character(unique(Cell_type$CellType))[i]]]
}

cell_isoform_1 <- t(data.frame(do.call(rbind,lapply(cov_tu_1_cell, function(x) rowSums(x>0))), row.names = as.character(unique(Cell_type$CellType))))
cell_s <- colSums(cell_isoform_1)

cell_isoform_1_plot <- as.matrix(apply(cell_isoform_1, 2, function(x) x/sum(x))*100)

barplot(cell_isoform_1_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_1_plot.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_1_plot, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = paste(colnames(cell_isoform_1_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_isoform_1_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()




cov_tu_2 <- cov_tu[,colSums(cov_tu>0)==2]
apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))
table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")))
# RABC|RB RABC|RBC  RABC|RO RAB|RABC   RAB|RB  RAB|RBC  RAC|RBC   RBC|RO
#       1       86        1       62       31        1        1        2
#  RB|RBC    RB|RO
#     129        1
# RABC|RAC RABC|RBC  RABC|RC  RABC|RO RAB|RABC   RAB|RB   RAB|RO   RA|RAB
#        4      132        9        3       89       49        2        3
#  RA|RABC    RA|RB   RA|RBC   RBC|RC   RBC|RO   RB|RBC    RB|RC    RB|RO
#       37        4        2       17        6      167        4       19
#    RC|RO
#        1

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Distribution of isoform_2.pdf")
barplot(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))), las=2)
dev.off()





cov_tu_2_top1 <- cov_tu_2[, apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")) == "RB|RBC"]

cov_tu_2_list <- list()
for(i in 1:length(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))))){
	cov_tu_2_list[[i]] <- subset.data.frame(cov_tu_2, select = apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")) == names(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))))[i])
}
names(cov_tu_2_list) <- names(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))))

lapply(cov_tu_2_list, function(x) table(Cell_type[names(x),"CellType"]))

cell_isoform_2 <- t(do.call(rbind, lapply(cov_tu_2_list, function(x) table(Cell_type[names(x),"CellType"]))))

cell_s <- colSums(cell_isoform_2)

cell_isoform_2_plot <- as.matrix(apply(cell_isoform_2, 2, function(x) x/sum(x))*100)

barplot(cell_isoform_2_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_2_plot.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_2_plot, las = 2, col = RColorBrewer::brewer.pal(n = 6, name = "Set3"), 
			 ylab = "%", names.arg = paste(colnames(cell_isoform_2_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.35,0), legend = row.names(cell_isoform_2_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 6, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()


cell_isoform_2 <- do.call(rbind, lapply(cov_tu_2_list, function(x) table(Cell_type[names(x),"CellType"])))
cell_s <- colSums(cell_isoform_2)
cell_isoform_2_plot <- as.matrix(apply(cell_isoform_2, 2, function(x) x/sum(x))*100)
barplot(cell_isoform_2_plot)

library(RColorBrewer)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_2_plot2.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_2_plot, las = 2, col = colorRampPalette(brewer.pal(n = 10, name = "Set3"))(17), 
			 ylab = "%", names.arg = paste(colnames(cell_isoform_2_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.35,0), legend = row.names(cell_isoform_2_plot),
	bty = "n", fill = colorRampPalette(brewer.pal(n = 10, name = "Set3"))(17), title = "Isoforms", xpd=TRUE)
dev.off()








#### bar chart of differents cell types CD45 major transcript



cov_merge[,row.names(tsne)] -> cov_tu

for(i in 1:ncol(cov_tu)){
  cov_tu[,i][-which.max(cov_tu[,i])] <- 0
}

cov_tu_1 <- cov_tu[,colSums(cov_tu>0)==1]
colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]
sum(colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"])
#[1] 211
cov_tu_1_Bcel <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]]
rowSums(cov_tu_1_Bcel>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   3    1  179    0    7   17    3    1


cov_tu_1_cell <- list()
for(i in 1:length(unique(Cell_type$CellType))){
  cov_tu_1_cell[[i]] <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType==as.character(unique(Cell_type$CellType))[i]]]
}

cell_isoform_1 <- t(data.frame(do.call(rbind,lapply(cov_tu_1_cell, function(x) rowSums(x>0))), row.names = as.character(unique(Cell_type$CellType))))
cell_s <- colSums(cell_isoform_1)

cell_isoform_1_plot <- as.matrix(apply(cell_isoform_1, 2, function(x) x/sum(x))*100)

barplot(cell_isoform_1_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_major_isoform_plot.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_1_plot, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
       ylab = "%", names.arg = paste(colnames(cell_isoform_1_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_isoform_1_plot),
  bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()




cov_merge[,row.names(Cell_type)] -> cov_tu
cov_tu_p <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu_p[is.na(cov_tu_p)] <- 0
major_isoform_p <- apply(cov_tu_p, 2, max)
identical(row.names(Cell_type), names(major_isoform_p))
#[1] TRUE
Cell_type$P_major <- major_isoform_p


ggplot(data = Cell_type, aes(x = major_isoform_p, colour = CellType))+
    geom_density()

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Density plot of major isoform percentage.pdf")
ggplot(data = Cell_type, aes(x = CellType, y = major_isoform_p, fill = CellType))+
    geom_violin(adjust = .2, show.legend = FALSE)+
    ylab("Percentage of major isoform")
dev.off()






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

cov_merge -> test

test[test<50] <- 0
cs <- colSums(test>0)==1

test <- test[, cs]
colnames(test[,which(test["RC",]>0)])

sashimi(junction = "1:198692374-198703297", cell = colnames(test[,which(test["RC",]>0)]), min = 1)


CD45_V2 <- read.table("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2/CD45_RABC.txt", sep = "\t")

colnames(CD45_V2)[which(CD45_V2["RAC",]>0)]
sashimi(junction = "1:198692374-198703297", cell = "UCB3.00517")
Cell_type[colnames(CD45_V2)[which(CD45_V2["RAC",]>0)],]
#                CellType Individual
# UCB3.00517 CD4+ T Cells       UCB3

colnames(CD45_V2)[which(CD45_V2["RC",]>0)]
# [1] "UCB3.01532" "UCB3.01579" "UCB4.00338"
CD45_V2[,which(CD45_V2["RC",]>0)]
Cell_type[colnames(CD45_V2)[which(CD45_V2["RC",]>0)],]
#                CellType Individual
# UCB3.01532 CD4+ T Cells       UCB3
# NA                 <NA>       <NA>
# UCB4.00338 CD8+ T Cells       UCB4

sashimi(junction = "1:198692374-198703297", cell = colnames(CD45_V2)[which(CD45_V2["RC",]>0)])


colnames(CD45_SJ)[which((CD45_SJ["SJ34",] > 0 & CD45_SJ["SJ46",] > 0 & CD45_SJ["SJ67",] > 0))]
# [1] "UCB3.00427" "UCB3.00486" "UCB3.00517"
sashimi(junction = "1:198692374-198703297", cell = colnames(CD45_SJ)[which((CD45_SJ["SJ34",] > 0 & CD45_SJ["SJ46",] > 0 & CD45_SJ["SJ67",] > 0))])





Basic_stat <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
isos <- as.data.frame(colSums(cov_tu>0))
colnames(isos) <- "Isoforms"
Basic_stat <- merge(x = Basic_stat, y = isos, by = 0, all.y = T)
Basic_stat$Isoforms <- as.factor(Basic_stat$Isoforms)


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Sequencing depth of isoform numbers.pdf")
ggplot(Basic_stat, aes(x = Isoforms, y = log10(Mapped_reads)))+
    geom_boxplot()

ggplot(Basic_stat, aes(x = Isoforms, y = log10(Input_reads)))+
    geom_boxplot()

dev.off()


ggplot(Basic_stat, aes(x = Isoforms, y = ReadsCounts))+
    geom_boxplot()

test <- Basic_stat
test$Isoforms <- as.numeric(test$Isoforms)
test[test$Isoforms == 1, "Isoforms"] <- 0
test[test$Isoforms != 0, "Isoforms"] <- 1
test$Isoforms <- factor(test$Isoforms, labels = c("0", ">0"))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Sequencing depth of isoform numbers.pdf")
ggplot(test, aes(x = Isoforms, y = log10(Mapped_reads)))+
    geom_boxplot()+
    ggtitle(paste("P.value = ",signif(t.test(test$Mapped_reads~test$Isoforms)$p.value,3), sep = ""))
dev.off()
signif(t.test(test$Mapped_reads~test$Isoforms)$p.value,3)




#### PCA of CD45



cov_merge[,row.names(Cell_type)] -> cov_tu

#cov_tu <- apply(cov_tu, 2, function(x) x/sum(x))
#cov_tu[is.na(cov_tu)] <- 0
#cov_tu[, row.names(Cell_type)] -> cov_tu
cov_tu <- cov_tu[, colSums(cov_tu)>0]

pca <- FactoMineR::PCA(t(cov_tu),ncp = 100,graph = F)
#pdf("/mnt/nfs_nas/Lab/People/tangchao/Project/Science/result/QC/PCs_percentage_off_cutoff_HF_IR.pdf")
barplot(pca$eig[,2][1:10],xlab = "PC1 to PC10",ylab = "Percentage of variance(%)")
dev.off()


pca_result <- data.frame(pca$svd$U, cell = colnames(cov_tu), stringsAsFactors=F)
pca_result <- merge(x=pca_result, y=Cell_type, by.x="cell", by.y=0)
library(ggplot2)
ggplot(pca_result,aes(x=X1,y=X2,col=CellType))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca$eig[,2][2],2),"%)",sep = ""))+
  geom_text(label="",colour="black",size=3)



#### Combat
library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(cov_tu),]
#Define center as batch
#mod <- model.matrix(~as.factor(sample.info$Cell_type), data=exprs)
pcs <- pca$svd$U[,1:8]
row.names(pcs) <- colnames(cov_tu)
colnames(pcs) <- paste("PC",1:8,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = NULL, mean.only = T)


library(Rtsne)
pca_Rtsne <- Rtsne(unique(t(combat1)), perplexity=300, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
row.names(tsne) <- row.names(unique(t(combat)))

tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/PCA_Rtsne_of_CD45.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()




combat <- sva::ComBat(dat=cov_tu, batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=cov_tu, batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL, mean.only = T)
pca_Rtsne <- Rtsne(unique(t(combat1)), perplexity=300, theta = 1)

tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
row.names(tsne) <- row.names(unique(t(combat)))

tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Rtsne_of_CD45.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()







multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(paste(filenames,"quant.sf",sep = "/"), function(x){
    tmp<-fread(x, header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "TPM",strsplit(x,"/")[[1]][9])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/xiaqiuqi/result/3039_cd45_result")

system.time(tpm_merge<-multmerge(path))

tpm_merge <- as.data.frame(tpm_merge)
row.names(tpm_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")
colnames(tpm_merge) <- substr(colnames(tpm_merge), nchar(colnames(tpm_merge))-9, nchar(colnames(tpm_merge)))



tpm_merge[,row.names(Cell_type)] -> tpm_tu

#tpm_tu <- apply(tpm_tu, 2, function(x) x/sum(x))
#tpm_tu[is.na(tpm_tu)] <- 0
#tpm_tu[, row.names(Cell_type)] -> tpm_tu
tpm_tu <- tpm_tu[, colSums(tpm_tu)>0]

pca <- FactoMineR::PCA(t(tpm_tu),ncp = 100,graph = F)
#pdf("/mnt/nfs_nas/Lab/People/tangchao/Project/Science/result/QC/PCs_percentage_off_cutoff_HF_IR.pdf")
barplot(pca$eig[,2][1:10],xlab = "PC1 to PC10",ylab = "Percentage of variance(%)")
dev.off()


pca_result <- data.frame(pca$svd$U, cell = colnames(tpm_tu), stringsAsFactors=F)
pca_result <- merge(x=pca_result, y=Cell_type, by.x="cell", by.y=0)
library(ggplot2)
ggplot(pca_result,aes(x=X1,y=X2,col=CellType))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca$eig[,2][2],2),"%)",sep = ""))+
  geom_text(label="",colour="black",size=3)



#### Combat
library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(tpm_tu),]
#Define center as batch
#mod <- model.matrix(~as.factor(sample.info$Cell_type), data=exprs)
pcs <- pca$svd$U[,1:8]
row.names(pcs) <- colnames(tpm_tu)
colnames(pcs) <- paste("PC",1:8,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = NULL, mean.only = T)


library(Rtsne)
pca_Rtsne <- Rtsne(unique(t(combat1)), perplexity=300, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
row.names(tsne) <- row.names(unique(t(combat)))

tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/PCA_Rtsne_of_CD45.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



combat <- sva::ComBat(dat=tpm_tu, batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=tpm_tu, batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL, mean.only = T)
pca_Rtsne <- Rtsne(unique(t(combat)), perplexity=300, theta = 1)

tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
row.names(tsne) <- row.names(unique(t(combat)))

tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Rtsne_of_CD45_use_TPM.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



