#### CCD45 of salmon
library(data.table)
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

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")
colnames(cov_merge) <- substr(colnames(cov_merge), nchar(colnames(cov_merge))-9, nchar(colnames(cov_merge)))

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

#cov_merge[cov_merge<5000] <- 0

cov_merge[,row.names(tsne)] -> cov_tu

cov_tu <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu[is.na(cov_tu)] <- 0
cov_tu[cov_tu<=0.1] <- 0


library(ggplot2)
t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[1,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab2[2,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab3[3,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab4[4,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab5[5,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab6[6,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab7[7,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab8[8,])), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v6/CD45_RABC_v6_TPM.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()







#### Only plot the major isoform




#cov_merge[,row.names(tsne)] -> cov_tu

for(i in 1:ncol(cov_tu)){
  cov_tu[,i][-which.max(cov_tu[,i])] <- 0
}


t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[1,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab2[2,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab3[3,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab4[4,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab5[5,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab6[6,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab7[7,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab8[8,])), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v6/CD45_RABC_v6_2_TPM.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()






cov_merge[,row.names(tsne)] -> cov_tu

cov_tu <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu[is.na(cov_tu)] <- 0
cov_tu[cov_tu<=0.1] <- 0

table(apply(cov_tu,2,function(x) sum(x>0)))
#    0    1    2    3    4    5    6
#  386 1223  797  425  190   17    1

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v6/Bar chart of cell isoform numbers TPM.pdf")
x = barplot(table(apply(cov_tu,2,function(x) sum(x>0))), ylim = c(0,1300),
  ylab = "No. Cells", col = "white", xlab = "No. Isoforms")
n <- table(apply(cov_tu,2,function(x) sum(x>0)))
text(x = x, y = n, labels = n, cex = .8, pos = 3)
dev.off()



#### Only plot the cells have only one isoform


for(i in 1:ncol(cov_tu)){
  if(sum(cov_tu[,i] > 0 )>1){
    cov_tu[,i] <- 0
  } 
}


t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[1,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab2[2,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab3[3,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab4[4,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab5[5,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab6[6,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab7[7,])), size = 1)+
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
  geom_point(aes(colour = as.numeric(CD45_R_r_tab8[8,])), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v6/CD45_RABC_v6_3_TPM.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()







cov_merge[,row.names(tsne)] -> cov_tu

cov_tu <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu[is.na(cov_tu)] <- 0
cov_tu[cov_tu<=0.1] <- 0


Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

apply(cov_tu[,row.names(Cell_type)], 2, function(x))

Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
#                   0   1   2   3   4   5   6
#  B Cells         92 247 150  89  24   4   1
#  CD4+ T Cells   211 676 446 227 107  11   0
#  CD8+ T Cells    53 233 162  82  47   2   0
#  Megakaryocytes   4  10   4   1   5   0   0
#  Monocytes       10  14  15  12   1   0   0
#  NK cells        16  43  20  14   6   0   0




cell_iso <- table(Cell_type$CellType, Cell_type$Isoforms)


#### Bar chart of cell isoform numbers


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell isoform numbers TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 2.1))
x <- barplot(cell_iso, col = RColorBrewer::brewer.pal(n = 6, name = "Pastel2"), ylim = c(0,1300),
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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell type isoform numbers TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_iso_plot*100, las = 2, col = RColorBrewer::brewer.pal(n = 7, name = "Set3"), 
       ylab = "%", names.arg = paste(colnames(cell_iso_plot), "(", colSums(cell_iso), ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_iso_plot),
  bty = "n", fill = RColorBrewer::brewer.pal(n = 7, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()

cell_iso <- table(Cell_type$Isoforms, Cell_type$CellType)[-1,]
cell_iso_plot <- apply(cell_iso,2,function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell type isoform numbers2 TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_iso_plot*100, las = 2, col = RColorBrewer::brewer.pal(n = 7, name = "Set3")[-1], 
       ylab = "%", names.arg = paste(colnames(cell_iso_plot), "(", colSums(cell_iso), ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_iso_plot),
  bty = "n", fill = RColorBrewer::brewer.pal(n = 7, name = "Set3")[-1], title = "Isoforms", xpd=TRUE)
dev.off()




#### Cells has more than one isoforms



cov_merge[,row.names(Cell_type)] -> cov_tu

Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
Cell_type[Cell_type$Isoforms<2,"Isoforms"] <- 0
Cell_type[Cell_type$Isoforms>0,"Isoforms"] <- 1
table(Cell_type$Isoforms)
#   0    1
# 396 2643



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
#    0    1
# 2734  305



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
    tmp<-fread(x, header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "TPM",strsplit(x,"/")[[1]][9])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

system.time(bulk_merge<-multmerge(path))

bulk_merge <- as.data.frame(bulk_merge)
row.names(bulk_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

mon_merge <- bulk_merge[,grep("mono", colnames(bulk_merge))]
tcl_merge <- bulk_merge[,grep("mono", colnames(bulk_merge), invert=T)]

tcl_merge <- apply(tcl_merge, 2, function(x) x/sum(x))
tcl_merge[is.na(tcl_merge)] <- 0
tcl_merge[tcl_merge<=0.1] <- 0

rowSums(tcl_merge>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   0   10   10    0    0   10    0    0

mon_merge <- apply(mon_merge, 2, function(x) x/sum(x))
mon_merge[is.na(mon_merge)] <- 0
mon_merge[mon_merge<=0.1] <- 0
rowSums(mon_merge>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   0    7    0    0   10    0    0   10


tcl_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD4+ T Cells"]]), tcl_merge)
tcl_merge <- apply(tcl_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cd4 sc and bulk TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(tcl_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
       ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,10)))
legend("topright", inset=c(-.2,0), legend = rownames(tcl_merge),
  bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()



mon_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Monocytes"]]), mon_merge)
mon_merge <- apply(mon_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of Monocytes sc and bulk TPM.pdf")
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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Box plot of Exon expression of different CD45 isoforms TPM.pdf")
boxplot(log10(I0+1),log10(I1+1),log10(I2+1),log10(I3+1),log10(I4+1),log10(I5+1),log10(I6+1),
  names=0:6, xlab = "No. Isoforms", ylab = "log10(exon expression)")
dev.off()




#### bar chart of differents cell types CD45 transcript



cov_merge[,row.names(tsne)] -> cov_tu

cov_tu <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu[is.na(cov_tu)] <- 0
cov_tu[cov_tu<=0.1] <- 0

cov_tu_1 <- cov_tu[,colSums(cov_tu>0)==1]
colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]
sum(colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"])
#[1] 247
cov_tu_1_Bcel <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]]
rowSums(cov_tu_1_Bcel>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   2    2  205    0    9   23    4    2


cov_tu_1_cell <- list()
for(i in 1:length(unique(Cell_type$CellType))){
  cov_tu_1_cell[[i]] <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType==as.character(unique(Cell_type$CellType))[i]]]
}

cell_isoform_1 <- t(data.frame(do.call(rbind,lapply(cov_tu_1_cell, function(x) rowSums(x>0))), row.names = as.character(unique(Cell_type$CellType))))
cell_s <- colSums(cell_isoform_1)

cell_isoform_1_plot <- as.matrix(apply(cell_isoform_1, 2, function(x) x/sum(x))*100)

barplot(cell_isoform_1_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_1_plot TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_1_plot, las = 2, col = RColorBrewer::brewer.pal(n = 7, name = "Set3"), 
       ylab = "%", names.arg = paste(colnames(cell_isoform_1_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_isoform_1_plot),
  bty = "n", fill = RColorBrewer::brewer.pal(n = 7, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()




cov_tu_2 <- cov_tu[,colSums(cov_tu>0)==2]
apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))
table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Distribution of isoform_2 TPM.pdf")
barplot(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))), las=2)
dev.off()





cov_tu_2_top1 <- cov_tu_2[, apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")) == "RB|RBC"]

cov_tu_2_list <- list()
for(i in 1:length(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))))){
  cov_tu_2_list[[i]] <- subset.data.frame(cov_tu_2, select = apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|")) == names(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))))[i])
}
names(cov_tu_2_list) <- names(table(apply(cov_tu_2, 2, function(x) paste(row.names(cov_tu_2)[x>0], collapse="|"))))

lapply(cov_tu_2_list, function(x) table(Cell_type[names(x),"CellType"]))

cell_isoform_2 <- t(do.call(rbind, lapply(cov_tu_2_list, function(x) table(Cell_type[colnames(x),"CellType"]))))

cell_s <- colSums(cell_isoform_2)

cell_isoform_2_plot <- as.matrix(apply(cell_isoform_2, 2, function(x) x/sum(x))*100)

barplot(cell_isoform_2_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_2_plot TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_2_plot, las = 2, col = RColorBrewer::brewer.pal(n = 6, name = "Set3"), 
       ylab = "%", names.arg = paste(colnames(cell_isoform_2_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.35,0), legend = row.names(cell_isoform_2_plot),
  bty = "n", fill = RColorBrewer::brewer.pal(n = 6, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()


cell_isoform_2 <- do.call(rbind, lapply(cov_tu_2_list, function(x) table(Cell_type[colnames(x),"CellType"])))
cell_s <- colSums(cell_isoform_2)
cell_isoform_2_plot <- as.matrix(apply(cell_isoform_2, 2, function(x) x/sum(x))*100)
barplot(cell_isoform_2_plot)

library(RColorBrewer)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_2_plot2 TPM.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_2_plot, las = 2, col = colorRampPalette(brewer.pal(n = 10, name = "Set3"))(20), 
       ylab = "%", names.arg = paste(colnames(cell_isoform_2_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.35,0), legend = row.names(cell_isoform_2_plot),
  bty = "n", fill = colorRampPalette(brewer.pal(n = 10, name = "Set3"))(20), title = "Isoforms", xpd=TRUE)
dev.off()

cell_isoform_2_plot[cell_isoform_2_plot<10] <- 0
cell_isoform_2_plot <- cell_isoform_2_plot[rowSums(cell_isoform_2_plot)>0,]

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat2/Bar chart of cell_isoform_2_plot2 TPM 2.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_2_plot, las = 2, col = brewer.pal(n = 9, name = "Set3"), 
       ylab = "%", names.arg = paste(colnames(cell_isoform_2_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.35,0), legend = row.names(cell_isoform_2_plot),
  bty = "n", fill = brewer.pal(n = 9, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()

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


















