#### CD45 Statistic
library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(paste(filenames,"t_data.ctab",sep = "/"), function(x){
    tmp<-fread(x, header=TRUE, select = 11, sep = "\t")
    setnames(tmp, "cov",strsplit(x,"/")[[1]][10])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result")

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_merge[cov_merge<10] <- 0

cov_merge[,row.names(tsne)] -> cov_tu

cov_tu[cov_tu<50] <- 0


Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

apply(cov_tu[,row.names(Cell_type)], 2, function(x))

Cell_type$Isoforms <- colSums(cov_tu[,row.names(Cell_type)]>0)
table(Cell_type$CellType, Cell_type$Isoforms)
#                   0   1   2   3   4
#  B Cells        317 249  30   6   5
#  CD4+ T Cells   762 592 209  58  57
#  CD8+ T Cells   257 213  62  29  18
#  Megakaryocytes  10   9   3   1   1
#  Monocytes       34  13   3   1   1
#  NK cells        41  45   8   1   4

cell_iso <- table(Cell_type$CellType, Cell_type$Isoforms)


#### Bar chart of cell isoform numbers


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cell isoform numbers.pdf")
par(mar=c(10.1, 4.1, 4.1, 2.1))
x <- barplot(cell_iso, col = RColorBrewer::brewer.pal(n = 6, name = "Pastel2"), ylim = c(0,1600),
			 ylab = "No. Cells", xlab = "No. Isoforms", ames.arg = colnames(cell_iso))
legend("topright", legend = rownames(cell_iso),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 6, name = "Pastel2"), title = "Cell Type", xpd=TRUE)

n <- colSums(cell_iso)
p <- round(n/sum(n)*100,2)
text(x = x, y = n, labels = paste(n,"\n",p, sep=""), cex = .8, pos = 3)
dev.off()



#### Bar chart of cell type isoform numbers



cell_iso <- table(Cell_type$Isoforms, Cell_type$CellType)
cell_iso_plot <- apply(cell_iso,2,function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cell type isoform numbers.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_iso_plot*100, las = 2, col = RColorBrewer::brewer.pal(n = 5, name = "Set3"), 
			 ylab = "%", names.arg = paste(colnames(cell_iso_plot), "(", colSums(cell_iso), ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_iso_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 5, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()

cell_iso <- table(Cell_type$Isoforms, Cell_type$CellType)[-1,]
cell_iso_plot <- apply(cell_iso,2,function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cell type isoform numbers2.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_iso_plot*100, las = 2, col = RColorBrewer::brewer.pal(n = 5, name = "Set3")[-1], 
			 ylab = "%", names.arg = paste(colnames(cell_iso_plot), "(", colSums(cell_iso), ")", sep = " "))
legend("topright", inset=c(-.2,0), legend = rownames(cell_iso_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 5, name = "Set3")[-1], title = "Isoforms", xpd=TRUE)
dev.off()








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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/CD45_RABC_v5.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()




#### Only plot the major isoform


cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[cov_tu<50] <- 0

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


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/CD45_RABC_v5_2.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()



#### Only plot the cells have only one isoform



cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[cov_tu<50] <- 0

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


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/CD45_RABC_v5_3.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()





#### Cells has tow isoforms



cov_merge[,row.names(Cell_type)] -> cov_tu
cov_tu[cov_tu<50] <- 0

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


























path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result_bulk/Tcell")
system.time(tcl_merge<-multmerge(path))

tcl_merge <- as.data.frame(tcl_merge)
row.names(tcl_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")
tcl_merge[tcl_merge<50] <- 0

#  tcl_merge
#          Tcell    Tcell    Tcell     Tcell     Tcell     Tcell      Tcell
#  RA     0.0000   0.0000   0.0000   0.00000   0.00000   0.00000    0.00000
#  RAB  291.7117 270.4915 301.7529 348.56516 354.85364 275.19565  588.04749
#  RABC 448.5605 954.9783 625.5502 754.70734 842.12384 685.38928 1337.64563
#  RAC    0.0000   0.0000   0.0000   0.00000   0.00000   0.00000    0.00000
#  RB     0.0000   0.0000   0.0000  50.50396  75.73828  61.29072   68.16885
#  RBC  166.4196 200.4341 201.0803 248.20618 325.43601 250.35371  326.51446
#  RC     0.0000   0.0000   0.0000   0.00000   0.00000   0.00000    0.00000
#  RO     0.0000   0.0000   0.0000   0.00000   0.00000   0.00000    0.00000
#           Tcell    Tcell    Tcell     Tcell
#  RA      0.0000   0.0000   0.0000   0.00000
#  RAB   513.1597 201.3109 219.8292 307.14709
#  RABC 1407.3457 520.2159 494.6178 761.75238
#  RAC     0.0000   0.0000   0.0000   0.00000
#  RB      0.0000   0.0000   0.0000  50.03038
#  RBC   109.7562 109.7448 188.8185 213.48045
#  RC      0.0000   0.0000   0.0000   0.00000
#  RO      0.0000   0.0000   0.0000   0.00000
#  rowSums(tcl_merge>0)
#    RA  RAB RABC  RAC   RB  RBC   RC   RO
#     0   11   11    0    5   11    0    0



path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result_bulk/Monocytes")
system.time(mon_merge<-multmerge(path))

mon_merge <- as.data.frame(mon_merge)
row.names(mon_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")
mon_merge[mon_merge<50] <- 0

#       Monocytes Monocytes Monocytes Monocytes Monocytes Monocytes Monocytes
#  RA     0.00000    0.0000   0.00000   0.00000   0.00000   0.00000   0.00000
#  RAB  203.72086  111.6921 149.03250 163.34731 124.15694 122.23420 131.10495
#  RABC 200.28134  145.8406 142.39612 189.57208 141.84233 143.43378 144.81928
#  RAC    0.00000    0.0000   0.00000   0.00000   0.00000   0.00000   0.00000
#  RB   612.22656 1054.1984 469.14236 707.59515 646.01996 900.87695 725.86114
#  RBC   92.64352  122.7750  63.72162  80.13504  72.82333  67.39078  96.08542
#  RC     0.00000    0.0000   0.00000   0.00000   0.00000   0.00000   0.00000
#  RO     0.00000  103.0294   0.00000   0.00000  83.45482  91.92748  58.14054
#       Monocytes Monocytes Monocytes
#  RA     0.00000   0.00000   0.00000
#  RAB  142.34024 117.39995 100.06690
#  RABC 125.78239 144.21730 143.69282
#  RAC    0.00000   0.00000   0.00000
#  RB   717.62921 591.85785 661.19898
#  RBC   66.54792  94.47496  84.79369
#  RC     0.00000   0.00000   0.00000
#  RO    53.38268  55.36517   0.00000



tcl_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="CD4+ T Cells"]]), tcl_merge)
tcl_merge <- apply(tcl_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cd4 sc and bulk.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(tcl_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,11)))
legend("topright", inset=c(-.2,0), legend = rownames(tcl_merge),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()



mon_merge <- cbind(rowSums(cov_tu[,row.names(Cell_type)[Cell_type$CellType=="Monocytes"]]), mon_merge)
mon_merge <- apply(mon_merge, 2, function(x) x/sum(x))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of Monocytes sc and bulk.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(mon_merge, las = 2, col = RColorBrewer::brewer.pal(n = 8, name = "Set3"), 
			 ylab = "%", names.arg = rep(c("SC","Bulk"), c(1,11)))
legend("topright", inset=c(-.2,0), legend = rownames(mon_merge),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 8, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()









#### Box plot of Exon expression of different CD45 isoforms 

cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[cov_tu<50] <- 0


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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Box plot of Exon expression of different CD45 isoforms.pdf")
boxplot(log10(I0+1),log10(I1+1),log10(I2+1),log10(I3+1),log10(I4+1), names=0:4, xlab = "No. Isoforms", ylab = "log10(exon expression)")
dev.off()




#### bar chart of differents cell types CD45 transcript

cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[cov_tu<50] <- 0

cov_tu_1 <- cov_tu[,colSums(cov_tu>0)==1]
colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]
sum(colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"])
#[1] 249
cov_tu_1_Bcel <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType=="B Cells"]]
rowSums(cov_tu_1_Bcel>0)
#  RA  RAB RABC  RAC   RB  RBC   RC   RO
#   0    5  215    0    6   23    0    0

cov_tu_1_cell <- list()
for(i in 1:length(unique(Cell_type$CellType))){
	cov_tu_1_cell[[i]] <- cov_tu_1[,colnames(cov_tu_1) %in% row.names(Cell_type)[Cell_type$CellType==as.character(unique(Cell_type$CellType))[i]]]
}

cell_isoform_1 <- t(data.frame(do.call(rbind,lapply(cov_tu_1_cell, function(x) rowSums(x>0))), row.names = as.character(unique(Cell_type$CellType))))
cell_s <- colSums(cell_isoform_1)

cell_isoform_1_plot <- as.matrix(apply(cell_isoform_1, 2, function(x) x/sum(x))*100)

barplot(cell_isoform_1_plot)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cell_isoform_1_plot.pdf")
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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Distribution of isoform_2.pdf")
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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cell_isoform_2_plot.pdf")
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

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat/Bar chart of cell_isoform_2_plot2.pdf")
par(mar=c(10.1, 4.1, 4.1, 8.1))
x <- barplot(cell_isoform_2_plot, las = 2, col = RColorBrewer::brewer.pal(n = 10, name = "Set3"), 
			 ylab = "%", names.arg = paste(colnames(cell_isoform_2_plot), "(", cell_s, ")", sep = " "))
legend("topright", inset=c(-.35,0), legend = row.names(cell_isoform_2_plot),
	bty = "n", fill = RColorBrewer::brewer.pal(n = 10, name = "Set3"), title = "Isoforms", xpd=TRUE)
dev.off()





#### PCA of CD45

cov_merge[,row.names(tsne)] -> cov_tu
cov_tu[cov_tu<50] <- 0

cov_tu <- apply(cov_tu, 2, function(x) x/sum(x))
cov_tu[is.na(cov_tu)] <- 0


pca <- FactoMineR::PCA(t(cov_tu),ncp = 10,graph = F)
pdf("/mnt/nfs_nas/Lab/People/tangchao/Project/Science/result/QC/PCs_percentage_off_cutoff_HF_IR.pdf")
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

library(ica)

imod <- icafast(t(cov_tu),2)

imod$Y
pca_result <- data.frame(imod$Y, cell = colnames(cov_tu), stringsAsFactors=F)
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








