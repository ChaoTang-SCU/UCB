
setwd("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4")
load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te_no_NA_and_normalized.RData")
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

as.data.frame(te[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),]) -> CD45_SJ
CD45_SJ[CD45_SJ<10] <- 0
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

sum(colSums(CD45_SJ[1:9,])>=10)
# [1] 2537
## sum of junction reads (except the cells only have SJ37) must more than 10.
#test <- as.data.frame(apply(CD45_SJ[1:9,],2, function(x) x/sum(x) <= 0.05))
#test[is.na(test)] <- TRUE
#test["SJ37",] <- FALSE
#CD45_SJ[as.matrix(test)] <- 0



CD45_R_r_tab <- CD45_SJ[1:8,]
row.names(CD45_R_r_tab) <- c("RO","RAC","RBC","RAB","RA","RC","RB","RABC")
CD45_R_r_tab[CD45_R_r_tab>0] <- 0

for (i in 1:ncol(CD45_SJ)) {
  if(CD45_SJ["SJ37",i] > 0){
    CD45_R_r_tab["RO",i] <- CD45_SJ["SJ37",i]
  }
  if(CD45_SJ["SJ46",i] > 0){
    CD45_R_r_tab["RAC",i] <- CD45_SJ["SJ46",i]
  }
  if(CD45_SJ["SJ35",i] > 0 & CD45_SJ["SJ56",i] > 0 & CD45_SJ["SJ67",i] > 0){
  	if(CD45_SJ["SJ34",i] > 0 & CD45_SJ["SJ45",i] > 0){
  		CD45_R_r_tab["RBC",i] <- CD45_SJ["SJ35",i]
  	}else{
  		CD45_R_r_tab["RBC",i] <- mean(c(CD45_SJ["SJ56",i], CD45_SJ["SJ67",i]))
  	}#yin wei dang RBC yu RB tong shi cun zai de shi hou, yong J35 ji suan PSI RBC hui bu dui.
  }
  if(CD45_SJ["SJ45",i] > 0 & CD45_SJ["SJ57",i] > 0 & CD45_SJ["SJ34",i] > 0){
  	if(CD45_SJ["SJ56",i] > 0 & CD45_SJ["SJ67",i] > 0){
  		CD45_R_r_tab["RAB",i] <- CD45_SJ["SJ57",i]
  	}else{
  		CD45_R_r_tab["RAB",i] <- mean(c(CD45_SJ["SJ34",i], CD45_SJ["SJ45",i]))
  	}    
  }
  if(CD45_SJ["SJ47",i] > 0 & CD45_SJ["SJ34",i] > 0){
    CD45_R_r_tab["RA",i] <- CD45_SJ["SJ47",i]
  }
  if(CD45_SJ["SJ36",i] > 0 & CD45_SJ["SJ67",i] > 0){
    CD45_R_r_tab["RC",i] <- CD45_SJ["SJ36",i]
  }
  if(CD45_SJ["SJ35",i] > 0 & CD45_SJ["SJ57",i] > 0){
  	if(sum(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),i])>0){
  		CD45_R_r_tab["RB",i] <- min(c(CD45_SJ["SJ35",i], CD45_SJ["SJ57",i]))
  		}else{
  			CD45_R_r_tab["RB",i] <- mean(c(CD45_SJ["SJ35",i], CD45_SJ["SJ57",i]))
  		}
  }# dang RB yu RAB huozhe RBC tongshi cunzai de shihou, xuyao jianqu RAB huozhe RBC de reads.
  # zheli women congshi zhizhong dou meiyou kaolv tongshi cunzai 4 zhong isoform de qingkuang, shi yinwei women houmian yao jiang 4 zhong isoform de qingkuang tiaochulai chongxin dingyi.
  if(CD45_SJ["SJ34",i] > 0 & CD45_SJ["SJ45",i] > 0 & CD45_SJ["SJ56",i] > 0 & CD45_SJ["SJ67",i] > 0){
  	if(CD45_SJ["SJ35",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ34",i], CD45_SJ["SJ45",i]))
}
if(CD45_SJ["SJ36",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ34",i], CD45_SJ["SJ45",i],CD45_SJ["SJ56",i]))
}
if(CD45_SJ["SJ46",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ45",i], CD45_SJ["SJ56",i]))
}
if(CD45_SJ["SJ47",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ45",i], CD45_SJ["SJ56",i], CD45_SJ["SJ67",i]))
}
if(CD45_SJ["SJ57",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ56",i], CD45_SJ["SJ67",i]))
}
if(CD45_SJ["SJ35",i] == 0 & CD45_SJ["SJ36",i] == 0 & CD45_SJ["SJ46",i] == 0 & CD45_SJ["SJ47",i] == 0 & CD45_SJ["SJ57",i] == 0){
    CD45_R_r_tab["RABC",i] <- (CD45_SJ["SJ34",i] + CD45_SJ["SJ45",i] + CD45_SJ["SJ56",i] + CD45_SJ["SJ67",i])/4
}
}
}





CD45_R_r_tab[,1:5]
table(colSums(CD45_R_r_tab>0))
   0    1    2    3    4    5    6
1449 1282  529   29  273   11    1

table(colSums(CD45_R_r_tab!=0))
   0    1    2    3    4    5    6
1449 1282  529   29  273   11    1


CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2][,1:5]
CD45_SJ[,colSums(CD45_R_r_tab!=0)==2][,1:5]

CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3][,1:5]
CD45_SJ[,colSums(CD45_R_r_tab!=0)==3][,1:5]

CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==4][,1:5]
CD45_SJ[,colSums(CD45_R_r_tab!=0)==4][,1:5]

recalculete <- colnames(CD45_R_r_tab)[c(colSums(CD45_R_r_tab!=0)==4 | colSums(CD45_R_r_tab!=0)==5 | colSums(CD45_R_r_tab!=0)==6)]
dim(CD45_R_r_tab[, recalculete])
# [1]   8 298
dim(CD45_SJ[, recalculete])
# [1]  10 298
rowSums(CD45_SJ[, recalculete])
#  SJ34  SJ35  SJ45  SJ46  SJ56  SJ57  SJ67  SJ47  SJ36  SJ37
# 80656 90594 72565     0 73570 78443 75796     0     0  1524
# 148669 218469 132351      0 183265 148432 192620      0    458   3511



for(i in 1:length(recalculete)){
	# -+-, -+-
	if(CD45_SJ["SJ35",recalculete[i]] > mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] > mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- mean(c(CD45_SJ["SJ35",recalculete[i]], CD45_SJ["SJ57",recalculete[i]]))
		CD45_R_r_tab["RAB",recalculete[i]] <- 0
		CD45_R_r_tab["RBC",recalculete[i]] <- 0
		CD45_R_r_tab["RABC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/4
	}
	# +-+, -+-
	if(CD45_SJ["SJ35",recalculete[i]] < mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] > mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- 0
		CD45_R_r_tab["RAB",recalculete[i]] <- (CD45_SJ["SJ35",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/3
		CD45_R_r_tab["RBC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ57",recalculete[i]])/3
		CD45_R_r_tab["RABC",recalculete[i]] <- 0
	}
	# -+-, +-+
	if(CD45_SJ["SJ35",recalculete[i]] > mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] < mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- 0
		CD45_R_r_tab["RAB",recalculete[i]] <- (CD45_SJ["SJ35",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/3
		CD45_R_r_tab["RBC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ57",recalculete[i]])/3
		CD45_R_r_tab["RABC",recalculete[i]] <- 0
	}
	# +-+, +-+
	if(CD45_SJ["SJ35",recalculete[i]] < mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] < mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- mean(c(CD45_SJ["SJ35",recalculete[i]], CD45_SJ["SJ57",recalculete[i]]))
		CD45_R_r_tab["RAB",recalculete[i]] <- 0
		CD45_R_r_tab["RBC",recalculete[i]] <- 0
		CD45_R_r_tab["RABC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/4
	}
}


table(colSums(CD45_R_r_tab!=0))
    0    1    2    3    4
 1416 1288  826   43    1
rowSums(CD45_R_r_tab>0)
  RO  RAC  RBC  RAB   RA   RC   RB RABC
 162    1  890  343    3    5  808  861


write.table(CD45_R_r_tab,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC.txt", sep = "\t")
write.table(CD45_SJ,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_SJ.txt", sep = "\t")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

library(ggplot2)
t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab2)[1])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC_v4.pdf", height = 5,width = 10)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, ncol = 4)
dev.off()



dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==1])

CD45_R_r_tab -> chao
chao[,colSums(CD45_R_r_tab>0)>1] <- 0
dim(chao[,colSums(chao>0)==1])
dim(chao[,colSums(chao>0)>1])
#[1] 8 0
chao <- chao[rowSums(chao)>0,]
dim(chao)

library(ggplot2)
t(t(chao)[order(t(chao)[,1]),]) -> chao1
f1 <- ggplot(data = tsne[colnames(chao1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao1[1,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(chao2)[1])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,2]),]) -> chao2
f2 <- ggplot(data = tsne[colnames(chao2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao2[2,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(chao2)[2])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,3]),]) -> chao3
f3 <- ggplot(data = tsne[colnames(chao3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao3[3,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(chao3)[3])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,4]),]) -> chao4
f4 <- ggplot(data = tsne[colnames(chao4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = round(log10(as.numeric(chao4[4,]+1)))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(chao4)[4])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,5]),]) -> chao5
f5 <- ggplot(data = tsne[colnames(chao5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao5[5,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(chao5)[5])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,6]),]) -> chao6
f6 <- ggplot(data = tsne[colnames(chao6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao6[6,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(chao6)[6])+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC_v4_only_one_type.pdf", height = 4,width = 8)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, ncol = 3)
dev.off()



#### Heat map
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)




library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC_v4_heatmap.pdf", width = 7, height = 7)
pheatmap(t(log10(CD45_R_r_tab[,row.names(Cell_type)]+1)), annotation_row = Cell_type, show_rownames = F)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC_v4_heatmap.2.pdf", width = 7, height = 7)
pheatmap(t(log10(CD45_R_r_tab[,row.names(Cell_type)]+1)), annotation_row = Cell_type, show_rownames = F, cluster_rows = F)
dev.off()

library(RColorBrewer)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC_v4_heatmap.3.pdf", width = 7, height = 7)
pheatmap(t(log10(CD45_R_r_tab[c(8,7,3,4,1),row.names(Cell_type)]+1)), 
  annotation_row = Cell_type, show_rownames = F, cluster_rows = F, cluster_cols = F,
  color = colorRampPalette(brewer.pal(n = 7, name = "Reds"), bias = .5)(100))
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_RABC_v4_heatmap.4.pdf", width = 7, height = 7)
pheatmap(t(log10(CD45_R_r_tab[c(8,7,3,4,1),row.names(Cell_type)]+1)), 
  annotation_row = Cell_type, show_rownames = F, cluster_rows = F, cluster_cols = F)
dev.off()







#### Exon expression assignment: ================================================================================================================

table(apply(CD45_R_r_tab[, colSums(CD45_R_r_tab)>0],2,function(x) paste(row.names(CD45_R_r_tab)[x>0], collapse = "|")))
           RA           RAB          RABC        RAB|RA      RAB|RABC
            2            84           474             1            38
       RAB|RB       RAC|RBC            RB           RBC       RBC|RAB
           66             1           342           315           141
     RBC|RABC        RBC|RB       RB|RABC       RC|RABC    RC|RB|RABC
          171           196           132             1             1
           RO        RO|RAB       RO|RABC     RO|RAB|RB         RO|RB
           65             4             6             3            32
       RO|RBC    RO|RBC|RAB   RO|RBC|RABC RO|RBC|RAB|RC     RO|RBC|RB
           13             4             7             1            17
   RO|RB|RABC      RO|RC|RB
            6             2

table(apply(CD45_R_r_tab[, colSums(CD45_R_r_tab>0)>1],2,function(x) paste(row.names(CD45_R_r_tab)[x>0], collapse = "|")))
       RAB|RA      RAB|RABC        RAB|RB       RAC|RBC       RBC|RAB
            1            38            66             1           141
     RBC|RABC        RBC|RB       RB|RABC       RC|RABC    RC|RB|RABC
          171           196           132             1             1
       RO|RAB       RO|RABC     RO|RAB|RB         RO|RB        RO|RBC
            4             6             3            32            13
   RO|RBC|RAB   RO|RBC|RABC RO|RBC|RAB|RC     RO|RBC|RB    RO|RB|RABC
            4             7             1            17             6
     RO|RC|RB
            2




#### Intron centric PSI: ========================================================================================================================
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

as.data.frame(psi_sj_same_start_table[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),row.names(tsne)]) -> psi
row.names(psi) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Intron_centric_PSI4.pdf", height = 5,width = 10)
for (i in 1:nrow(psi)){
  tmp_data <- data.frame(tsne, PSI = as.numeric(psi[i,row.names(tsne)]))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(psi)[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()


for (i in 1:nrow(psi)){
  tmp_data <- data.frame(tsne, PSI = as.numeric(psi[i,row.names(tsne)]))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  assign(paste("p",i,sep = ""),ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
           geom_point(aes(colour = PSI), size = .6)+
           scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
           ggtitle(row.names(psi)[i])+
           theme(panel.background = element_blank(),
                 legend.position = "none",
                 axis.title = element_blank(),
                 axis.text= element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 panel.border = element_blank()))
}

p10 <- ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = .6)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90")+
      ggtitle(row.names(psi)[i])+
      theme(panel.background = element_blank(),
            legend.title = element_blank(),
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Intron_centric_PSI5.pdf", height = 6,width = 10)
library(cowplot)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 4)
dev.off()




for (i in 1:nrow(psi)){
  tmp_data <- data.frame(tsne, PSI = as.numeric(psi[i,row.names(tsne)]))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  assign(paste("p",i,sep = ""),ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
           geom_point(aes(colour = PSI), size = .6)+
           scale_colour_gradient2(low = "grey80", high = "#EF3B2C", na.value = "grey90", midpoint = 0.8, guide = guide_legend(title = "PSI"))+
           ggtitle(row.names(psi)[i])+
           theme(panel.background = element_blank(),
                 legend.position = "none",
                 axis.title = element_blank(),
                 axis.text= element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 panel.border = element_blank()))
}

p10 <- ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = .6)+
      scale_colour_gradient2(low = "grey80", high = "#EF3B2C", na.value = "grey90", midpoint = 0.8)+
      ggtitle(row.names(psi)[i])+
      theme(panel.background = element_blank(),
            legend.title = element_blank(),
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Intron_centric_PSI6.pdf", height = 6,width = 10)
library(cowplot)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 4)
dev.off()




