
#### CDC45 in bulk
load('/mnt/data1/projects/UCB/data_ucb/BP_ucb_bulk_SJ/psi_same_start_end_SJ.RData')
load("/mnt/data1/projects/UCB/results_ucb/tangchao/bulk_gene_and_psi/RData/bulk_gene_Seurat.RData") ## just for the sample information
library(data.table)
SJ34 = "chr1:198692374-198696711"
SJ35 = "chr1:198692374-198699563"
SJ45 = "chr1:198696910-198699563"
SJ46 = "chr1:198696910-198702386"
SJ56 = "chr1:198699705-198702386"
SJ57 = "chr1:198699705-198703297"
SJ67 = "chr1:198702531-198703297"
SJ47 = "chr1:198696910-198703297"
SJ36 = "chr1:198692374-198702386"
SJ37 = "chr1:198692374-198703297"

psi_sj_same_start[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),] -> CD45_SJ
cd45 <- cbind(t(CD45_SJ),data_plot[colnames(CD45_SJ),])

tang <- aggregate(cd45[,1:10],by=list(cd45$abb),function(x) {mean(x,na.rm = T)})
row.names(tang) <- as.character(tang[,1])
tang <- tang[,-1]
tang[is.na(tang)] <- -1

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_BULK_intron-centric_PSI.pdf")
pheatmap(tang)
dev.off()


#### version2 from chen
psi_sj_same_start[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),] -> CD45_SJ
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")
cd45 <- cbind(t(CD45_SJ),data_plot[colnames(CD45_SJ),])
cd45_2 <- as.matrix(cd45[,1:10])
row.names(cd45_2) <- cd45$abb

tang2 <- cd45_2
tang2[is.na(tang2)] <- -1

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_BULK_intron-centric_PSI.2.pdf", width = 7, height = 12)
pheatmap(tang2)
dev.off()

tang2 <- cd45_2
tang2[is.na(tang2)] <- -0.1

library(RColorBrewer)
breaksList = seq(-0.1,1,0.1)
library(RColorBrewer)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(10))

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_BULK_intron-centric_PSI.3.pdf", width = 7, height = 12)
pheatmap(tang2, color = colors, breaks = breaksList)
dev.off()



tang2 <- cd45_2
tang2[is.na(tang2)] <- 0

library(RColorBrewer)
breaksList = seq(0,1,0.1)
library(RColorBrewer)
colors = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(10)

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_BULK_intron-centric_PSI.4.pdf", width = 7, height = 12)
pheatmap(tang2, color = colors, breaks = breaksList)
dev.off()


#### Compare with single cell
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
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

CD45_sc_SJ <- psi_sj_same_start_table[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),]
row.names(CD45_sc_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

CD45_sc_SJ <- CD45_sc_SJ[, row.names(Cell_type)]
CD45_sc_SJ[is.na(CD45_sc_SJ)] <- 0

library(RColorBrewer)
breaksList = seq(0,1,0.1)
library(RColorBrewer)
colors = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(10)

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_SC_intron-centric_PSI.pdf", width = 7, height = 12)
pheatmap(t(CD45_sc_SJ), color = colors, breaks = breaksList, annotation_row = Cell_type, show_rownames = F)
dev.off()



CD45_sc_SJ <- psi_sj_same_start_table[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),]
row.names(CD45_sc_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")
CD45_sc_SJ <- CD45_sc_SJ[, row.names(Cell_type)]
CD45_sc_SJ[is.na(CD45_sc_SJ)] <- -0.1

breaksList = c(-0.1, -0.05, 0, seq(0.15,1.05,0.1))
colors = c(grey.colors(1,start = .8), brewer.pal(3,"Greens")[1:2], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(9))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_SC_intron-centric_PSI.2.pdf", width = 7, height = 12)
pheatmap(t(CD45_sc_SJ), color = colors, breaks = breaksList, annotation_row = Cell_type, show_rownames = F)
dev.off()

ann_colors = list(
  Individual = c(UCB1 = "#7570B3", UCB3 = "#E7298A", UCB4 = "#66A61E"),
  CellType = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_SC_intron-centric_PSI.3.pdf", width = 7, height = 6)
pheatmap(t(CD45_sc_SJ), color = colors, breaks = breaksList, annotation_row = Cell_type, annotation_colors = ann_colors, show_rownames = F, cluster_rows = F)
dev.off()







#### SJ =========================================================================================================================================
library("data.table")

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames = list.files(path=mypath, pattern="tab", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp <- unique(fread(x,header=TRUE))
    setkey(tmp,event_id)
    return(tmp)})
  Reduce(function(x,y) {merge(x,y,all=T,by="event_id")}, datalist)
}
path<-file.path("/mnt/data1/projects/UCB/data_ucb/BP_ucb_bulk_SJ")
system.time(te<-multmerge(path))

save(te, file="/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/Bulk_SJ_merged_raw_te.RData")

te_table <- data.frame(as.data.frame(te)[,-1],row.names = te[[1]])

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

FastRemoveMissingValues(te_table)

te_colsum <- colSums(te_table)

for(i in 1:ncol(te_table)){
  te_table[,i] <- (te_table[,i]/te_colsum[i])*mean(te_colsum)
}
save(te_table, file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/Bulk_SJ_merged_raw_te_Lib_normalized.RData")


te_table[c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37),] -> CD45_SJ
CD45_SJ[CD45_SJ<10] <- 0
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

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

table(colSums(CD45_R_r_tab>0))
 0  1  2  3  4  5  6
 7  2  7  9  6 22  6

recalculete <- colnames(CD45_R_r_tab)[c(colSums(CD45_R_r_tab!=0)==4 | colSums(CD45_R_r_tab!=0)==5 | colSums(CD45_R_r_tab!=0)==6)]
dim(CD45_R_r_tab[, recalculete])


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

table(colSums(CD45_R_r_tab>0))
 0  1  2  3  4
 7  2 13 31  6



write.table(CD45_R_r_tab,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/Bulk_CD45_RABC.txt", sep = "\t")
write.table(CD45_SJ,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/Bulk_CD45_SJ.txt", sep = "\t")


cd45 <- cbind(t(CD45_R_r_tab),data_plot[colnames(CD45_R_r_tab),])

chao <- aggregate(cd45[,1:8],by=list(cd45$abb),function(x) {mean(x,na.rm = T)})

row.names(chao) <- as.character(chao[,1])
chao <- chao[,-1]

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_BULK_RABC.pdf")
pheatmap(log10(chao+1))
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/CD45_BULK_RABC2.pdf")
pheatmap(log10(chao+1), cluster_rows = F)
dev.off()






