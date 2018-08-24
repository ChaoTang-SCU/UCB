#### Cell Type Specific AS events
#################################################################################################################################################
#### SE #########################################################################################################################################
#################################################################################################################################################


load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
PSI <- SE_psi[, 13:ncol(SE_psi)]
row.names(PSI) <- SE_psi$loci
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)
PSI <- PSI[, row.names(Cell_type)]
dim(PSI)
# [1] 17670  3039
rs <- rowSums(!is.na(PSI))
sum(rs == 0)
#[1] 53
PSI_tu <- PSI[rs>=1, ]

#Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
#  x <- x[!is.na(x)]
#  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",])] > 0.5)#
#  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)#

#  D = sum(x >= 0.15)
#  n = length(x) - D
#  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",]))
#  pval = phyper(k, D, n, N, lower.tail=FALSE)
#  if(k==0) {
#    adj_pval <- pval
#    } else {
#    adj_pval <- pval * k
#    }
#    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
#  return(enrichment)
#}))#

#Bcell_Sep <- as.data.frame(Bcell_Sep)
#Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
#Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
#Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
#Bcell_Sep <- na.omit(Bcell_Sep)
#dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
##[1] 5436    9
#dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
##[1] 401    9
#hist(Bcell_Sep$pval)#

#dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .2 & Bcell_Sep$Cell_p > 4*Bcell_Sep$Other_p, ])
## [1] 266    9#

#dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 1.5*Bcell_Sep$Other_p & Bcell_Sep$k > 60, ])#
#
#

#PSI_sub <- PSI[row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > Bcell_Sep$Other_p & Bcell_Sep$k > 60, ]), ]#
#

#pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell_Sep_SE_vioplot2.pdf", , width=9, height=6)
#for(i in 1:nrow(PSI_sub)){
#  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
#  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
#        geom_violin()+
#        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
#        geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
#        ggtitle(row.names(PSI_sub[i,]))+
#        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
#                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
#}
#dev.off()#
#
#
#
#
#

#Tcel <- row.names(Cell_type[Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells",])
#Cell_type <- Cell_type[Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells",]#

#PSI_tu <- PSI[rs>=1, Tcel]#

#CD4_Sep <- t(apply(PSI_tu, 1, function(x){
#  x <- x[!is.na(x)]
#  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="CD4+ T Cells",])] > 0.5)#
#  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)#

#  D = sum(x >= 0.15)
#  n = length(x) - D
#  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="CD4+ T Cells",]))
#  pval = phyper(k, D, n, N, lower.tail=FALSE)
#  if(k==0) {
#    adj_pval <- pval
#    } else {
#    adj_pval <- pval * k
#    }
#    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
#  return(enrichment)
#}))#
#

#CD4_Sep <- as.data.frame(CD4_Sep)
#CD4_Sep$Cell_p <- CD4_Sep$k/CD4_Sep$N
#CD4_Sep$Bg_p <- CD4_Sep$D/(CD4_Sep$D + CD4_Sep$n)
#CD4_Sep$Other_p <- (CD4_Sep$D - CD4_Sep$k)/(CD4_Sep$D + CD4_Sep$n - CD4_Sep$N)
#CD4_Sep <- na.omit(CD4_Sep)
#dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
##[1] 3565    9
#dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
##[1] 630    9
#hist(CD4_Sep$pval)#

#dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .2 & CD4_Sep$Cell_p > 4*CD4_Sep$Other_p, ])
## [1] 390    9#

#dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > CD4_Sep$Other_p & CD4_Sep$k > 168, ])
## [1] 69  9#
#

#PSI_sub <- PSI[row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > CD4_Sep$Other_p & CD4_Sep$k > 168, ]), Tcel]#

#pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/CD4_Sep_CD4vsCD8_vioplot.pdf", , width=9, height=6)
#for(i in 1:nrow(PSI_sub)){
#  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
#  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
#        geom_violin()+
#        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
#        geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
#        ggtitle(row.names(PSI_sub[i,]))+
#        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
#                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
#}
#dev.off()#
#
#
#
#
#
#

#CD8_Sep <- t(apply(PSI_tu, 1, function(x){
#  x <- x[!is.na(x)]
#  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="CD8+ T Cells",])] > 0.5)#
#  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)#

#  D = sum(x >= 0.15)
#  n = length(x) - D
#  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="CD8+ T Cells",]))
#  pval = phyper(k, D, n, N, lower.tail=FALSE)
#  if(k==0) {
#    adj_pval <- pval
#    } else {
#    adj_pval <- pval * k
#    }
#    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
#  return(enrichment)
#}))#
#

#CD8_Sep <- as.data.frame(CD8_Sep)
#CD8_Sep$Cell_p <- CD8_Sep$k/CD8_Sep$N
#CD8_Sep$Bg_p <- CD8_Sep$D/(CD8_Sep$D + CD8_Sep$n)
#CD8_Sep$Other_p <- (CD8_Sep$D - CD8_Sep$k)/(CD8_Sep$D + CD8_Sep$n - CD8_Sep$N)
#CD8_Sep <- na.omit(CD8_Sep)
#dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
##[1] 5318    9
#dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
##[1] 403    9
#hist(CD8_Sep$pval)#

#dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .2 & CD8_Sep$Cell_p > 4*CD8_Sep$Other_p, ])
## [1] 256    9#

#dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > CD8_Sep$Other_p & CD8_Sep$k > 168, ])
## [1] 142  9#

#PSI_sub <- PSI[row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > CD8_Sep$Other_p & CD8_Sep$k > 168, ]), Tcel]#

#pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/CD8_Sep_CD8vsCD4_vioplot.pdf", , width=9, height=6)
#for(i in 1:nrow(PSI_sub)){
#  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
#  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
#        geom_violin()+
#        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
#        geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
#        ggtitle(row.names(PSI_sub[i,]))+
#        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
#                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
#}
#dev.off()


#### CD4 VS. CD8 ===================


PSI_tu_CD4 <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "CD4+ T Cells",])]
PSI_tu_CD8 <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "CD8+ T Cells",])]

psi_mean_4 <- rowMeans(PSI_tu_CD4, na.rm=T)
psi_mean_8 <- rowMeans(PSI_tu_CD8, na.rm=T)

psi_number_4 <- rowSums(!is.na(PSI_tu_CD4))
psi_number_8 <- rowSums(!is.na(PSI_tu_CD8))

test <- data.frame(cbind(psi_mean_4,psi_mean_8,psi_number_4,psi_number_8))

sum(psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1)
#[1] 3842
sum(psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1 )
test <- test[psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1, ]
dim(na.omit(test))
#[1] 3842    4
max(test$psi_mean_4 - test$psi_mean_8)
#[1] 0.1849871
test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]

Cell_type <- Cell_type[Cell_type$CellType == "CD4+ T Cells" | Cell_type$CellType == "CD8+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/CD4vsCD8_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
        geom_violin()+
        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
        #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
        ggtitle(row.names(PSI_sub[i,]))+
        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/CD4vsCD8_featureplot.pdf", , width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()
  
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/CD4vsCD8_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()




#### B Cell VS. CD8 ===================



Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

PSI_tu_BCell <- PSI[row.names(PSI_tu), row.names(Cell_type[Cell_type$CellType == "B Cells",])]
psi_mean_B <- rowMeans(PSI_tu_BCell, na.rm=T)
psi_number_B <- rowSums(!is.na(PSI_tu_BCell))

test <- data.frame(cbind(psi_mean_B,psi_mean_8,psi_number_B,psi_number_8))
test <- test[psi_number_B > ncol(PSI_tu_BCell)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1, ]
dim(test)
#[1] 3510    4
dim(na.omit(test))
#[1] 3510    4
abs(test$psi_mean_B - test$psi_mean_8)

test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]


Cell_type <- Cell_type[Cell_type$CellType == "B Cells" | Cell_type$CellType == "CD8+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/BcellvsCD8_vioplot.pdf", width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
        geom_violin()+
        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
        #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
        ggtitle(row.names(PSI_sub[i,]))+
        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/BcellvsCD8_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/BcellvsCD8_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()




#### B Cell VS. CD4 ===================



Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

test <- data.frame(cbind(psi_mean_B,psi_mean_4,psi_number_B,psi_number_4))
test <- test[psi_number_B > ncol(PSI_tu_BCell)*.1 & psi_number_4 > ncol(PSI_tu_CD4)*.1, ]
dim(test)
#[1] 3360    4
dim(na.omit(test))
#[1] 3360    4
abs(test$psi_mean_B - test$psi_mean_4)

test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]


Cell_type <- Cell_type[Cell_type$CellType == "B Cells" | Cell_type$CellType == "CD4+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/BcellvsCD4_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
          geom_violin()+
          #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
          #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
          ggtitle(row.names(PSI_sub[i,]))+
          scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()



tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/BcellvsCD4_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/BcellvsCD4_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()



PSI_tu_BCell["1:198692374-198694017-198696910-198699563", ]
sashimi(junction = "1:198692374-198699563", cell = c("UCB4.01164", "UCB4.01435", "UCB4.01408", "UCB4.01357"), 
  outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/")

PSI_tu_BCell["1:198692374-198699563-198699705-198703297", ]
sashimi(junction = "1:198692374-198699563", cell = c("UCB4.01164", "UCB4.01190", "UCB4.01226", "UCB4.01248"), 
  outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/")

PSI_tu_BCell["1:198699705-198702386-198702531-198703297", ]
sashimi(junction = "1:198699705-198703297", cell = c("UCB4.01435", "UCB4.01327", "UCB4.01357"), 
  outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/")



#################################################################################################################################################
#### A3SS & A5SS ################################################################################################################################
#################################################################################################################################################



load("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi_new.RData")
load("/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_psi_new.RData")

PSI3 <- A3SS_psi[, 15:ncol(A3SS_psi)]
PSI5 <- A5SS_psi[, 15:ncol(A5SS_psi)]
row.names(PSI3) <- A3SS_psi$as2
row.names(PSI5) <- A5SS_psi$as2
PSI <- rbind(PSI3, PSI5)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

PSI <- PSI[, row.names(Cell_type)]
dim(PSI)

rs <- rowSums(!is.na(PSI))
sum(rs == 0)
# [1] 13

PSI_tu <- PSI[rs>=1, ]


PSI_tu_CD4 <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "CD4+ T Cells",])]
PSI_tu_CD8 <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "CD8+ T Cells",])]
PSI_tu_B <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "B Cells",])]

psi_mean_4 <- rowMeans(PSI_tu_CD4, na.rm=T)
psi_mean_8 <- rowMeans(PSI_tu_CD8, na.rm=T)
psi_mean_B <- rowMeans(PSI_tu_B, na.rm=T)

psi_number_4 <- rowSums(!is.na(PSI_tu_CD4))
psi_number_8 <- rowSums(!is.na(PSI_tu_CD8))
psi_number_B <- rowSums(!is.na(PSI_tu_B))


#### CD4 VS. CD8 ===================


test <- data.frame(cbind(psi_mean_4,psi_mean_8,psi_number_4,psi_number_8))

sum(psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1)
#[1] 608
test <- test[psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1, ]
dim(na.omit(test))
#[1] 608    4
max(test$psi_mean_4 - test$psi_mean_8)
#[1] 0.1381762
test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[Cell_type$CellType == "CD4+ T Cells" | Cell_type$CellType == "CD8+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/CD4vsCD8_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
        geom_violin()+
        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
        #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
        ggtitle(row.names(PSI_sub[i,]))+
        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/CD4vsCD8_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/CD4vsCD8_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()




#### B Cell VS. CD4 ===================



test <- data.frame(cbind(psi_mean_B,psi_mean_4,psi_number_B,psi_number_4))
test <- test[psi_number_B > ncol(PSI_tu_BCell)*.1 & psi_number_4 > ncol(PSI_tu_CD4)*.1, ]
dim(test)
#[1] 507    4
abs(test$psi_mean_B - test$psi_mean_4)

test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[Cell_type$CellType == "B Cells" | Cell_type$CellType == "CD4+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/BcellvsCD4_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
          geom_violin()+
          #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
          #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
          ggtitle(row.names(PSI_sub[i,]))+
          scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/BcellvsCD4_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/BcellvsCD4_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()




PSI_tu_B["X:19746442-19836124_X:19746442-19836150",]


sashimi(junction = "X:19746442-19836150", cell = c("UCB4.01141", "UCB4.00386"), 
  outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/")

sashimi(junction = "X:19836124-19836150", cell = c("UCB4.01141", "UCB4.00386"), 
  outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/")



#### B Cell VS. CD8 ===================



test <- data.frame(cbind(psi_mean_B,psi_mean_8,psi_number_B,psi_number_8))
test <- test[psi_number_B > ncol(PSI_tu_BCell)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1, ]
dim(test)
#[1] 543    4
abs(test$psi_mean_B - test$psi_mean_8)

test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[Cell_type$CellType == "B Cells" | Cell_type$CellType == "CD8+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/BcellvsCD8_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
        geom_violin()+
        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
        #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
        ggtitle(row.names(PSI_sub[i,]))+
        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/BcellvsCD8_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/BcellvsCD8_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()



#################################################################################################################################################
#### AFE & ALE ##################################################################################################################################
#################################################################################################################################################



load("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi_new.RData")
load("/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_psi_new.RData")

PSIf <- AFE_psi[, 12:ncol(AFE_psi)]
PSIl <- ALE_psi[, 12:ncol(ALE_psi)]
row.names(PSIf) <- AFE_psi$as2
row.names(PSIl) <- ALE_psi$as2
PSI <- rbind(PSIf, PSIl)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

PSI <- PSI[, row.names(Cell_type)]
dim(PSI)
# [1] 4263  3039

rs <- rowSums(!is.na(PSI))
sum(rs == 0)
# [1] 3

PSI_tu <- PSI[rs>=1, ]


PSI_tu_CD4 <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "CD4+ T Cells",])]
PSI_tu_CD8 <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "CD8+ T Cells",])]
PSI_tu_B <- PSI_tu[, row.names(Cell_type[Cell_type$CellType == "B Cells",])]

psi_mean_4 <- rowMeans(PSI_tu_CD4, na.rm=T)
psi_mean_8 <- rowMeans(PSI_tu_CD8, na.rm=T)
psi_mean_B <- rowMeans(PSI_tu_B, na.rm=T)

psi_number_4 <- rowSums(!is.na(PSI_tu_CD4))
psi_number_8 <- rowSums(!is.na(PSI_tu_CD8))
psi_number_B <- rowSums(!is.na(PSI_tu_B))


#### CD4 VS. CD8 ===================


test <- data.frame(cbind(psi_mean_4,psi_mean_8,psi_number_4,psi_number_8))

sum(psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1)
#[1] 1202
test <- test[psi_number_4 > ncol(PSI_tu_CD4)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1, ]
dim(na.omit(test))
#[1] 1202    4
max(test$psi_mean_4 - test$psi_mean_8)
#[1] 0.1381762
test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[Cell_type$CellType == "CD4+ T Cells" | Cell_type$CellType == "CD8+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/CD4vsCD8_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
        geom_violin()+
        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
        #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
        ggtitle(row.names(PSI_sub[i,]))+
        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/CD4vsCD8_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_4 - test$psi_mean_8)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/CD4vsCD8_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()


#### B Cell VS. CD4 ===================



test <- data.frame(cbind(psi_mean_B,psi_mean_4,psi_number_B,psi_number_4))
test <- test[psi_number_B > ncol(PSI_tu_BCell)*.1 & psi_number_4 > ncol(PSI_tu_CD4)*.1, ]
dim(test)
#[1] 1063    4
abs(test$psi_mean_B - test$psi_mean_4)

test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[Cell_type$CellType == "B Cells" | Cell_type$CellType == "CD4+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/BcellvsCD4_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
          geom_violin()+
          #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
          #geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
          ggtitle(row.names(PSI_sub[i,]))+
          scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/BcellvsCD4_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_4)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/BcellvsCD4_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()



#### B Cell VS. CD8 ===================



test <- data.frame(cbind(psi_mean_B,psi_mean_8,psi_number_B,psi_number_8))
test <- test[psi_number_B > ncol(PSI_tu_BCell)*.1 & psi_number_8 > ncol(PSI_tu_CD8)*.1, ]
dim(test)
#[1] 1116    4
dim(na.omit(test))
#[1] 1116    4
abs(test$psi_mean_B - test$psi_mean_8)

test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[Cell_type$CellType == "B Cells" | Cell_type$CellType == "CD8+ T Cells", ]
Cell_type$CellType <- as.factor(as.character(Cell_type$CellType))
PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/BcellvsCD8_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI_sub)){
  plot_tab <- data.frame(psi = as.numeric(PSI_sub[i,]), Cell = Cell_type$CellType)
  print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
        geom_violin()+
        #geom_point(color="blue", alpha=.5, na.rm = TRUE)+
        geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
        ggtitle(row.names(PSI_sub[i,]))+
        scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
                        table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()


tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/BcellvsCD8_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()

PSI_sub <- PSI[row.names(test[tail(order(abs(test$psi_mean_B - test$psi_mean_8)),10),]), ]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/BcellvsCD8_featureplot2.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "#00FFFF", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()


PSI_tu["13:40982283-40982857_13:40982283-41019227",]


sashimi(junction = "13:40982283-41019227", cell = c("UCB4.00409", "UCB4.00842", "UCB4.00457", "UCB4.00549"), 
  outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/")




