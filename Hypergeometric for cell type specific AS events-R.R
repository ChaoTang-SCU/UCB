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
Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
  x <- x[!is.na(x)]
  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",])] > 0.15)#
  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)

  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",]))
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
Bcell_Sep -> Bcell_Sep2
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
#[1] 7728    9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
#[1] 2463    9
hist(Bcell_Sep$pval)


PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 0
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="B Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
#[1] 385   9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
#[1] 219   9
hist(Bcell_Sep$pval)



CD4_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD4+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

CD4_Sep <- as.data.frame(CD4_Sep)
CD4_Sep$Cell_p <- CD4_Sep$k/CD4_Sep$N
CD4_Sep$Bg_p <- CD4_Sep$D/(CD4_Sep$D + CD4_Sep$n)
CD4_Sep$Other_p <- (CD4_Sep$D - CD4_Sep$k)/(CD4_Sep$D + CD4_Sep$n - CD4_Sep$N)

head(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ],20)
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
# [1] 237   9
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
# [1] 27  9
hist(CD4_Sep$pval)


CD8_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

CD8_Sep <- as.data.frame(CD8_Sep)
CD8_Sep$Cell_p <- CD8_Sep$k/CD8_Sep$N
CD8_Sep$Bg_p <- CD8_Sep$D/(CD8_Sep$D + CD8_Sep$n)
CD8_Sep$Other_p <- (CD8_Sep$D - CD8_Sep$k)/(CD8_Sep$D + CD8_Sep$n - CD8_Sep$N)

head(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ],20)
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
# [1] 434    9
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
# [1] 23   9
hist(CD8_Sep$pval)



MK_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Megakaryocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="Megakaryocytes")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

MK_Sep <- as.data.frame(MK_Sep)
MK_Sep$Cell_p <- MK_Sep$k/MK_Sep$N
MK_Sep$Bg_p <- MK_Sep$D/(MK_Sep$D + MK_Sep$n)
MK_Sep$Other_p <- (MK_Sep$D - MK_Sep$k)/(MK_Sep$D + MK_Sep$n - MK_Sep$N)
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
# [1] 309    9
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
# [1] 302    9
hist(MK_Sep$pval)



Mono_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Monocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="Monocytes")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Mono_Sep <- as.data.frame(Mono_Sep)
Mono_Sep$Cell_p <- Mono_Sep$k/Mono_Sep$N
Mono_Sep$Bg_p <- Mono_Sep$D/(Mono_Sep$D + Mono_Sep$n)
Mono_Sep$Other_p <- (Mono_Sep$D - Mono_Sep$k)/(Mono_Sep$D + Mono_Sep$n - Mono_Sep$N)
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
# [1] 260    9
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
# [1] 251    9
hist(Mono_Sep$pval)



NK_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="NK cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="NK cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

NK_Sep <- as.data.frame(NK_Sep)
NK_Sep$Cell_p <- NK_Sep$k/NK_Sep$N
NK_Sep$Bg_p <- NK_Sep$D/(NK_Sep$D + NK_Sep$n)
NK_Sep$Other_p <- (NK_Sep$D - NK_Sep$k)/(NK_Sep$D + NK_Sep$n - NK_Sep$N)
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
# [1] 234    9
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
# [1] 167   9
hist(NK_Sep$pval)



Tcell_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Tcell_Sep <- as.data.frame(Tcell_Sep)
Tcell_Sep$Cell_p <- Tcell_Sep$k/Tcell_Sep$N
Tcell_Sep$Bg_p <- Tcell_Sep$D/(Tcell_Sep$D + Tcell_Sep$n)
Tcell_Sep$Other_p <- (Tcell_Sep$D - Tcell_Sep$k)/(Tcell_Sep$D + Tcell_Sep$n - Tcell_Sep$N)

head(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ],20)
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
# [1] 444    9
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
# [1] 226   9
hist(Tcell_Sep$pval)


save(Bcell_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep, Tcell_Sep, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_sep_SE.RData")



sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])


length(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK))
# [1] 1215
length(unique(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)))
# [1] 1087

c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK) -> juncs


PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 1215 3039


te_sub[is.na(te_sub)] <- -0.2

cells <- c("B Cells", "T Cells", "CD4+ T Cells", "CD8+ T Cells", "Megakaryocytes", "Monocytes", "NK cells")
sj_anno <- matrix(rep(cells, c(length(sj_Bcell),length(sj_Tcell),length(sj_CD4),length(sj_CD8),length(sj_MK),length(sj_Mono),length(sj_NK))), ncol = 1)
rownames(sj_anno) <- row.names(te_sub)
#sj_anno <- data.frame(sj_anno[juncs,])
colnames(sj_anno) <- "Cell"
sj_anno <- as.data.frame(sj_anno)
sj_anno$Cell <- factor(sj_anno$Cell, levels = cells)

dim(sj_anno)
# [1] 1215    1
identical(row.names(sj_anno), row.names(te_sub))
# [1] TRUE


library(RColorBrewer)
library(pheatmap)
breaksList = seq(-0.2,1,0.1)
library(RColorBrewer)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(11))
#colors = c(rep(brewer.pal(3,"Greens")[1], 4), brewer.pal(,"Reds")[4:9])
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(11))

brewer.pal(6,"BrBG")
ann_colors = list(
  Individual = c(UCB1 = "#FACA13", UCB3 = "#AACBE2", UCB4 = "#F96566"),
  CellType = c("#FF7F00","#4DAF4A","#377EB8","#F781BF","#A65628","#984EA3"),
  Cell = c("#FF7F00","#E41A1C","#4DAF4A","#377EB8","#F781BF","#A65628","#984EA3")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))
names(ann_colors[[3]]) <- c(names(table(Cell_type$CellType))[1],"T Cells", names(table(Cell_type$CellType))[-1])

cumsum(as.numeric(table(sj_anno$Cell)))[1:6]


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_specific_SE_heatmap.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()



juncs <- names(which(table(juncs)==1))
PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 972 3039
te_sub[is.na(te_sub)] <- -0.2

#sj_anno2 <- subset.data.frame(sj_anno, row.names(sj_anno) %in% row.names(te_sub))
sj_anno2 <- data.frame(Cell = sj_anno[row.names(te_sub), ], row.names = row.names(te_sub))
identical(row.names(sj_anno2), row.names(te_sub))

te_sub <- te_sub[row.names(sj_anno2)[order(sj_anno2$Cell)], ]

sj_anno2 <- data.frame(Cell = sj_anno2[order(sj_anno2$Cell),], row.names = row.names(sj_anno2)[order(sj_anno2$Cell)])

cumsum(as.numeric(table(sj_anno2$Cell)))[1:6]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_specific_SE_heatmap2.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno2$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()



PSI_sub <- PSI[sj_Bcell, ]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell_Sep_SE_vioplot.pdf",  width=9, height=6)
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




















which(row.names(te_sub) %in% SE_psi[SE_psi$SE_type == "novel_exon-skipping", "loci"])
row.names(te_sub)[which(row.names(te_sub) %in% SE_psi[SE_psi$SE_type == "novel_exon-skipping", "loci"])]
#[1] "12:51284788-51287883-51288030-51288107"
#[2] "12:51284788-51287883-51288033-51288107"
sj_anno[which(row.names(te_sub) %in% SE_psi[SE_psi$SE_type == "novel_exon-skipping", "loci"]), ]
#[1] Megakaryocytes Megakaryocytes
table(SE_psi[SE_psi$loci %in% row.names(te_sub), "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other
#                   818                     2                   267


load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_tu <- te_sub[which(row.names(te_sub) %in% SE_psi[SE_psi$SE_type == "novel_exon-skipping", "loci"]),]

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

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab1
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[2,])), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[2])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/MK_specific_novel_SE_heatmap.pdf")
f1
f2
dev.off()


sashimi <- function(junction, cell, left_expand = 1000, right_expand = 1000, tsv = "/mnt/data5/BGI/UCB/tangchao/data/BAM_tsv/", outDir = "~/", min = 10, ann_height = 4){
    if(dir.exists(outDir) == FALSE) dir.create(outDir)
    dir.create("/home/tangchao/test2")
    Coordinates <- data.frame(do.call(rbind,strsplit(junction, split = "[:-]")),stringsAsFactors=F)
    Coordinates$X2 <- abs(as.numeric(Coordinates$X2) - left_expand)
    Coordinates$X3 <- abs(as.numeric(Coordinates$X3) + right_expand)
    window <- paste(Coordinates$X1, ":", Coordinates$X2, "-", Coordinates$X3, sep = "")
    if(length(cell) == 1){
        Dir_out <- paste(outDir, cell, "_", junction, sep = "")
        work <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi.sh ", window, cell, Dir_out, min, ann_height, "')")
        eval(parse(text = work))
    }else{
        Dir_out <- paste(outDir, junction, sep = "")
        tmpdir <- "/home/tangchao/test2"
        work1 <- paste("system('", "cat ", paste(tsv,cell,".tsv",sep = "", collapse = " "), " > ", tmpdir,"/tmp.tsv", "')", sep = "")
        work2 <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi2.sh ", window, paste(tmpdir,"/tmp.tsv", sep = ""), Dir_out, min, ann_height, "')")
        work3 <- paste("system('", "rm",paste(tmpdir,"/tmp.tsv", "')", sep = ""))
        eval(parse(text = work1))
        eval(parse(text = work2))
        eval(parse(text = work3))
    }
}


PSI[row.names(te_sub)[which(row.names(te_sub) %in% SE_psi[SE_psi$SE_type == "novel_exon-skipping", "loci"])], row.names(Cell_type[Cell_type$CellType == "Megakaryocytes",])]


sashimi(junction = "12:51284788-51288107", cell = c("UCB3.00069", "UCB3.01025", "UCB3.02215", "UCB3.02111"), 
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
# [1] 6086  3039

rs <- rowSums(!is.na(PSI))
sum(rs == 0)
# [1] 13

PSI_tu <- PSI[rs>=1, ]
Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
  x <- x[!is.na(x)]
  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",])] > 0.15)#
  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)

  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",]))
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
Bcell_Sep -> Bcell_Sep2
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
#[1] 3193    9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
#[1] 1049    9
hist(Bcell_Sep$pval)


PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 0
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="B Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
#[1] 89   9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
#[1] 62   9
hist(Bcell_Sep$pval)



CD4_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD4+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

CD4_Sep <- as.data.frame(CD4_Sep)
CD4_Sep$Cell_p <- CD4_Sep$k/CD4_Sep$N
CD4_Sep$Bg_p <- CD4_Sep$D/(CD4_Sep$D + CD4_Sep$n)
CD4_Sep$Other_p <- (CD4_Sep$D - CD4_Sep$k)/(CD4_Sep$D + CD4_Sep$n - CD4_Sep$N)

head(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ],20)
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
# [1] 46   9
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
# [1] 9  9
hist(CD4_Sep$pval)


CD8_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

CD8_Sep <- as.data.frame(CD8_Sep)
CD8_Sep$Cell_p <- CD8_Sep$k/CD8_Sep$N
CD8_Sep$Bg_p <- CD8_Sep$D/(CD8_Sep$D + CD8_Sep$n)
CD8_Sep$Other_p <- (CD8_Sep$D - CD8_Sep$k)/(CD8_Sep$D + CD8_Sep$n - CD8_Sep$N)

head(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ],20)
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
# [1] 68    9
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
# [1] 7   9
hist(CD8_Sep$pval)



MK_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Megakaryocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="Megakaryocytes")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

MK_Sep <- as.data.frame(MK_Sep)
MK_Sep$Cell_p <- MK_Sep$k/MK_Sep$N
MK_Sep$Bg_p <- MK_Sep$D/(MK_Sep$D + MK_Sep$n)
MK_Sep$Other_p <- (MK_Sep$D - MK_Sep$k)/(MK_Sep$D + MK_Sep$n - MK_Sep$N)
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
# [1] 136    9
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
# [1] 127    9
hist(MK_Sep$pval)



Mono_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Monocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="Monocytes")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Mono_Sep <- as.data.frame(Mono_Sep)
Mono_Sep$Cell_p <- Mono_Sep$k/Mono_Sep$N
Mono_Sep$Bg_p <- Mono_Sep$D/(Mono_Sep$D + Mono_Sep$n)
Mono_Sep$Other_p <- (Mono_Sep$D - Mono_Sep$k)/(Mono_Sep$D + Mono_Sep$n - Mono_Sep$N)
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
# [1] 75    9
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
# [1] 68    9
hist(Mono_Sep$pval)



NK_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="NK cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="NK cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

NK_Sep <- as.data.frame(NK_Sep)
NK_Sep$Cell_p <- NK_Sep$k/NK_Sep$N
NK_Sep$Bg_p <- NK_Sep$D/(NK_Sep$D + NK_Sep$n)
NK_Sep$Other_p <- (NK_Sep$D - NK_Sep$k)/(NK_Sep$D + NK_Sep$n - NK_Sep$N)
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
# [1] 81    9
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
# [1] 61   9
hist(NK_Sep$pval)



Tcell_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Tcell_Sep <- as.data.frame(Tcell_Sep)
Tcell_Sep$Cell_p <- Tcell_Sep$k/Tcell_Sep$N
Tcell_Sep$Bg_p <- Tcell_Sep$D/(Tcell_Sep$D + Tcell_Sep$n)
Tcell_Sep$Other_p <- (Tcell_Sep$D - Tcell_Sep$k)/(Tcell_Sep$D + Tcell_Sep$n - Tcell_Sep$N)

head(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ],20)
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
# [1] 68    9
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
# [1] 32   9
hist(Tcell_Sep$pval)


save(Bcell_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep, Tcell_Sep, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/All_Cell_sep_ASS.RData")



sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])



length(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK))
# [1] 366
length(unique(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)))
# [1] 315

c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK) -> juncs


PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 366 3039


te_sub[is.na(te_sub)] <- -0.2

cells <- c("B Cells", "T Cells", "CD4+ T Cells", "CD8+ T Cells", "Megakaryocytes", "Monocytes", "NK cells")
sj_anno <- matrix(rep(cells, c(length(sj_Bcell),length(sj_Tcell),length(sj_CD4),length(sj_CD8),length(sj_MK),length(sj_Mono),length(sj_NK))), ncol = 1)
rownames(sj_anno) <- row.names(te_sub)
#sj_anno <- data.frame(sj_anno[juncs,])
colnames(sj_anno) <- "Cell"
sj_anno <- as.data.frame(sj_anno)
sj_anno$Cell <- factor(sj_anno$Cell, levels = cells)

dim(sj_anno)
# [1] 366    1
identical(row.names(sj_anno), row.names(te_sub))
# [1] TRUE


library(RColorBrewer)
library(pheatmap)
breaksList = seq(-0.2,1,0.1)
library(RColorBrewer)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(11))
#colors = c(rep(brewer.pal(3,"Greens")[1], 4), brewer.pal(,"Reds")[4:9])

brewer.pal(6,"BrBG")
ann_colors = list(
  Individual = c(UCB1 = "#7570B3", UCB3 = "#E7298A", UCB4 = "#66A61E"),
  CellType = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E"),
  Cell = c("#8C510A","#7570B3","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))
names(ann_colors[[3]]) <- c(names(table(Cell_type$CellType))[1],"T Cells", names(table(Cell_type$CellType))[-1])

cumsum(as.numeric(table(sj_anno$Cell)))[1:6]


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/All_Cell_specific_ASS_heatmap.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()



juncs <- names(which(table(juncs)==1))
PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 972 3039
te_sub[is.na(te_sub)] <- -0.2

#sj_anno2 <- subset.data.frame(sj_anno, row.names(sj_anno) %in% row.names(te_sub))
sj_anno2 <- data.frame(Cell = sj_anno[row.names(te_sub), ], row.names = row.names(te_sub))
identical(row.names(sj_anno2), row.names(te_sub))

te_sub <- te_sub[row.names(sj_anno2)[order(sj_anno2$Cell)], ]

sj_anno2 <- data.frame(Cell = sj_anno2[order(sj_anno2$Cell),], row.names = row.names(sj_anno2)[order(sj_anno2$Cell)])

cumsum(as.numeric(table(sj_anno2$Cell)))[1:6]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/All_Cell_specific_ASS_heatmap2.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno2$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()











#################################################################################################################################################
#### AFE & ALE ################################################################################################################################
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

Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
  x <- x[!is.na(x)]
  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",])] > 0.15)#
  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)

  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",]))
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
Bcell_Sep -> Bcell_Sep2
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
#[1] 1052    9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
#[1] 422    9
hist(Bcell_Sep$pval)


PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 0
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

Bcell_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="B Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
#[1] 59   9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
#[1] 37   9
hist(Bcell_Sep$pval)



CD4_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD4+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

CD4_Sep <- as.data.frame(CD4_Sep)
CD4_Sep$Cell_p <- CD4_Sep$k/CD4_Sep$N
CD4_Sep$Bg_p <- CD4_Sep$D/(CD4_Sep$D + CD4_Sep$n)
CD4_Sep$Other_p <- (CD4_Sep$D - CD4_Sep$k)/(CD4_Sep$D + CD4_Sep$n - CD4_Sep$N)

head(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ],20)
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
# [1] 37   9
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
# [1] 2  9
hist(CD4_Sep$pval)


CD8_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

CD8_Sep <- as.data.frame(CD8_Sep)
CD8_Sep$Cell_p <- CD8_Sep$k/CD8_Sep$N
CD8_Sep$Bg_p <- CD8_Sep$D/(CD8_Sep$D + CD8_Sep$n)
CD8_Sep$Other_p <- (CD8_Sep$D - CD8_Sep$k)/(CD8_Sep$D + CD8_Sep$n - CD8_Sep$N)

head(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ],20)
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
# [1] 95    9
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
# [1] 8   9
hist(CD8_Sep$pval)



MK_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Megakaryocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="Megakaryocytes")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

MK_Sep <- as.data.frame(MK_Sep)
MK_Sep$Cell_p <- MK_Sep$k/MK_Sep$N
MK_Sep$Bg_p <- MK_Sep$D/(MK_Sep$D + MK_Sep$n)
MK_Sep$Other_p <- (MK_Sep$D - MK_Sep$k)/(MK_Sep$D + MK_Sep$n - MK_Sep$N)
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
# [1] 83    9
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
# [1] 82    9
hist(MK_Sep$pval)



Mono_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Monocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="Monocytes")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Mono_Sep <- as.data.frame(Mono_Sep)
Mono_Sep$Cell_p <- Mono_Sep$k/Mono_Sep$N
Mono_Sep$Bg_p <- Mono_Sep$D/(Mono_Sep$D + Mono_Sep$n)
Mono_Sep$Other_p <- (Mono_Sep$D - Mono_Sep$k)/(Mono_Sep$D + Mono_Sep$n - Mono_Sep$N)
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
# [1] 76    9
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
# [1] 74    9
hist(Mono_Sep$pval)



NK_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="NK cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="NK cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

NK_Sep <- as.data.frame(NK_Sep)
NK_Sep$Cell_p <- NK_Sep$k/NK_Sep$N
NK_Sep$Bg_p <- NK_Sep$D/(NK_Sep$D + NK_Sep$n)
NK_Sep$Other_p <- (NK_Sep$D - NK_Sep$k)/(NK_Sep$D + NK_Sep$n - NK_Sep$N)
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
# [1] 53    9
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
# [1] 38   9
hist(NK_Sep$pval)



Tcell_Sep <- t(apply(PSI_tu, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))

Tcell_Sep <- as.data.frame(Tcell_Sep)
Tcell_Sep$Cell_p <- Tcell_Sep$k/Tcell_Sep$N
Tcell_Sep$Bg_p <- Tcell_Sep$D/(Tcell_Sep$D + Tcell_Sep$n)
Tcell_Sep$Other_p <- (Tcell_Sep$D - Tcell_Sep$k)/(Tcell_Sep$D + Tcell_Sep$n - Tcell_Sep$N)

head(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ],20)
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
# [1] 91    9
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
# [1] 40   9
hist(Tcell_Sep$pval)


save(Bcell_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep, Tcell_Sep, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/All_Cell_sep_AFLE.RData")


sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])



length(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK))
# [1] 281
length(unique(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)))
# [1] 253

c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK) -> juncs


PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 281 3039


te_sub[is.na(te_sub)] <- -0.2
cells <- c("B Cells", "T Cells", "CD4+ T Cells", "CD8+ T Cells", "Megakaryocytes", "Monocytes", "NK cells")
sj_anno <- matrix(rep(cells, c(length(sj_Bcell),length(sj_Tcell),length(sj_CD4),length(sj_CD8),length(sj_MK),length(sj_Mono),length(sj_NK))), ncol = 1)
rownames(sj_anno) <- row.names(te_sub)
#sj_anno <- data.frame(sj_anno[juncs,])
colnames(sj_anno) <- "Cell"
sj_anno <- as.data.frame(sj_anno)
sj_anno$Cell <- factor(sj_anno$Cell, levels = cells)

dim(sj_anno)
# [1] 281    1
identical(row.names(sj_anno), row.names(te_sub))
# [1] TRUE


library(RColorBrewer)
library(pheatmap)
breaksList = seq(-0.2,1,0.1)
library(RColorBrewer)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(11))
#colors = c(rep(brewer.pal(3,"Greens")[1], 4), brewer.pal(,"Reds")[4:9])

brewer.pal(6,"BrBG")
ann_colors = list(
  Individual = c(UCB1 = "#7570B3", UCB3 = "#E7298A", UCB4 = "#66A61E"),
  CellType = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E"),
  Cell = c("#8C510A","#7570B3","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))
names(ann_colors[[3]]) <- c(names(table(Cell_type$CellType))[1],"T Cells", names(table(Cell_type$CellType))[-1])

cumsum(as.numeric(table(sj_anno$Cell)))[1:6]


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/All_Cell_specific_AFLE_heatmap.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()



juncs <- names(which(table(juncs)==1))
PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 972 3039
te_sub[is.na(te_sub)] <- -0.2

#sj_anno2 <- subset.data.frame(sj_anno, row.names(sj_anno) %in% row.names(te_sub))
sj_anno2 <- data.frame(Cell = sj_anno[row.names(te_sub), ], row.names = row.names(te_sub))
identical(row.names(sj_anno2), row.names(te_sub))

te_sub <- te_sub[row.names(sj_anno2)[order(sj_anno2$Cell)], ]

sj_anno2 <- data.frame(Cell = sj_anno2[order(sj_anno2$Cell),], row.names = row.names(sj_anno2)[order(sj_anno2$Cell)])

cumsum(as.numeric(table(sj_anno2$Cell)))[1:6]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/All_Cell_specific_AFELE_heatmap2.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno2$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()




#################################################################################################################################################
#### MXE ########################################################################################################################################
#################################################################################################################################################


load("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi_new.RData")
PSI <- MXE_psi[, 21:ncol(MXE_psi)]
MXE_psi$ID <- paste("MXE", 1:nrow(MXE_psi), sep="")
row.names(PSI) <- MXE_psi$ID
rs <- rowSums(!is.na(PSI))
sum(rs == 0)
#[1] 40

PSI <- PSI[rowSums(!is.na(PSI))>1, row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Cell_Sep/All_MXE_vioplot.pdf", , width=9, height=6)
for(i in 1:nrow(PSI)){
	plot_tab <- data.frame(psi = as.numeric(PSI[i,]), Cell = Cell_type$CellType)
	print(ggplot(data = plot_tab, aes(x = Cell, y = psi, fill = Cell))+
				geom_violin()+
				#geom_point(color="blue", alpha=.5, na.rm = TRUE)+
				geom_jitter(height = 0, color="blue", alpha=.5, na.rm = TRUE)+
				ggtitle(row.names(PSI[i,]))+
				scale_x_discrete(labels=paste(names(table(na.omit(plot_tab)$Cell)), "\n", 
											  table(na.omit(plot_tab)$Cell), "/", table(plot_tab$Cell), sep="")))
}
dev.off()








