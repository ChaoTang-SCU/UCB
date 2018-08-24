#### AS Sorting out

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

#### SE =========================================================================================================================================

load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
PSI <- SE_psi[, row.names(Cell_type)]
row.names(PSI) <- SE_psi$loci
rs <- rowSums(!is.na(PSI))
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0.0    17.0    76.0   256.3   273.8  3039.0
sum(rs==0)
#[1] 53
sum(rs==1)
#[1] 513
PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)

table(SE_psi[SE_psi$loci %in% names(rs)[rs>=1], "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5272                  1859                  2684
table(SE_psi[SE_psi$loci %in% names(rs)[rs==0], "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5088                  1529                  1185

SE_psi <- SE_psi[SE_psi$loci %in% row.names(PSI), ]
dim(SE_psi)
#[1] 17617  3586
table(SE_psi$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                 10360                  3388                  3869

table(SE_psi[SE_psi$loci %in% names(rs)[rs>=1], "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other                                     
#                  5272                  1859                  2684                                     
table(SE_psi[SE_psi$loci %in% names(rs)[rs==0], "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other                                     
#                  5088                  1529                  1185                                     

Pseudo_SE_psi <- SE_psi[SE_psi$loci %in% names(rs)[rs==0], ]
Real_SE_psi <- SE_psi[SE_psi$loci %in% names(rs)[rs>=1], ]

table(Real_SE_psi$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5272                  1859                  2684


#### A3SS =======================================================================================================================================


load("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi_new.RData")
colnames(A3SS_psi)[2:3] <- c("sj1", "sj2")
colnames(A3SS_psi)[7] <- "counts_start/end"

PSI <- A3SS_psi[, row.names(Cell_type)]
row.names(PSI) <- A3SS_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   0.000   0.000   2.791   1.000 575.000

sum(rs==0)
#[1] 2056
sum(rs==1)
#[1] 470
sum(rs>=1)
#[1] 987

Pseudo_A3SS_psi <- A3SS_psi[A3SS_psi$as2 %in% names(rs)[rs==0], ]
Real_A3SS_psi <- A3SS_psi[A3SS_psi$as2 %in% names(rs)[rs>=1], ]
dim(Pseudo_A3SS_psi)
# [1] 2056 3588
dim(Real_A3SS_psi)
# [1]  987 3588



#### A5SS =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_psi_new.RData")
colnames(A5SS_psi)[2:3] <- c("sj1", "sj2")
colnames(A5SS_psi)[7] <- "counts_start/end"

PSI <- A5SS_psi[, row.names(Cell_type)]
row.names(PSI) <- A5SS_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   0.000   0.000   3.632   1.000 960.000

sum(rs==0)
#[1] 2107
sum(rs==1)
#[1] 470
sum(rs>=1)
#[1] 923

Pseudo_A5SS_psi <- A5SS_psi[A5SS_psi$as2 %in% names(rs)[rs==0], ]
Real_A5SS_psi <- A5SS_psi[A5SS_psi$as2 %in% names(rs)[rs>=1], ]
dim(Pseudo_A5SS_psi)
# [1] 2107 3588
dim(Real_A5SS_psi)
# [1]  923 3588



table(Real_A3SS_psi$strand)
#   +   -   *
# 577 410   0
table(Real_A5SS_psi$strand)
#   +   -   *
# 392 531   0

table(Pseudo_A3SS_psi$strand)
#    +    -    *
# 1282  774    0
table(Pseudo_A5SS_psi$strand)
#   +    -    *
# 789 1318    0

## Real_A3SS_Parsed
577 + 531 = 1108
## Real_A5SS_Parsed
410 + 392 = 802

## Pseudo_A3SS_Parsed
1282 + 1318 = 2600
## Pseudo_A5SS_Parsed
774 + 789 = 1563

ASS_psi <- rbind(A3SS_psi,A5SS_psi)

Real_A3SS_PSI <- rbind(Real_A3SS_psi[Real_A3SS_psi$strand == "+", ], Real_A5SS_psi[Real_A5SS_psi$strand == "-", ])
Real_A5SS_PSI <- rbind(Real_A3SS_psi[Real_A3SS_psi$strand == "-", ], Real_A5SS_psi[Real_A5SS_psi$strand == "+", ])

Pseudo_A3SS_PSI <- rbind(Pseudo_A3SS_psi[Pseudo_A3SS_psi$strand == "+", ], Pseudo_A5SS_psi[Pseudo_A5SS_psi$strand == "-", ])
Pseudo_A5SS_PSI <- rbind(Pseudo_A3SS_psi[Pseudo_A3SS_psi$strand == "-", ], Pseudo_A5SS_psi[Pseudo_A5SS_psi$strand == "+", ])

A3SS_psi <- rbind(Real_A3SS_PSI, Pseudo_A3SS_PSI)
A5SS_psi <- rbind(Real_A5SS_PSI, Pseudo_A5SS_PSI)
dim(A3SS_psi)
#[1] 3708 3588
dim(A5SS_psi)
#[1] 2365 3588




#### AFE =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi_new.RData")

PSI <- AFE_psi[, row.names(Cell_type)]
row.names(PSI) <- AFE_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  	   0       0       1      16       5    2600

sum(rs==0)
#[1] 1093
sum(rs==1)
#[1] 537
sum(rs>=1)
#[1] 1780

Pseudo_AFE_psi <- AFE_psi[AFE_psi$as2 %in% names(rs)[rs==0], ]
dim(Pseudo_AFE_psi)
#[1] 1093 3585
Real_AFE_psi <- AFE_psi[AFE_psi$as2 %in% names(rs)[rs>=1], ]
dim(Real_AFE_psi)
#[1] 1780 3585




#### ALE =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_psi_new.RData")

PSI <- ALE_psi[, row.names(Cell_type)]
row.names(PSI) <- ALE_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.00    0.00    1.00   11.97    3.00 1646.00

sum(rs==0)
#[1] 644
sum(rs==1)
#[1] 245
sum(rs>=1)
#[1] 743

Pseudo_ALE_psi <- ALE_psi[ALE_psi$as2 %in% names(rs)[rs==0], ]
dim(Pseudo_ALE_psi)
#[1]  644 3585
Real_ALE_psi <- ALE_psi[ALE_psi$as2 %in% names(rs)[rs>=1], ]
dim(Real_ALE_psi)
#[1]  743 3585




#### MXE ========================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi_new.RData")

PSI <- MXE_psi[, row.names(Cell_type)]
row.names(PSI) <- paste(MXE_psi$as2, MXE_psi$as2.1, sep = "|")
rs <- rowSums(!is.na(PSI))
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.00    0.00    2.00   36.99   17.00  596.00
sum(rs>=1)
#[1] 75
PSI <- PSI[rs!=0, ]
rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   0.000   0.000   1.453   1.000  16.000
sum(rs==0)
#[1] 47
sum(rs==1)
#[1] 11
sum(rs>=1)
#[1] 28

Pseudo_MXE_psi <- MXE_psi[paste(MXE_psi$as2, MXE_psi$as2.1, sep = "|") %in% names(rs)[rs==0], ]
dim(Pseudo_MXE_psi)
#[1]   47 3594
Real_MXE_psi <- MXE_psi[paste(MXE_psi$as2, MXE_psi$as2.1, sep = "|") %in% names(rs)[rs>=1], ]
dim(Real_MXE_psi)
#[1]   28 3594

table(rowSums(is.na(Real_MXE_psi[,19:20])))
#  0  1  2
# 17  5  6
table(rowSums(is.na(Pseudo_MXE_psi[,19:20])))
#  0  1  2
# 14 18 15

test <- as.matrix(Real_MXE_psi[, row.names(Cell_type)])
row.names(test) <- paste(Real_MXE_psi$chr, ":", Real_MXE_psi$start, "-", Real_MXE_psi$end, sep = "")
test[is.na(test)] <- 0

cells <- vector()
for(i in 1:nrow(test)){
	cells[i] <- names(which.min(abs(test[i, ] - .5)))
}


sashimi <- function(junction, cell, left_expand = 1000, right_expand = 1000, tsv = "/mnt/data5/BGI/UCB/tangchao/data/BAM_tsv/", outDir = "~/", min = 10, ann_height = 4){
  if(dir.exists(outDir) == FALSE) dir.create(outDir)
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
    tmpdir <- tempdir()
    work1 <- paste("system('", "cat ", paste(tsv,cell,".tsv",sep = "", collapse = " "), " > ", tmpdir,"/tmp.tsv", "')", sep = "")
    work2 <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi2.sh ", window, paste(tmpdir,"/tmp.tsv", sep = ""), Dir_out, min, ann_height, "')")
    work3 <- paste("system('", "rm",paste(tmpdir,"/tmp.tsv", "')", sep = ""))
    eval(parse(text = work1))
    eval(parse(text = work2))
    eval(parse(text = work3))
  }
}

for(i in 1:nrow(test)){
	sashimi(junction = row.names(test)[i], cell = cells[i], outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/Real_MXE/", ann_height = 10)
}


test <- as.matrix(Pseudo_MXE_psi[, row.names(Cell_type)])
row.names(test) <- paste(Pseudo_MXE_psi$chr, ":", Pseudo_MXE_psi$start, "-", Pseudo_MXE_psi$end, sep = "")
test[is.na(test)] <- -1

cells <- list()
for(i in 1:nrow(test)){
	if(sum(0:1 %in% names(table(test[i, ])))==2){
		cells[[i]] <- c(sample(names(which(test[i, ] == 0)),1), sample(names(which(test[i, ] == 1)),1))
	}else{
		cells[[i]] <- NA
	}	
}

for(i in 1:nrow(test)){
	if(!is.na(cells[[i]])){
		sashimi(junction = row.names(test)[i], cell = cells[[i]], outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/Pseudo_MXE/", ann_height = 10)
	}
}




dim(Pseudo_SE_psi)
#[1] 7802 3586
dim(Pseudo_A3SS_psi)
#[1] 2056 3588
dim(Pseudo_A5SS_psi)
#[1] 2107 3588
dim(Pseudo_AFE_psi)
#[1] 1093 3585
dim(Pseudo_ALE_psi)
#[1]  644 3585
dim(Pseudo_MXE_psi)
#[1]   47 3594

dim(Real_SE_psi)
#[1] 9815 3586
dim(Real_A3SS_psi)
#[1]  987 3588
dim(Real_A5SS_psi)
#[1]  923 3588
dim(Real_AFE_psi)
#[1] 1780 3585
dim(Real_ALE_psi)
#[1]  743 3585
dim(Real_MXE_psi)
#[1]   28 3594

table(Real_SE_psi$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5272                  1859                  2684
table(Pseudo_SE_psi$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5088                  1529                  1185







#### Cell-type-specific AS events ===============================================================================================================
#### Cell-type-specific AS events ===============================================================================================================
#### Cell-type-specific AS events ===============================================================================================================
#### SE -----------------------------------------------------------------------------------------------------------------------------------------


load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_sep_SE.RData")

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))
# [1] 219  27  23 302 251 167 226

table(SE_psi[SE_psi$loci %in% sj_Bcell, "SE_type"])
# exon-skipping_exactly                 Other
#                   154                    65
table(SE_psi[SE_psi$loci %in% sj_Tcell, "SE_type"])
# exon-skipping_exactly                 Other
#                   158                    68
table(SE_psi[SE_psi$loci %in% sj_MK, "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other
#                   259                     2                    41
table(SE_psi[SE_psi$loci %in% sj_Mono, "SE_type"])
# exon-skipping_exactly                 Other
#                   197                    54
table(SE_psi[SE_psi$loci %in% sj_NK, "SE_type"])
# exon-skipping_exactly                 Other
#                  115                    52


table(Pseudo_SE_psi[Pseudo_SE_psi$loci %in% sj_Bcell, "SE_type"])
# exon-skipping_exactly                 Other
#                    22                     1
table(Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell, "SE_type"])
# exon-skipping_exactly                 Other
#                  132                    64


Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell & Real_SE_psi$SE_type == "Other", "sj13"]

dim(Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell & Real_SE_psi$SE_type == "Other", row.names(Cell_type)[Cell_type$CellType == "B Cells"]])


test <- as.matrix(Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell & Real_SE_psi$SE_type == "Other", row.names(Cell_type)[Cell_type$CellType == "B Cells"]])
row.names(test) <- Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell & Real_SE_psi$SE_type == "Other", "sj13"]
test[is.na(test)] <- 0

cells <- vector()
for(i in 1:nrow(test)){
	cells[i] <- names(which.min(abs(test[i, ] - .5)))
}

for(i in 1:nrow(test)){
	sashimi(junction = row.names(test)[i], 
			cell = cells[i], 
			outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell/Real_SE/Other/Sashimi/", 
			ann_height = 10)
}


test <- as.matrix(Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell & Real_SE_psi$SE_type == "exon-skipping_exactly", row.names(Cell_type)[Cell_type$CellType == "B Cells"]])
row.names(test) <- Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell & Real_SE_psi$SE_type == "exon-skipping_exactly", "sj13"]
test[is.na(test)] <- 0

cells <- vector()
for(i in 1:nrow(test)){
	cells[i] <- names(which.min(abs(test[i, ] - .5)))
}

for(i in 1:nrow(test)){
	sashimi(junction = row.names(test)[i], 
			cell = cells[i], 
			outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell/Real_SE/exon-skipping_exactly/Sashimi/", 
			ann_height = 10)
}





background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/DSU2/Document/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

SE_psi[SE_psi$loci %in% sj_Bcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"]


library(clusterProfiler)
library(org.Hs.eg.db)##### data base target to human-geneIDs

#### B cell

Bcel_SE_BP <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Bcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Bcel_SE_CC <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Bcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Bcel_SE_MF <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Bcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell/GO/all_Bcel_Sep_SE_GO.pdf",height = 10,width = 15)
print(dotplot(Bcel_SE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcel_SE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_SE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_SE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Bcel_SE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_SE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_SE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Bcel_SE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_SE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


#### T cell

Tcel_SE_BP <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Tcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Tcel_SE_CC <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Tcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Tcel_SE_MF <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Tcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Tcell/GO/all_Tcel_Sep_SE_GO.pdf",height = 10,width = 15)
print(dotplot(Tcel_SE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Tcel_SE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_SE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_SE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Tcel_SE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_SE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_SE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Tcel_SE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_SE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### MK cell

MK_SE_BP <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_MK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
MK_SE_CC <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_MK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
MK_SE_MF <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_MK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/MK/GO/all_MK_Sep_SE_GO.pdf",height = 10,width = 15)
print(dotplot(MK_SE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_SE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_SE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_SE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(MK_SE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_SE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_SE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(MK_SE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_SE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### Monocytes 

Mono_SE_BP <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Mono & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Mono_SE_CC <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Mono & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Mono_SE_MF <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_Mono & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Monocytes/GO/all_Mono_Sep_SE_GO.pdf",height = 10,width = 15)
print(dotplot(Mono_SE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_SE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_SE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_SE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Mono_SE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_SE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_SE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Mono_SE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_SE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### Monocytes 

NK_SE_BP <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_NK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
NK_SE_CC <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_NK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
NK_SE_MF <- enrichGO(gene          = SE_psi[SE_psi$loci %in% sj_NK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/NK/GO/all_NK_Sep_SE_GO.pdf",height = 10,width = 15)
print(dotplot(NK_SE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_SE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_SE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_SE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(NK_SE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_SE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_SE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(NK_SE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_SE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




all_gene <- list(Bcell = SE_psi[SE_psi$loci %in% sj_Bcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"], 
                 Tcell = SE_psi[SE_psi$loci %in% sj_Tcell & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"], 
                 MK = SE_psi[SE_psi$loci %in% sj_MK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                 Monocytes = SE_psi[SE_psi$loci %in% sj_Mono & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"],
                 NK = SE_psi[SE_psi$loci %in% sj_NK & SE_psi$SE_type == "exon-skipping_exactly", "gene_id"])

CompareGO_BP <- compareCluster(geneCluster   = all_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = all_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = all_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/GO/all_cell_Sep_SE_CompareGO.pdf", height = 16, width = 10)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "SE GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "SE GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "SE GO CC"))
dev.off()


CompareGO_BP2 <- simplify(CompareGO_BP, cutoff=0.7, by="p.adjust", select_fun=min)
CompareGO_CC2 <- simplify(CompareGO_CC, cutoff=0.7, by="p.adjust", select_fun=min)
CompareGO_MF2 <- simplify(CompareGO_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/GO/all_cell_Sep_SE_CompareGO_simplify.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP2,showCategory = 10,font.size = 16,title = "SE GO BP"))
print(dotplot(CompareGO_MF2,showCategory = 10,font.size = 16,title = "SE GO MF"))
print(dotplot(CompareGO_CC2,showCategory = 10,font.size = 16,title = "SE GO CC"))
dev.off()



CompareGO_BP <- compareCluster(geneCluster   = all_gene, 
                               #C      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = all_gene, 
                               #universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = all_gene, 
                               #universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/GO/all_cell_Sep_SE_CompareGO_universe.pdf", height = 16, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "SE GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "SE GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "SE GO CC"))
dev.off()





Bcel_SE_BP2 <- simplify(Bcel_SE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Bcel_SE_CC2 <- simplify(Bcel_SE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Bcel_SE_MF2 <- simplify(Bcel_SE_MF, cutoff=0.7, by="p.adjust", select_fun=min)


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell/GO/all_Bcel_Sep_SE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Bcel_SE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcel_SE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_SE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_SE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Bcel_SE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_SE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_SE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Bcel_SE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_SE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



#### T cell Sepcific T cell differentiation(GO:0030217) SE PSI and gene expression

Tcel_Sep_Tcel_differentiation_gene <- Tcel_SE_BP@result[Tcel_SE_BP@result$Description == "T cell differentiation","geneID"]
Tcel_Sep_Tcel_differentiation_gene <- unlist(strsplit(Tcel_Sep_Tcel_differentiation_gene,"/"))

dim(SE_psi[SE_psi$gene_id %in% Tcel_Sep_Tcel_differentiation_gene,])
#[1]   49 3586
test <- SE_psi[SE_psi$gene_id %in% Tcel_Sep_Tcel_differentiation_gene,]
dim(test[test$loci %in% sj_Tcell, ])
#[1]   34 3586
test <- test[test$loci %in% sj_Tcell, ]
row.names(test) <- paste(test$loci,"@",test$gene_id, sep = "")
test <- test[, row.names(Cell_type)]

test2 <- cbind(Cell_type, t(test))

TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/UCB_rsem_TPM_mat.txt", header = T, row.names = 1)
unique(do.call(rbind,strsplit(colnames(test2), "@"))[,2])[-(1:2)]

TPM <- TPM[unique(do.call(rbind,strsplit(colnames(test2), "@"))[,2])[-(1:2)], row.names(Cell_type)]
dim(TPM)
#[1]   13 3039
test3 <- cbind(Cell_type, t(TPM))

library(ggplot2)
#pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Tcell/GO/Tcel_Sep_Tcel_differentiation_SE_PSI_vioplot.pdf")
#ggplot(data = test2, aes(x = CellType, y = test2[,3]))+
#		geom_violin()+
#		ggtitle(colnames(test2)[3])+
#		labs(y = "PSI")#

#genei <- unique(do.call(rbind,strsplit(colnames(test2), "@"))[,2])[3]#

#ggplot(data = test3, aes(x = CellType, y = log10(1+test3[,genei])))+
#		geom_violin()+
#		ggtitle(genei)+
#		labs(y = "log10(TPM)")#

#dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Tcell/GO/Tcel_Sep_Tcel_differentiation_SE_PSI_vioplot_4.pdf")
for(i in 3:ncol(test2)){
  print(ggplot(data = test2, aes(x = CellType, y = test2[,i]))+
          geom_violin(adjust = 0.1, scale = "area")+
          geom_jitter(height = 0, size = .8, alpha = .5)+
          ggtitle(colnames(test2)[i])+
          labs(y = "PSI")+
          scale_x_discrete(breaks = names(table(test2[which(test2[,i]>0.15),1])),
                           labels = paste(names(table(test2[which(test2[,i]>0.15),1])),"\n" , round(table(test2[which(test2[,i]>0.15),1])/table(Cell_type$CellType),2), sep = "")))


genei <- do.call(rbind,strsplit(colnames(test2), "@"))[,2][i]
print(ggplot(data = test3, aes(x = CellType, y = log10(1+test3[,genei])))+
        geom_violin()+
        ggtitle(genei)+
        labs(y = "log10(TPM)")+
		scale_x_discrete(breaks = names(table(test3[which(test3[,genei]>=1),1])),
                labels = paste(names(table(test3[which(test3[,genei]>=1),1])),"\n" , round(table(test3[which(test3[,genei]>=1),1])/table(Cell_type$CellType),2), sep = "")))

}
dev.off()






#### B cell Sepcific B cell activation(GO:0042113) SE PSI and gene expression

Bcel_Sep_Bcel_activation_gene <- Bcel_SE_BP@result[Bcel_SE_BP@result$Description == "B cell activation","geneID"]
Bcel_Sep_Bcel_activation_gene <- unlist(strsplit(Bcel_Sep_Bcel_activation_gene,"/"))

dim(SE_psi[SE_psi$gene_id %in% Bcel_Sep_Bcel_activation_gene,])
#[1]   33 3586
test <- SE_psi[SE_psi$gene_id %in% Bcel_Sep_Bcel_activation_gene,]
dim(test[test$loci %in% sj_Bcell, ])
#[1]   27 3586
test <- test[test$loci %in% sj_Bcell, ]
row.names(test) <- paste(test$loci,"@",test$gene_id, sep = "")
test <- test[, row.names(Cell_type)]

test2 <- cbind(Cell_type, t(test))
dim(test2)
#[1] 3039   29

TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/UCB_rsem_TPM_mat.txt", header = T, row.names = 1)
unique(do.call(rbind,strsplit(colnames(test2), "@"))[,2])[-(1:2)]

TPM <- TPM[unique(do.call(rbind,strsplit(colnames(test2), "@"))[,2])[-(1:2)], row.names(Cell_type)]
dim(TPM)
#[1]   12 3039
test3 <- cbind(Cell_type, t(TPM))

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Bcell/GO/Bcel_Sep_Bcel_activation_SE_PSI_vioplot.pdf")
for(i in 3:ncol(test2)){
  print(ggplot(data = test2, aes(x = CellType, y = test2[,i]))+
          geom_violin()+
          geom_jitter(height = 0)+
          ggtitle(colnames(test2)[i])+
          labs(y = "PSI")+
          scale_x_discrete(breaks = names(table(test2[which(test2[,i]>0.15),1])),
                           labels = paste(names(table(test2[which(test2[,i]>0.15),1])),"\n" , round(table(test2[which(test2[,i]>0.15),1])/table(Cell_type$CellType),2), sep = "")))


genei <- do.call(rbind,strsplit(colnames(test2), "@"))[,2][i]
print(ggplot(data = test3, aes(x = CellType, y = log10(1+test3[,genei])))+
        geom_violin()+
        ggtitle(genei)+
        labs(y = "log10(TPM)"))
}
dev.off()



#### Pseudo vs. Real

load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_sep_SE.RData")

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))
# [1] 385 237 434 309 260 234 444

dim(Real_SE_psi[Real_SE_psi$loci %in% sj_Bcell, row.names(Cell_type)[Cell_type$CellType == "B Cells"]])
#[1] 357 607
dim(Pseudo_SE_psi[Pseudo_SE_psi$loci %in% sj_Bcell, row.names(Cell_type)[Cell_type$CellType == "B Cells"]])
#[1]  28 607

dim(Real_SE_psi[Real_SE_psi$loci %in% sj_Tcell, row.names(Cell_type)[Cell_type$CellType == "B Cells"]])
#[1] 435 607
dim(Pseudo_SE_psi[Pseudo_SE_psi$loci %in% sj_Tcell, row.names(Cell_type)[Cell_type$CellType == "B Cells"]])
#[1] 9 607


row.names(Real_SE_psi) <- Real_SE_psi$loci
row.names(Pseudo_SE_psi) <- Pseudo_SE_psi$loci

Real_SE_Bcel <- Real_SE_psi[, row.names(Cell_type)[Cell_type$CellType == "B Cells"]]
Pseudo_SE_Bcel <- Pseudo_SE_psi[, row.names(Cell_type)[Cell_type$CellType == "B Cells"]]

dim(Real_SE_Bcel)
#[1] 9815  607
sum(rowSums(!is.na(Real_SE_Bcel))>1)
#[1] 9127

dim(Pseudo_SE_Bcel)
#[1] 7802  607
sum(rowSums(!is.na(Pseudo_SE_Bcel))>1)
#[1] 5192

Real_SE_psi[Real_SE_psi$loci %in% row.names(Real_SE_Bcel)[rowSums(!is.na(Real_SE_Bcel))>1],"gene_id"]
unique(Real_SE_psi[Real_SE_psi$loci %in% row.names(Real_SE_Bcel)[rowSums(!is.na(Real_SE_Bcel))>1],"gene_id"])
Real_SE_Bcel_gene <- unique(Real_SE_psi[Real_SE_psi$loci %in% row.names(Real_SE_Bcel)[rowSums(!is.na(Real_SE_Bcel))>1],"gene_id"])[-2]

unique(Pseudo_SE_psi[Pseudo_SE_psi$loci %in% row.names(Pseudo_SE_Bcel)[rowSums(!is.na(Pseudo_SE_Bcel))>1],"gene_id"])
Pseudo_SE_Bcel_gene <- unique(Pseudo_SE_psi[Pseudo_SE_psi$loci %in% row.names(Pseudo_SE_Bcel)[rowSums(!is.na(Pseudo_SE_Bcel))>1],"gene_id"])[-1]


all_gene = list(Real = Real_SE_Bcel_gene, Pseudo = Pseudo_SE_Bcel_gene)
CompareGO_BP <- compareCluster(geneCluster   = all_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = all_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = all_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)


CompareGO_BP2 <- simplify(CompareGO_BP, cutoff=0.7, by="p.adjust", select_fun=min)
CompareGO_CC2 <- simplify(CompareGO_CC, cutoff=0.7, by="p.adjust", select_fun=min)
CompareGO_MF2 <- simplify(CompareGO_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/GO/Bcell_SE_Real_Pseudo_CompareGO_simplify.pdf", height = 16, width = 15)
print(dotplot(CompareGO_BP2,showCategory = 10,font.size = 16,title = "SE GO BP"))
print(dotplot(CompareGO_MF2,showCategory = 10,font.size = 16,title = "SE GO MF"))
print(dotplot(CompareGO_CC2,showCategory = 10,font.size = 16,title = "SE GO CC"))
dev.off()














#################################################################################################################################################
#### A3SS & A5SS ================================================================================================================================
#################################################################################################################################################





load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Cell_Sep/All_Cell_sep_ASS.RData")

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))
# [1]  62   9   7 127  68  61  32

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))
# [1]  89  46  68 136  75  81  68


dim(A3SS_psi[A3SS_psi$as2 %in% sj_Bcell, ])
#[1]   66 3588
dim(A5SS_psi[A5SS_psi$as2 %in% sj_Bcell, ])
#[1]   23 3588

test3 <- A3SS_psi[A3SS_psi$as2 %in% sj_Bcell, row.names(Cell_type)]
row.names(test3) <- test3$as2
test5 <- A5SS_psi[A5SS_psi$as2 %in% sj_Bcell, row.names(Cell_type)]
row.names(test5) <- test5$as2
test <- rbind(test3, test5)



background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/DSU2/Document/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

union(A3SS_psi[A3SS_psi$as2 %in% sj_Bcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Bcell, "gene_id"])

library(clusterProfiler)
library(org.Hs.eg.db)##### data base target to human-geneIDs

#### B cell

Bcel_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Bcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Bcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Bcel_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Bcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Bcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Bcel_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Bcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Bcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/Bcell/GO/all_Bcel_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(Bcel_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcel_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Bcel_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Bcel_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


Bcel_ASS_BP2 <- simplify(Bcel_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Bcel_ASS_CC2 <- simplify(Bcel_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Bcel_ASS_MF2 <- simplify(Bcel_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/Bcell/GO/all_Bcel_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Bcel_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcel_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Bcel_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Bcel_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



#### T cell



Tcel_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Tcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Tcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Tcel_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Tcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Tcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Tcel_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Tcell, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Tcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/Tcell/GO/all_Tcel_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(Tcel_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Tcel_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Tcel_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Tcel_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


Tcel_ASS_BP2 <- simplify(Tcel_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Tcel_ASS_CC2 <- simplify(Tcel_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Tcel_ASS_MF2 <- simplify(Tcel_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/Tcell/GO/all_Tcel_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Tcel_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Tcel_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Tcel_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Tcel_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### MK cell



MK_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_MK, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_MK, "gene_id"]),
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
MK_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_MK, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_MK, "gene_id"]),
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
MK_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_MK, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_MK, "gene_id"]),
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/MK/GO/all_MK_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(MK_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(MK_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(MK_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


MK_ASS_BP2 <- simplify(MK_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
MK_ASS_CC2 <- simplify(MK_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
MK_ASS_MF2 <- simplify(MK_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/MK/GO/all_MK_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(MK_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(MK_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(MK_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()





#### Mono cell



Mono_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Mono, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Mono, "gene_id"]),
                      universe      = background,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
Mono_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Mono, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Mono, "gene_id"]),
                      universe      = background,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
Mono_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_Mono, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_Mono, "gene_id"]),
                      universe      = background,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/Mono/GO/all_Mono_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(Mono_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Mono_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Mono_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


Mono_ASS_BP2 <- simplify(Mono_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Mono_ASS_CC2 <- simplify(Mono_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Mono_ASS_MF2 <- simplify(Mono_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/Mono/GO/all_Mono_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Mono_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Mono_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Mono_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()







#### NK cell



NK_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_NK, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_NK, "gene_id"]),
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
NK_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_NK, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_NK, "gene_id"]),
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
NK_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_NK, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_NK, "gene_id"]),
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/NK/GO/all_NK_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(NK_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(NK_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(NK_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


NK_ASS_BP2 <- simplify(NK_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
NK_ASS_CC2 <- simplify(NK_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
NK_ASS_MF2 <- simplify(NK_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/NK/GO/all_NK_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(NK_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(NK_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(NK_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



#### CD4 cell



CD4_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_CD4, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_CD4, "gene_id"]),
                      universe      = background,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
CD4_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_CD4, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_CD4, "gene_id"]),
                      universe      = background,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
CD4_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_CD4, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_CD4, "gene_id"]),
                      universe      = background,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/CD4/GO/all_CD4_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(CD4_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD4_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD4_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD4_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


CD4_ASS_BP2 <- simplify(CD4_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
CD4_ASS_CC2 <- simplify(CD4_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
CD4_ASS_MF2 <- simplify(CD4_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/CD4/GO/all_CD4_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(CD4_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD4_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD4_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD4_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()







#### CD8 cell



CD8_ASS_BP <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_CD8, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_CD8, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CD8_ASS_CC <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_CD8, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_CD8, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CD8_ASS_MF <- enrichGO(gene          = union(A3SS_psi[A3SS_psi$as2 %in% sj_CD8, "gene_id"],A5SS_psi[A5SS_psi$as2 %in% sj_CD8, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/CD8/GO/all_CD8_Sep_ASS_GO.pdf",height = 10,width = 15)
print(dotplot(CD8_ASS_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD8_ASS_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_ASS_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_ASS_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD8_ASS_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_ASS_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_ASS_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD8_ASS_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_ASS_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


CD8_ASS_BP2 <- simplify(CD8_ASS_BP, cutoff=0.7, by="p.adjust", select_fun=min)
CD8_ASS_CC2 <- simplify(CD8_ASS_CC, cutoff=0.7, by="p.adjust", select_fun=min)
CD8_ASS_MF2 <- simplify(CD8_ASS_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/ASS/Cell_Sep/CD8/GO/all_CD8_Sep_ASS_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(CD8_ASS_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD8_ASS_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_ASS_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_ASS_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD8_ASS_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_ASS_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_ASS_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD8_ASS_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_ASS_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()






#################################################################################################################################################
#### AFE & ALE ################################################################################################################################
#################################################################################################################################################


load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Cell_Sep/All_Cell_sep_AFLE.RData")

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))
# [1]  37  2  8 82 74 38 40

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))
# [1] 59 37 95 83 76 53 91



dim(AFE_psi[AFE_psi$as2 %in% sj_Bcell, ])
#[1]   34 3588
dim(ALE_psi[ALE_psi$as2 %in% sj_Bcell, ])
#[1]   25 3588


background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/DSU2/Document/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

union(AFE_psi[AFE_psi$as2 %in% sj_Bcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Bcell, "gene_id"])

library(clusterProfiler)
library(org.Hs.eg.db)##### data base target to human-geneIDs

#### B cell

Bcel_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Bcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Bcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Bcel_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Bcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Bcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Bcel_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Bcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Bcell, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/Bcell/GO/all_Bcel_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(Bcel_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcel_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Bcel_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Bcel_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


Bcel_AFLE_BP2 <- simplify(Bcel_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Bcel_AFLE_CC2 <- simplify(Bcel_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Bcel_AFLE_MF2 <- simplify(Bcel_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/Bcell/GO/all_Bcel_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Bcel_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcel_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Bcel_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Bcel_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Bcel_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcel_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



#### T cell


Tcel_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Tcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Tcell, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
Tcel_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Tcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Tcell, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
Tcel_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Tcell, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Tcell, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/Tcell/GO/all_Tcel_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(Tcel_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Tcel_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Tcel_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Tcel_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


Tcel_AFLE_BP2 <- simplify(Tcel_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Tcel_AFLE_CC2 <- simplify(Tcel_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Tcel_AFLE_MF2 <- simplify(Tcel_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/Tcell/GO/all_Tcel_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Tcel_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Tcel_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Tcel_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Tcel_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Tcel_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Tcel_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



#### MK cell


MK_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_MK, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_MK, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
MK_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_MK, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_MK, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
MK_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_MK, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_MK, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/MK/GO/all_MK_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(MK_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(MK_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(MK_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


MK_AFLE_BP2 <- simplify(MK_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
MK_AFLE_CC2 <- simplify(MK_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
MK_AFLE_MF2 <- simplify(MK_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/MK/GO/all_MK_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(MK_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(MK_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(MK_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(MK_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### Mono cell


Mono_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Mono, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Mono, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Mono_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Mono, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Mono, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
Mono_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_Mono, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_Mono, "gene_id"]),
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/Mono/GO/all_Mono_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(Mono_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Mono_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Mono_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


Mono_AFLE_BP2 <- simplify(Mono_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
Mono_AFLE_CC2 <- simplify(Mono_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
Mono_AFLE_MF2 <- simplify(Mono_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/Mono/GO/all_Mono_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(Mono_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(Mono_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(Mono_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(Mono_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### NK cell


NK_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_NK, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_NK, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
NK_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_NK, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_NK, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
NK_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_NK, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_NK, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/NK/GO/all_NK_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(NK_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(NK_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(NK_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


NK_AFLE_BP2 <- simplify(NK_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
NK_AFLE_CC2 <- simplify(NK_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
NK_AFLE_MF2 <- simplify(NK_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/NK/GO/all_NK_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(NK_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(NK_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(NK_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(NK_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### CD4 cell


CD4_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_CD4, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_CD4, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
CD4_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_CD4, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_CD4, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
CD4_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_CD4, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_CD4, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/CD4/GO/all_CD4_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(CD4_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD4_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD4_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD4_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


CD4_AFLE_BP2 <- simplify(CD4_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
CD4_AFLE_CC2 <- simplify(CD4_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
CD4_AFLE_MF2 <- simplify(CD4_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/CD4/GO/all_CD4_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(CD4_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD4_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD4_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD4_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD4_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



#### CD8 cell


CD8_AFLE_BP <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_CD8, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_CD8, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
CD8_AFLE_CC <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_CD8, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_CD8, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
CD8_AFLE_MF <- enrichGO(gene          = union(AFE_psi[AFE_psi$as2 %in% sj_CD8, "gene_id"],ALE_psi[ALE_psi$as2 %in% sj_CD8, "gene_id"]),
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/CD8/GO/all_CD8_Sep_AFLE_GO.pdf",height = 10,width = 15)
print(dotplot(CD8_AFLE_BP,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD8_AFLE_BP,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_AFLE_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_AFLE_CC,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD8_AFLE_CC,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_AFLE_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_AFLE_MF,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD8_AFLE_MF,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_AFLE_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


CD8_AFLE_BP2 <- simplify(CD8_AFLE_BP, cutoff=0.7, by="p.adjust", select_fun=min)
CD8_AFLE_CC2 <- simplify(CD8_AFLE_CC, cutoff=0.7, by="p.adjust", select_fun=min)
CD8_AFLE_MF2 <- simplify(CD8_AFLE_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFLE/Cell_Sep/CD8/GO/all_CD8_Sep_AFLE_GO_simplify.pdf",height = 10,width = 15)
print(dotplot(CD8_AFLE_BP2,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD8_AFLE_BP2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_AFLE_BP2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_AFLE_CC2,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(CD8_AFLE_CC2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_AFLE_CC2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(CD8_AFLE_MF2,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(CD8_AFLE_MF2,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_AFLE_MF2,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()









#### AS type GO =================================================================================================================================



SE_real_gene <- unique(Real_SE_psi[!is.na(Real_SE_psi$gene_id), "gene_id"])
AFE_real_gene <- unique(Real_AFE_psi[!is.na(Real_AFE_psi$gene_id), "gene_id"])
ALE_real_gene <- unique(Real_ALE_psi[!is.na(Real_ALE_psi$gene_id), "gene_id"])

A3SS_real_gene <- unique(Real_A3SS_PSI$gene_id)
A5SS_real_gene <- unique(Real_A5SS_PSI$gene_id)




unlist(lapply(list(SE_real_gene, A3SS_real_gene, A5SS_real_gene, AFE_real_gene, ALE_real_gene),length))
# [1] 3044  933  714 1258  627

AS_gene <- list(SE = SE_real_gene, A3SS = A3SS_real_gene, A5SS = A5SS_real_gene, AFE = AFE_real_gene, ALE = ALE_real_gene)

CompareGO_BP <- compareCluster(geneCluster   = AS_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_AS_type_CompareGO.pdf", height = 16, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()

CompareGO_BP2 <- simplify(CompareGO_BP, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_AS_type_CompareGO_simplify.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP2,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()

CompareGO_BP_l4 <- gofilter(CompareGO_BP, level = 4)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_AS_type_CompareGO_L4.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP_l4,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()

CompareGO_BP_l3 <- gofilter(CompareGO_BP, level = 3)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_AS_type_CompareGO_L3.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP_l3,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()




SE_real_BP <- enrichGO(gene          = SE_real_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
A3SS_real_BP <- enrichGO(gene          = A3SS_real_gene,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A5SS_real_BP <- enrichGO(gene          = A5SS_real_gene,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
AFE_real_BP <- enrichGO(gene          = AFE_real_gene,
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
ALE_real_BP <- enrichGO(gene          = ALE_real_gene,
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

SE_real_BP@result -> SE
A3SS_real_BP@result -> A3SS
A5SS_real_BP@result -> A5SS
AFE_real_BP@result -> AFE
ALE_real_BP@result -> ALE


myGOslim <- function(x){
  requireNamespace("topGO") || stop("package topGO is required")
  groupGOTerms <- rvcheck::get_fun_from_pkg("topGO", "groupGOTerms")
  annFUN.gene2GO <- rvcheck::get_fun_from_pkg("topGO", "annFUN.gene2GO")
  if (!class(x) %in% c("gseaResult", "enrichResult")) {
    stop("x should be output of gseGO or enrichGO...")
  }
  gene2GO <- inverseList(x@geneSets)
  if (is(x, "gseaResult")) {
    ont <- x@setType
    allgenes <- x@geneList
    core_genes <- unique(unlist(geneInCategory(x)))
    allgenes[!names(allgenes) %in% core_genes] <- -1
    allgenes[core_genes] <- 1
  }
  else {
    ont <- x@ontology
    universe <- x@universe
    allgenes <- numeric(length(universe))
    names(allgenes) <- universe
    allgenes[x@gene] <- 1
  }
  selector <- function(scores) return(scores == 1)
  if (!ont %in% c("BP", "MF", "CC")) {
    stop("ontology should be one of 'BP', 'MF' or 'CC'...")
  }
  pvalue <- x@result$p.adjust
  names(pvalue) <- x@result$ID
  groupGOTerms()
  GOdata <- new("topGOdata", description = "clusterProfiler enrichment results", 
                ontology = ont, 
                allGenes = allgenes, 
                geneSel = selector, 
                annot = annFUN.gene2GO, 
                gene2GO = gene2GO)
  result <- buildLevels(graph(GOdata))
  return(result)
}

SE_DAG <- myGOslim(SE_real_BP)
A3SS_DAG <- myGOslim(A3SS_real_BP)
A5SS_DAG <- myGOslim(A5SS_real_BP)
AFE_DAG <- myGOslim(AFE_real_BP)
ALE_DAG <- myGOslim(ALE_real_BP)


getNoOfLevels(SE_DAG)
getNoOfLevels(A3SS_DAG)
getNoOfLevels(A5SS_DAG)
getNoOfLevels(AFE_DAG)
getNoOfLevels(ALE_DAG)

SE_DAG$level2nodes$`5`
A3SS_DAG$level2nodes$`5`
A5SS_DAG$level2nodes$`5`
AFE_DAG$level2nodes$`5`
ALE_DAG$level2nodes$`5`

lapply(FUN = length, list(SE_DAG$level2nodes$`2`,
                          A3SS_DAG$level2nodes$`2`,
                          A5SS_DAG$level2nodes$`2`,
                          AFE_DAG$level2nodes$`2`,
                          ALE_DAG$level2nodes$`2`))

length(Reduce(intersect, list(SE_DAG$level2nodes$`2`,
                              A3SS_DAG$level2nodes$`2`,
                              A5SS_DAG$level2nodes$`2`,
                              AFE_DAG$level2nodes$`2`,
                              ALE_DAG$level2nodes$`2`)))


select(GO.db, keys = SE_DAG$level2nodes$`2`, columns = "TERM")

SE$Description[1]
Term(GOTERM[[SE$ID[1]]])

select(GO.db, keys = GOBPANCESTOR[[SE$ID[1]]][GOBPANCESTOR[[SE$ID[1]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")
select(GO.db, keys = GOBPANCESTOR[[SE$ID[2]]][GOBPANCESTOR[[SE$ID[2]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")
select(GO.db, keys = GOBPANCESTOR[[SE$ID[12]]][GOBPANCESTOR[[SE$ID[12]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")




#### enrichment at level 3

# SE[SE$ID %in% SE_DAG$level2nodes$`3`,]$Description
# # [1] "multi-organism metabolic process"    "antigen processing and presentation" "protein folding"                    
# # [4] "methylation"
# A3SS[A3SS$ID %in% A3SS_DAG$level2nodes$`3`,]$Description
# # [1] "multi-organism metabolic process"    "antigen processing and presentation"
# A5SS[A5SS$ID %in% A5SS_DAG$level2nodes$`3`,]$Description
# # [1] "multi-organism metabolic process"    "antigen processing and presentation" "multi-organism localization"  
# AFE[AFE$ID %in% AFE_DAG$level2nodes$`3`,]$Description
# # [1] "multi-organism metabolic process" "protein folding" 
# ALE[ALE$ID %in% ALE_DAG$level2nodes$`3`,]$Description
# # [1] "multi-organism metabolic process"    "antigen processing and presentation"

SE[SE$ID %in% SE_DAG$level2nodes$`3`,]$Description
# [1] "antigen processing and presentation" "protein folding"
# [3] "leukocyte proliferation"             "methylation"
# [5] "multi-organism localization"         "chromosome segregation"
A3SS[A3SS$ID %in% A3SS_DAG$level2nodes$`3`,]$Description
# [1] "multi-organism localization" "leukocyte proliferation"
A5SS[A5SS$ID %in% A5SS_DAG$level2nodes$`4`,]$Description
# [1] "ribonucleoprotein complex biogenesis"
# [2] "adaptive immune response"
AFE[AFE$ID %in% AFE_DAG$level2nodes$`3`,]$Description
# [1] "process utilizing autophagic mechanism"
# [2] "protein folding"
# [3] "multi-organism localization"
# [4] "leukocyte proliferation"
# [5] "antigen processing and presentation"
ALE[ALE$ID %in% ALE_DAG$level2nodes$`4`,]$Description
# [1] "ribonucleoprotein complex biogenesis"
# [2] "generation of precursor metabolites and energy"
# [3] "antigen processing and presentation of exogenous antigen"



SE_L3 <-SE[SE$ID %in% SE_DAG$level2nodes$`3`,]
A3SS_L3 <- A3SS[A3SS$ID %in% A3SS_DAG$level2nodes$`3`,]
A5SS_L3 <- A5SS[A5SS$ID %in% A5SS_DAG$level2nodes$`4`,]
AFE_L3 <- AFE[AFE$ID %in% AFE_DAG$level2nodes$`3`,]
ALE_L3 <- ALE[ALE$ID %in% ALE_DAG$level2nodes$`4`,]

#ALE_L3$Description[3] <- AFE_L3$Description[5]
#ALE_L3$ID[3] <- AFE_L3$ID[5]

length(unique(c(SE_L3$ID, A3SS_L3$ID, A5SS_L3$ID, AFE_L3$ID, ALE_L3$ID)))
mat_l3 <- matrix(0, nrow = 5, ncol = length(unique(c(SE_L3$ID, A3SS_L3$ID, A5SS_L3$ID, AFE_L3$ID, ALE_L3$ID))))
row.names(mat_l3) <- c("SE","A3SS","A5SS","AFE","ALE")
colnames(mat_l3) <- unique(c(SE_L3$Description, A3SS_L3$Description, A5SS_L3$Description, AFE_L3$Description, ALE_L3$Description))

for(i in 1:nrow(SE_L3)){
  mat_l3[1, SE_L3[i,"Description"]] <- SE_L3[i,"Count"]
}

for(i in 1:nrow(A3SS_L3)){
  mat_l3[2, A3SS_L3[i,"Description"]] <- A3SS_L3[i,"Count"]
}

for(i in 1:nrow(A5SS_L3)){
  mat_l3[3, A5SS_L3[i,"Description"]] <- A5SS_L3[i,"Count"]
}

for(i in 1:nrow(AFE_L3)){
  mat_l3[4, AFE_L3[i,"Description"]] <- AFE_L3[i,"Count"]
}

for(i in 1:nrow(ALE_L3)){
  mat_l3[5, ALE_L3[i,"Description"]] <- ALE_L3[i,"Count"]
}


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/AS_GO_Count_circos_at_levels3.pdf")
#circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2), start.degree = 230, gap.after = c(rep(2, nrow(mat_l3)-1), 10, rep(2, ncol(mat_l3)-1), 10))#
par(mar = c(6,5,2,1))
circos.par(start.degree = 240, gap.after = c(rep(2, nrow(mat_l3)-1), 60, rep(2, ncol(mat_l3)-1), 60))# 块儿与块儿之间的空格
chordDiagram(mat_l3, annotationTrack = "grid", transparency = 0.5, 
             grid.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF",colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat_l3))),
             row.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF"),
             preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% get.all.sector.index()[1:5]){
    circos.text(mean(xlim), 0.3*nchar(sector.name), sector.name, niceFacing = T, facing = "clockwise")
  }
}, bg.border = NA)

circos.clear()

library(ComplexHeatmap)
lgd_points = Legend(at = colnames(mat_l3), type = "points", 
                    legend_gp = gpar(col = colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat_l3))), 
                    background = colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat_l3)),
                    title_position = "topleft", 
                    title = "GO Term", 
                    nrow = ncol(mat_l3))

lgd_list_vertical = packLegend(lgd_points, direction = "vertical")

pushViewport(viewport(x = unit(100, "mm"), y = unit(28, "mm"), 
                      width = grobWidth(lgd_list_vertical), 
                      height = grobHeight(lgd_list_vertical)))
grid.draw(lgd_list_vertical)
upViewport()
dev.off()


#### Pseudo -------------------------------------------------------------------------------------------------------------------------------------


SE_pseudo_gene <- unique(Pseudo_SE_psi[!is.na(Pseudo_SE_psi$gene_id), "gene_id"])
AFE_pseudo_gene <- unique(Pseudo_AFE_psi[!is.na(Pseudo_AFE_psi$gene_id), "gene_id"])
ALE_pseudo_gene <- unique(Pseudo_ALE_psi[!is.na(Pseudo_ALE_psi$gene_id), "gene_id"])

A3SS_pseudo_gene <- unique(Pseudo_A3SS_PSI$gene_id)
A5SS_pseudo_gene <- unique(Pseudo_A5SS_PSI$gene_id)




unlist(lapply(list(SE_pseudo_gene, A3SS_pseudo_gene, A5SS_pseudo_gene, AFE_pseudo_gene, ALE_pseudo_gene),length))
#[1] 3685 2142 1404 1055  622

AS_gene <- list(SE = SE_pseudo_gene, A3SS = A3SS_pseudo_gene, A5SS = A5SS_pseudo_gene, AFE = AFE_pseudo_gene, ALE = ALE_pseudo_gene)

CompareGO_BP <- compareCluster(geneCluster   = AS_gene, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_Pseudo_AS_type_CompareGO.pdf", height = 16, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()

CompareGO_BP2 <- simplify(CompareGO_BP, cutoff=0.7, by="p.adjust", select_fun=min)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_Pseudo_AS_type_CompareGO_simplify.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP2,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()

CompareGO_BP_l4 <- gofilter(CompareGO_BP, level = 4)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_Pseudo_AS_type_CompareGO_L4.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP_l4,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()

CompareGO_BP_l3 <- gofilter(CompareGO_BP, level = 3)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/GO/all_Pseudo_AS_type_CompareGO_L3.pdf", height = 10, width = 15)
print(dotplot(CompareGO_BP_l3,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()




SE_pseudo_BP <- enrichGO(gene          = SE_pseudo_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
A3SS_pseudo_BP <- enrichGO(gene          = A3SS_pseudo_gene,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A5SS_pseudo_BP <- enrichGO(gene          = A5SS_pseudo_gene,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
AFE_pseudo_BP <- enrichGO(gene          = AFE_pseudo_gene,
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
ALE_pseudo_BP <- enrichGO(gene          = ALE_pseudo_gene,
                        universe      = background,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

SE_pseudo_BP@result -> SE
A3SS_pseudo_BP@result -> A3SS
A5SS_pseudo_BP@result -> A5SS
AFE_pseudo_BP@result -> AFE
ALE_pseudo_BP@result -> ALE


SE_DAG <- myGOslim(SE_pseudo_BP)
A3SS_DAG <- myGOslim(A3SS_pseudo_BP)
A5SS_DAG <- myGOslim(A5SS_pseudo_BP)
AFE_DAG <- myGOslim(AFE_pseudo_BP)
ALE_DAG <- myGOslim(ALE_pseudo_BP)


getNoOfLevels(SE_DAG)
getNoOfLevels(A3SS_DAG)
getNoOfLevels(A5SS_DAG)
getNoOfLevels(AFE_DAG)
getNoOfLevels(ALE_DAG)

SE_DAG$level2nodes$`5`
A3SS_DAG$level2nodes$`5`
A5SS_DAG$level2nodes$`5`
AFE_DAG$level2nodes$`5`
ALE_DAG$level2nodes$`5`

lapply(FUN = length, list(SE_DAG$level2nodes$`2`,
                          A3SS_DAG$level2nodes$`2`,
                          A5SS_DAG$level2nodes$`2`,
                          AFE_DAG$level2nodes$`2`,
                          ALE_DAG$level2nodes$`2`))

length(Reduce(intersect, list(SE_DAG$level2nodes$`2`,
                              A3SS_DAG$level2nodes$`2`,
                              A5SS_DAG$level2nodes$`2`,
                              AFE_DAG$level2nodes$`2`,
                              ALE_DAG$level2nodes$`2`)))


select(GO.db, keys = SE_DAG$level2nodes$`2`, columns = "TERM")

SE$Description[1]
Term(GOTERM[[SE$ID[1]]])

select(GO.db, keys = GOBPANCESTOR[[SE$ID[1]]][GOBPANCESTOR[[SE$ID[1]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")
select(GO.db, keys = GOBPANCESTOR[[SE$ID[2]]][GOBPANCESTOR[[SE$ID[2]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")
select(GO.db, keys = GOBPANCESTOR[[SE$ID[12]]][GOBPANCESTOR[[SE$ID[12]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")


SE[SE$ID %in% SE_DAG$level2nodes$`3`,]$Description
# [1] "chromosome segregation"
# [2] "process utilizing autophagic mechanism"
# [3] "membrane docking"
# [4] "methylation"
A3SS[A3SS$ID %in% A3SS_DAG$level2nodes$`3`,]$Description
# [1] "chromosome segregation" "membrane docking"
A5SS[A5SS$ID %in% A5SS_DAG$level2nodes$`4`,]$Description
# [1] "ribonucleoprotein complex biogenesis"
AFE[AFE$ID %in% AFE_DAG$level2nodes$`3`,]$Description
ALE[ALE$ID %in% ALE_DAG$level2nodes$`4`,]$Description
# [1] "ribonucleoprotein complex biogenesis"
# [2] "generation of precursor metabolites and energy"
# [3] "antigen processing and presentation of exogenous antigen"


#### Novel SE -----------------------------------------------------------------------------------------------------------------------------------


dim(Real_SE_psi[Real_SE_psi$SE_type == "novel_exon-skipping",])
#[1] 1859 3586
test <- Real_SE_psi[Real_SE_psi$SE_type == "novel_exon-skipping",]
test <- test[,grep("UCB", colnames(test))]
apply(test,1,function(x){table(substr(colnames(test)[!is.na(x)],1,4))})
table(unlist(lapply(apply(test,1,function(x){table(substr(colnames(test)[!is.na(x)],1,4))}), length)))
#   1    2    3
#  48   86 1725
table(unlist(lapply(apply(test,1,function(x){table(substr(colnames(test)[!is.na(x)],1,4))}), length)))/nrow(test)*100
#        1         2         3
# 2.582033  4.626143 92.791824

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Novel_SE_individual_specific.pdf", height = 10, width = 15)
x = barplot(table(unlist(lapply(apply(test,1,function(x){table(substr(colnames(test)[!is.na(x)],1,4))}), length)))/nrow(test)*100, ylim = c(0,100))
text(x,table(unlist(lapply(apply(test,1,function(x){table(substr(colnames(test)[!is.na(x)],1,4))}), length)))/nrow(test)*100, labels=table(unlist(lapply(apply(test,1,function(x){table(substr(colnames(test)[!is.na(x)],1,4))}), length))),pos=3)
dev.off()












