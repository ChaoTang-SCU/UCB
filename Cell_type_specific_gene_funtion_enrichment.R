load("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi_new.RData")
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

dim(MXE_psi[,row.names(Cell_type)])
[1]  116 3039
MXE_psi[,row.names(Cell_type)] -> tab
tab[is.na(tab)] <- -.2
library(pheatmap)
library(RColorBrewer)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/All_MXE_heatmap.pdf")
pheatmap(tab[rowSums(tab!=-0.2)>0, ], annotation_col = Cell_type, show_rownames = F, show_colnames = F, cluster_cols = F,
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 4)(100))
dev.off()

pheatmap(tab, annotation_col = Cell_type, show_rownames = T, show_colnames = F, cluster_cols = F,
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 4)(100))

tab[c("8283315", "8283316", "48972"),1:20]
MXE_psi[c("8283315", "8283316", "48972"),1:30]
# gene_id1: ENSG00000163534
# gene: FCRL1
# function: May function as an activating coreceptor in B-cells. May function in B-cells activation and differentiation.
# URL: https://www.uniprot.org/uniprot/Q96LA6

# gene_id1: ENSG00000175857
# gene: GAPT
# function: Negatively regulates B-cell proliferation following stimulation through the B-cell receptor. 
# 			May play an important role in maintenance of marginal zone (MZ) B-cells (By similarity).
# URL: https://www.uniprot.org/uniprot/Q8N292

# UCB3.00779_5/58491542-58494246.pdf
# UCB4.01017_1/157797133-157798160.pdf


MXE_psi[c("8283315", "8283316"), row.names(Cell_type)] -> test
sum(colSums(!is.na(test))>0)
[1] 219
test[,colSums(!is.na(test))>0] -> test
pheatmap(test, annotation_col = Cell_type[colnames(test),], show_rownames = F, show_colnames = F, cluster_cols = F,
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 4)(100))
pheatmap(test, annotation_col = Cell_type[colnames(test),], show_rownames = F, show_colnames = F, cluster_cols = F)




sashimi(junction = "1:157797133-157798160", cell = "UCB3.00501", outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/Real_MXE/")
sashimi(junction = "1:157797133-157798160", cell = "UCB1.00060", outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/Real_MXE/")

MXE_psi["8283315",  row.names(Cell_type)] -> test
row.names(data.frame(na.omit(t(test))))
sashimi(junction = "1:157797133-157798160", cell = row.names(data.frame(na.omit(t(test)))), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/Real_MXE/", min = 1)


load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_sep_SE.RData")

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))

load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")

Bsep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_Bcell, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_Bcell, "gene_id"]))]
CD4sep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_CD4, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_CD4, "gene_id"]))]
CD8sep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_CD8, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_CD8, "gene_id"]))]
MKsep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_MK, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_MK, "gene_id"]))]
Monosep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_Mono, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_Mono, "gene_id"]))]
NKsep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_NK, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_NK, "gene_id"]))]
Tsep_SE_gene <- unique(SE_psi[SE_psi$loci %in% sj_Tcell, "gene_id"])[!is.na(unique(SE_psi[SE_psi$loci %in% sj_Tcell, "gene_id"]))]
mapply(length, list(Bsep_SE_gene, CD4sep_SE_gene, CD8sep_SE_gene, MKsep_SE_gene, Monosep_SE_gene, NKsep_SE_gene, Tsep_SE_gene))


list(`B Cell` = Bsep_SE_gene, 
	 `CD4+ T Cell` = CD4sep_SE_gene, 
	 `CD8+ T Cell` = CD8sep_SE_gene, 
	 MK = MKsep_SE_gene, 
	 Monocytes = Monosep_SE_gene, 
	 NK = NKsep_SE_gene, 
	 `T Cell` = Tsep_SE_gene) -> gene_list

library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/DSU2/Document/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])


CompareGO_BP <- compareCluster(geneCluster   = gene_list, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = gene_list, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = gene_list, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Cell_type_specific_CompareGO.pdf", height = 16, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "SE GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "SE GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "SE GO CC"))
dev.off()




#### cell-type-specific express gene function enrichment




library(gdata)
sheet1 <- read.xls("/mnt/data5/BGI/UCB/ExpMat_NewID/cluster specific genes.xlsx", sheet=1, stringsAsFactors = F, header = T)
sheet1 <- sheet1[,c(2,7,8)]
sheet2 <- read.xls("/mnt/data5/BGI/UCB/ExpMat_NewID/cluster specific genes.xlsx", sheet=2, stringsAsFactors = F)
sheet2$cluster <- "T cells"
sheet2 <- sheet2[,c(3,7,1)]
colnames(sheet1) <- colnames(sheet2)
gene <- rbind(sheet1,sheet2)
Gene_info <- read.delim("/mnt/data5/BGI/UCB/ExpMat_NewID/Gene_info.txt", header = T, stringsAsFactors = F)

gene <- merge(x = gene, y = Gene_info, by.x = "gene", by.y = "GeneName", all.x = T)


with(gene,GeneID[cluster=="Monocytes"])

gene_list <- list("B cells" = with(gene,GeneID[cluster=="B cells"]),
				  "T cells" = with(gene,GeneID[cluster=="T cells"]),
				  "Megakaryocyte" = with(gene,GeneID[cluster=="Megakaryocyte"]),
				  "Monocytes" = with(gene,GeneID[cluster=="Monocytes"]),
				  "NK cells" = with(gene,GeneID[cluster=="NK cells"]))

library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/DSU2/Document/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

CompareGO_BP <- compareCluster(geneCluster   = gene_list, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = gene_list, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = gene_list, 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/gene_expression/Cell_type_specific/Cell_type_specific_gene_CompareGO_BP.pdf", height = 14, width = 18)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "GO BP"))
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/gene_expression/Cell_type_specific/Cell_type_specific_gene_CompareGO_MF.pdf", height = 10, width = 16)
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "GO MF"))
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/gene_expression/Cell_type_specific/Cell_type_specific_gene_CompareGO_CC.pdf", height = 10, width = 16)
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "GO CC"))
dev.off()


gene_enid_list <- lapply(gene_list, function(x) {bitr(x, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID})
background_enid <- bitr(background, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

CompareKEGG <- compareCluster(geneCluster  = gene_enid_list,
							  universe 	   = background_enid,
							  fun          = "enrichKEGG",
							  organism	   = "hsa", 
							  pvalueCutoff = 0.05)
pdf("/mnt/data5/BGI/UCB/tangchao/gene_expression/Cell_type_specific/Cell_type_specific_gene_CompareKEGG.pdf", height = 8, width = 13)
print(dotplot(CompareKEGG,font.size = 16,title = "KEGG",showCategory = 20))
dev.off()

CompareDO <- compareCluster(geneCluster   = gene_enid_list,
							fun 		  = "enrichDO",
				  	  		ont           = "DO", 
				  	  		pvalueCutoff  = 0.05, 
				  	  		pAdjustMethod = "BH",
				  	  		#universe 	  = background_enid,
				  	  		minGSSize     = 5, 
				  	  		maxGSSize     = 500, 
				  	  		qvalueCutoff  = 0.05,
				  	  		readable      = FALSE)
pdf("/mnt/data5/BGI/UCB/tangchao/gene_expression/Cell_type_specific/Cell_type_specific_gene_CompareDO.pdf", height = 8, width = 12)
print(dotplot(CompareDO,font.size = 16,title = "SE KEGG"))
dev.off()




















