library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
background <- as.character(read.table("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

BG_BP_ego <- enrichGO(gene = sample(background, 3000, replace = F),
                        #universe = background,
                        OrgDb  = org.Hs.eg.db,
                        keytype= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
BG_BP_ego@result


tang <- list()
for(i in 1:100){
  print(i)
  BG_BP_ego <- enrichGO(gene = sample(background, 3000, replace = F),
                        #universe = background,
                        OrgDb  = org.Hs.eg.db,
                        keytype= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
  tang[[i]] <- BG_BP_ego@result
}

table(unlist(lapply(tang,nrow)))
0  1  2  3  4  5  6  9 
84  6  4  1  2  1  1  1 

tang2 <- do.call(rbind, tang[lapply(tang,nrow)>0])


save(tang, tang2, file = "/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/background_random_selecction.RData")

