## SE only in one cell
load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_type_and_PSI.RData")
SE_OneCell <- all_cells_SE_type_and_PSI[as.numeric(rowSums(!is.na(all_cells_SE_type_and_PSI[,12:ncol(all_cells_SE_type_and_PSI)]))) == 1, ]
SE_OneCell_exactly <- SE_OneCell[SE_OneCell$AStype == "exon-skipping_exactly", ]
SE_OneCell_intergenic <- SE_OneCell[SE_OneCell$AStype == "intergenic-splicing_site", ]

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
nsj <- sj[sj$novel == "Y","sj"]

sum(do.call(rbind,strsplit(as.character(SE_OneCell_intergenic$left), split="_"))[,1] %in% nsj)
#[1] 633
sum(do.call(rbind,strsplit(as.character(SE_OneCell_intergenic$right), split="_"))[,2] %in% nsj)
#[1] 683
sum(do.call(rbind,strsplit(as.character(SE_OneCell_intergenic$left), split="_"))[,1] %in% nsj & do.call(rbind,strsplit(as.character(SE_OneCell_intergenic$right), split="_"))[,2] %in% nsj)
#[1] 439

SE_OneCell_intergenic <- SE_OneCell_intergenic[do.call(rbind,strsplit(as.character(SE_OneCell_intergenic$left), split="_"))[,1] %in% nsj & do.call(rbind,strsplit(as.character(SE_OneCell_intergenic$right), split="_"))[,2] %in% nsj, ]
write.csv(SE_OneCell_intergenic, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/SE_OneCell_intergenic_novel.csv")

psi <- as.numeric(apply(SE_OneCell_intergenic[,12:ncol(SE_OneCell_intergenic)],1,function(x) x[!is.na(x)]))
cell <- as.character(apply(SE_OneCell_intergenic[,12:ncol(SE_OneCell_intergenic)],1,function(x) names(x[!is.na(x)])))

tab <- data.frame(SE_OneCell_intergenic[,1:11], psi = psi, cell = cell, stringsAsFactors = F)
sum(tab$ASlen > 50 & tab$ASlen < 1000 & tab$psi < 0.8 & tab$psi > 0.2)
#[1] 123
tab_tu <- tab[tab$ASlen > 1 & tab$ASlen < 1000 & tab$psi < 0.8 & tab$psi > 0.2, ]

tab_tu$range <- paste(do.call(rbind,strsplit(tab_tu$loci,"-"))[,1],"-",do.call(rbind,strsplit(tab_tu$loci,"-"))[,4], sep = "")

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


sashimi(junction = "14:22905454-22906194", cell = "UCB1.00003", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "6:107746890-107749683", cell = "UCB1.00003", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "7:74246764-74249011", cell = "UCB1.00012", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "22:39376198-39377597", cell = "UCB1.00013", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "12:109939248-109953091", cell = "UCB1.00010", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "3:69071618-69077797", cell = "UCB1.00012", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "15:44499814-44514597", cell = "UCB1.00012", 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)



for(i in 1:nrow(tab_tu)){
	print(paste(i, "of", nrow(tab_tu)))
	sashimi(junction = tab_tu$range[i], cell = tab_tu$cell[i], 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SE_OneCell/novel", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
}




