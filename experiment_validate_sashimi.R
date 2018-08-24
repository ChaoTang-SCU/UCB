
load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te_no_NA_and_normalized.RData")
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/interested_novel_SJ_and_max_cell.RData")
novel_sj_maxreads -> novel_sj_morethan_onepeople
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/interested_novel_SJ_and_max_cell_only_in_one_people.RData")
novel_sj_maxreads -> novel_sj_onepeople
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/cell_type_specific_unannotated_SJ.RData")
novel_sj_maxreads -> cell_type_sep
all_sj <- unique(c(cell_type_sep$SJ, novel_sj_onepeople$SJ, novel_sj_morethan_onepeople$SJ))
te_tu <- te[all_sj,]
rm("te")
gc()

all_sj_top3 <- t(apply(te_tu,1,function(x) sort(x, decreasing = T)[1:3]))
all_sj_top3_cell <- t(apply(te_tu,1,function(x) names(sort(x, decreasing = T))[1:3]))


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


for(i in 1:nrow(all_sj_top3_cell)){
  print(i)
  if(row.names(all_sj_top3_cell)[i] %in% cell_type_sep$SJ){
  		sashimi(junction = row.names(all_sj_top3_cell)[i], cell = all_sj_top3_cell[i,], 
    			outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/cell_type_sep/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
  }
  if(row.names(all_sj_top3_cell)[i] %in% novel_sj_morethan_onepeople$SJ){
  		sashimi(junction = row.names(all_sj_top3_cell)[i], cell = all_sj_top3_cell[i,], 
    			outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/novel_sj_morethan_onepeople/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
  }
  if(row.names(all_sj_top3_cell)[i] %in% novel_sj_onepeople$SJ){
  		sashimi(junction = row.names(all_sj_top3_cell)[i], cell = all_sj_top3_cell[i,], 
    			outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/novel_sj_onepeople/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
  }

}

for(i in 1:nrow(all_sj_top3_cell)){
  print(i)
  if(row.names(all_sj_top3_cell)[i] %in% novel_sj_onepeople$SJ){
  		sashimi(junction = row.names(all_sj_top3_cell)[i], cell = all_sj_top3_cell[i,], 
    			outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/novel_sj_onepeople/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
  }

}




load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
test <- cbind(all_sj_top3_cell,all_sj_top3)
test_allinfo <- merge(x=sj,y=test,by.x="sj",by.y="row.names",all.y=T)
write.csv(test_allinfo, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/SJ_information_for_experiment_vilidate.csv")












load(file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/bulk_cd8_sep_sj.RData")
load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te_no_NA_and_normalized.RData")

te_bulksep_tu <- te[bulk_na_sj,]
rm("te")
gc()

all_sj_top3 <- t(apply(te_bulksep_tu,1,function(x) sort(x, decreasing = T)[1:3]))
all_sj_top3_cell <- t(apply(te_bulksep_tu,1,function(x) names(sort(x, decreasing = T))[1:3]))

for(i in 1:nrow(all_sj_top3_cell)){
  print(i)
  sashimi(junction = row.names(all_sj_top3_cell)[i], cell = all_sj_top3_cell[i,], 
    	  outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/BULK_CD8_Sep/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
}


all_sj_top3_cell[all_sj_top3_cell%in%all_sj]

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", row.names = 1, header = F, sep = "\t", stringsAsFactors=F)

Cell_type[names(te_bulksep_tu[1,te_bulksep_tu[1,]>0])[names(te_bulksep_tu[1,te_bulksep_tu[1,]>0])%in%row.names(Cell_type)],]
table(Cell_type[names(te_bulksep_tu[1,te_bulksep_tu[1,]>0])[names(te_bulksep_tu[1,te_bulksep_tu[1,]>0])%in%row.names(Cell_type)],])

apply(te_bulksep_tu, 1, function(x) table(Cell_type[names(x[x>0])[names(x[x>0])%in%row.names(Cell_type)],]))











#### CD45 validate ==============================================================================================================================



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

cov_tu[,1:4]
one_isoform <- cov_tu[,colSums(cov_tu>0)==1]
one_isoform_todo <- one_isoform[,names(tail(sort(colSums(one_isoform)),10))]

two_isoform <- cov_tu[,colSums(cov_tu>100)==2 & colSums(cov_tu>0)==2]
two_isoform_todo <- two_isoform[,names(tail(sort(colSums(two_isoform)),10))]

three_isoform <- cov_tu[,colSums(cov_tu>100)==3 & colSums(cov_tu>0)==3]
three_isoform_todo <- three_isoform[,names(tail(sort(colSums(three_isoform)),10))]

four_isoform <- cov_tu[,colSums(cov_tu>100)==4 & colSums(cov_tu>0)==4]
four_isoform_todo <- four_isoform[,names(tail(sort(colSums(four_isoform)),10))]

five_isoform <- cov_tu[,colSums(cov_tu>0)==5]
five_isoform_todo <- five_isoform[,names(tail(sort(colSums(five_isoform)),10))]



sashimi(junction = "1:198692374-198703297", cell = colnames(one_isoform_todo), 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_one_isoform/", 
		ann_height = 5, left_expand = 100, right_expand = 100, min = 1)
write.csv(one_isoform_todo, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_one_isoform/one_isoform_todo.csv")

sashimi(junction = "1:198692374-198703297", cell = colnames(two_isoform_todo), 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_two_isoform/", 
		ann_height = 5, left_expand = 100, right_expand = 100, min = 1)
write.csv(two_isoform_todo, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_two_isoform/two_isoform_todo.csv")

sashimi(junction = "1:198692374-198703297", cell = colnames(three_isoform_todo), 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_three_isoform/", 
		ann_height = 5, left_expand = 100, right_expand = 100, min = 1)
write.csv(three_isoform_todo, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_three_isoform/three_isoform_todo.csv")

sashimi(junction = "1:198692374-198703297", cell = colnames(four_isoform_todo), 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_four_isoform/", 
		ann_height = 5, left_expand = 100, right_expand = 100, min = 1)
write.csv(four_isoform_todo, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_four_isoform/four_isoform_todo.csv")

sashimi(junction = "1:198692374-198703297", cell = colnames(five_isoform_todo), 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_five_isoform/", 
		ann_height = 5, left_expand = 100, right_expand = 100, min = 1)
write.csv(five_isoform_todo, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/CD45_five_isoform/five_isoform_todo.csv")



sashimi(junction = "2:68394114-68397248", cell = row.names(Cell_type[Cell_type$CellType == "Megakaryocytes",]), 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/MK_interesting/", 
		ann_height = 4, left_expand = 100, right_expand = 100, min = 1)

sashimi(junction = "2:68394114-68398248", cell = row.names(Cell_type[Cell_type$CellType == "Megakaryocytes",])[1], 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/MK_interesting/", 
		ann_height = 4, left_expand = 100, right_expand = 100, min = 1)

for(i in 1:length(row.names(Cell_type[Cell_type$CellType == "Megakaryocytes",]))){
	sashimi(junction = "2:68394114-68398248", cell = row.names(Cell_type[Cell_type$CellType == "Megakaryocytes",])[i], 
		outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/MK_interesting/", 
		ann_height = 4, left_expand = 100, right_expand = 100, min = 1)
}













