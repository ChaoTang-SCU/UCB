#### Novel splicing junction only in one cell:
load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
nsj <- sj[sj$novel == "Y","sj"]
load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
te <- te[name %in% nsj]
te <- te[rowSums(!is.na(te[,2:3575]))==1]

tab <- data.frame(as.data.frame(te[,-1]), row.names = te[[1]])
tabsum <- rowSums(tab, na.rm=T)
tab <- tab[tabsum>3000,]
tab <- as.data.frame(apply(tab,1,function(x)names(which(!is.na(x)))))
colnames(tab) <- "cell"
tab$cell <- substr(tab$cell, 1, 10)

save(tab, file = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/novel_sj_onecell/novel_sj_onecell.RData")
sj <- sj[sj$sj %in% row.names(tab),]
write.csv(sj, "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/novel_sj_onecell/novel_sj_onecell.csv")

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


for(i in 1:nrow(tab)){
  print(i)
  		sashimi(junction = row.names(tab)[i], cell = tab[i,1], 
    			outDir = "/mnt/data5/BGI/UCB/tangchao/for_vilidate/sashimi/novel_sj_onecell/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
}

