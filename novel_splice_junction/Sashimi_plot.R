
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
dim(sj)

sj <- sj[sj$exon_Shannon_entropy >= 1.8 & 
           !is.na(sj$exon_Relative_entropy) &
           sj$ReadSum >= 10 & 
           sj$overhang >= 10, ]
dim(sj)
# [1] 457282     37

sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)

load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj_merge.RData")

reads <- data.frame(as.data.frame(te)[,-1], row.names = te[[1]], stringsAsFactors=F)

reads[sj[sj$classification=="Novel",]$sj,] -> reads_tu

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t")
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

reads_tu[1:5,1:5]
colnames(reads_tu) <- substr(colnames(reads_tu), 1, 10)

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

FastRemoveMissingValues(reads_tu)
reads_tu[1:5,1:5]

reads_sd <- apply(reads_tu, 1, sd)
reads_sc <- apply(reads_tu, 1, function(x) sum(x>0))
reads_sr <- apply(reads_tu, 1, sum)

summary(reads_sd)
summary(reads_sc)
summary(reads_sr)

data <- log10(as.matrix(reads_tu[reads_sd >= 10 & reads_sc >= 10 & reads_sr >= 15000, row.names(Cell_type)])+1)
Cell_type -> annotation2
annotation2$ReadSum <- colSums(data)
annotation2 <- annotation2[order(annotation2$CellType, annotation2$Individual, annotation2$ReadSum),]
annotation3 <- subset.data.frame(annotation2, select = c(2,1))


data.frame(apply(data,1,function(x) colnames(data)[which.max(x)])) -> sashimi_tab
colnames(sashimi_tab) <- "Cell"



sashimi <- function(junction, cell, left_expand = 1000, right_expand = 1000, outDir = "~/"){
	Coordinates <- data.frame(do.call(rbind,strsplit(junction, split = "[:-]")),stringsAsFactors=F)
	Coordinates$X2 <- abs(as.numeric(Coordinates$X2) - left_expand)
	Coordinates$X3 <- abs(as.numeric(Coordinates$X3) + right_expand)
	window <- paste(Coordinates$X1, ":", Coordinates$X2, "-", Coordinates$X3, sep = "")
	Dir_out <- paste(outDir, cell, "_", junction, sep = "")

	work <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi.sh ", window, cell, Dir_out, "')")
	eval(parse(text = work))

}

sashimi(junction = "1:8940213-8945824", cell = "UCB3.02218")

sashimi(junction = "1:8940213-8945824", cell = "UCB3.02218", left_expand = 1000, right_expand = 1000)




for(i in 1:nrow(sashimi_tab)){
	sashimi(junction = row.names(sashimi_tab)[i], cell = sashimi_tab[i,1])
}

sashimi(junction = "1:8940213-8945824", cell = "UCB3.02218", left_expand = 1000, right_expand = 1000, outDir = "/mnt/data5/BGI/UCB/tangchao/")






#### sashimi2
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




sashimi(junction = "1:198692000-198703000", cell = c("UCB1.00070", "UCB4.01425", "UCB3.00040"))

sashimi(junction = "1:198694141-198696711", cell = c("UCB1.00070", "UCB3.01382", "UCB4.01056", "UCB4.01425", "UCB3.00040", "UCB3.01106", "UCB3.01429"))
sashimi(junction = "1:198694141-198699563", cell = c("UCB1.00070", "UCB3.01102", "UCB3.00465", "UCB3.00596", "UCB3.00481", "UCB3.01382", "UCB4.01056"))
sashimi(junction = "1:198729766-198731616", cell = c("UCB1.00070", "UCB4.01229", "UCB4.01285", "UCB4.01108", "UCB3.00820", "UCB4.00670", "UCB3.00062"))

sashimi(junction = "1:198692374-198703297", cell = "UCB4.00338")


sashimi(junction = "12:69353875-69934754", cell = "UCB1.00016")



python /mnt/data5/BGI/UCB/tangchao/soft/ggsashimi/sashimi-plot.py \
                -b /mnt/data5/BGI/UCB/tangchao/data/BAM_tsv/UCB1.00016.tsv \
                -c 12:69353875-69934754 \
                -g /mnt/data5/BGI/UCB/tangchao/MY.GRCh38.87.for.ggsashimi.gtf \
                -M 11 \
                --alpha 0.3 \
                --base-size=20 \
                --ann-height=12 \
                --height=4 \
                --width=18 \
                -o  ~/UCB1.00016 \
                -C 3



sashimi(junction = "12:102075095-102075193", cell = c("UCB3.00792","UCB1.00060"))
sashimi(junction = "12:53461144-53462480", cell = c("UCB4.01506","UCB4.01510","UCB4.01515","UCB4.01520"))
sashimi(junction = "10:100262064-100280123", cell = c("UCB1.00047", "UCB1.00048", "UCB1.00069"))




