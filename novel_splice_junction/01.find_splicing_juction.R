
#### novel splice sites 
## Fri Mar 23 19:20:05 CST 2018
## By Tang Chao

setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS")

library("data.table")

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames = list.files(path=mypath, pattern="tab", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp <- unique(fread(x,header=FALSE, select = c(1,2,3,7,9)))
    tmp <- tmp[V1 %in% c(1:22,"X","Y") & V9>2 & V7>1,]
    tmp$name <- paste(tmp[[1]],":",tmp[[2]],"-",tmp[[3]],sep="")
    tmp <- tmp[,c("name","V7")]
    setnames(tmp, "V7",tail(strsplit(x,"/")[[1]],n=1))
    setkey(tmp,name)
    return(tmp)})
  Reduce(function(x,y) {merge(x,y,all=T,by="name")}, datalist)
}
mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134")
system.time(te<-multmerge(mypath))


mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134")
filenames = list.files(path=mypath, pattern="tab", full.names=TRUE)
for(i in 1:length(filenames)){
	if(i%%100 == 0) print(paste(i, "of", length(filenames), "--", date()))
	if(i == 1){
		tmp <- unique(fread(filenames[i],header=FALSE, select = c(1,2,3,7,9)))
    	tmp <- tmp[V1 %in% c(1:22,"X","Y") & V9>2 & V7>1,]
    	tmp$name <- paste(tmp[[1]],":",tmp[[2]],"-",tmp[[3]],sep="")
    	tmp <- tmp[,c("name","V7")]
    	setnames(tmp, "V7",tail(strsplit(filenames[i],"/")[[1]],n=1))
    	setkey(tmp,name)
    	te <- tmp
	}else{
		tmp <- unique(fread(filenames[i],header=FALSE, select = c(1,2,3,7,9)))
    	tmp <- tmp[V1 %in% c(1:22,"X","Y") & V9>2 & V7>1,]
    	tmp$name <- paste(tmp[[1]],":",tmp[[2]],"-",tmp[[3]],sep="")
    	tmp <- tmp[,c("name","V7")]
    	setnames(tmp, "V7",tail(strsplit(filenames[i],"/")[[1]],n=1))
    	setkey(tmp,name)
    	te <- merge(x = te, y = tmp, all=T, by="name")
	}
	gc()	
}
save(te, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/sj_merge_.RData")


SJ_dic <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/SJ_dic.txt", header = F, stringsAsFactors = F)
table(SJ_dic$V4)/nrow(SJ_dic)
#          0          1
# 0.90389341 0.09610659
SJ_tu <- SJ_dic[SJ_dic$V1 %in% te[[1]],]
table(SJ_tu$V4)
#      0      1
# 693522 176529
table(SJ_tu$V4)/nrow(SJ_tu)
#        0        1
# 0.797105 0.202895


sjdb <- read.table("/mnt/data1/reference/ensembl/human/star_index/sjdbList.out.tab", header = F, stringsAsFactors = F)
sjdb$name <- paste(sjdb$V1, ":", sjdb$V2, "-", sjdb$V3, sep = "")

sum(SJ_tu$V1 %in% sjdb$name) == nrow(SJ_tu[SJ_tu$V4 == 1,])
[1] TRUE
## It's true that only 20 percent of the SJs have been annotated in sjdb.


data.frame(do.call(rbind,strsplit(SJ_tu$V1, split = "[:-]")), stringsAsFactors=F) -> loc
sj <- data.frame(sj = SJ_tu$V1, chr = loc$X1, start = loc$X2, end = loc$X3, strand = SJ_tu$V2, motif = SJ_tu$V3, annotation = SJ_tu$V4, stringsAsFactors = F)

## Read sum and Cell sum 
te_tab <- data.frame(as.data.frame(te)[,-1], row.names = te[[1]], stringsAsFactors = F)

te_tab[as.character(sj$sj), ] -> te_tab_tu

sj$ReadSum <- rowSums(te_tab_tu, na.rm = T)
sj$CellSum <- rowSums(!is.na(te_tab_tu))


save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/SJ_tu.RData")
write.table(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/SJ_tu.txt", row.names = F, quote = F, sep = "\t")












