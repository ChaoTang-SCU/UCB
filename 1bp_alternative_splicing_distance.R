
#### 1 bp alternative splicing distancce select ===================================================================================================================

require(parallel)
library(dplyr)
library(plyr)
library(pryr)
require(parallel)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(reshape2)
library(vioplot)

depar <- par()

## path and name of SJ merged RData
pfsj <- file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")

## path of figure/file output
pffo <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/")

## path of RData output
pfro <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/RData/")

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")


#### 1.load data --------------------------------------------------------------------------------------------------------------------------------

load(pfsj)


class(te[[1]])
#[1] "character"
length(te[[1]])
#[1] 2069161


#### 2. Identify alternative splicing events ----------------------------------------------------------------------------------------------------


junc.names=do.call(rbind,strsplit(te[[1]],split="[:-]"))
#alnum-- Letters and Numbers；digit -- Numbers
colnames(junc.names) <- c("chr","start","end")
rownames(junc.names) <- te[[1]]
junc.names=data.frame(junc.names,stringsAsFactors=F)
junc.names$start=as.integer(junc.names$start)
junc.names$end=as.integer(junc.names$end)
junc.names <- junc.names[junc.names$chr %in% c(1:22,"X","Y"),]
junc.names=junc.names[order(junc.names$chr,junc.names$start,junc.names$end),]
junc.names$names=rownames(junc.names)


junc.as=do.call(c,mclapply(unique(junc.names$chr),function(chr) {
  junc.chr=junc.names[junc.names$chr==chr,]
  same.start=dlply(junc.chr,c("start"),function(x) x)
  same.end=dlply(junc.chr,c("end"),function(x) {if(nrow(x)>1) {return(x)} else {return(NULL)}})
  print(chr)
  return(c(same.start,same.end[!sapply(same.end,is.null)]))
},mc.cores=10))


juncs <- junc.as

junc.as <- junc.as[sapply(junc.as,nrow)>1]
names(junc.as) <- sapply(junc.as,function(x) paste(x$names,collapse="_"))

## junction just has two alternative starts/ends:
junc.as.2=junc.as[sapply(junc.as,nrow)==2]

table(sapply(junc.as.2,nrow))

lapply(junc.as.2, function(x) {sum(diff(x$start),diff(x$end))}) -> junc_len
summary(as.numeric(junc_len))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      1    2511   13085   33290   33946 1935525
sum(as.numeric(junc_len)==1)
# [1] 1253
length(as.numeric(junc_len)==1)
# [1] 362642
# 1253/362642
# [1] 0.003455198
hist(as.numeric(junc_len)[as.numeric(junc_len)<1000],breaks = 100, xlab="Junction Length", main = "Histogram of junction length")

data.frame(table(as.numeric(junc_len))[order(table(as.numeric(junc_len)), decreasing=T)][1:40], AllJuns = length(as.numeric(junc_len)), 
	Pro = round(as.numeric(table(as.numeric(junc_len))[order(table(as.numeric(junc_len)), decreasing=T)][1:40])/length(as.numeric(junc_len))*100,2))
   Var1 Freq AllJuns  Pro
1     3 2106  362642 0.58
2     4 1966  362642 0.54
3     2 1726  362642 0.48
4     1 1253  362642 0.35
5     5 1035  362642 0.29
6     6  937  362642 0.26
7     7  522  362642 0.14
8     9  482  362642 0.13
9     8  461  362642 0.13
10   12  437  362642 0.12
11   10  397  362642 0.11
12   15  371  362642 0.10
13   11  344  362642 0.09
14   21  324  362642 0.09
15   18  321  362642 0.09
16  135  313  362642 0.09
17   14  306  362642 0.08
18   13  305  362642 0.08
19   16  300  362642 0.08
20   32  299  362642 0.08
21   24  290  362642 0.08
22   19  284  362642 0.08
23   20  283  362642 0.08
24   17  282  362642 0.08
25   22  259  362642 0.07
26   27  255  362642 0.07
27   23  253  362642 0.07
28  134  248  362642 0.07
29   25  239  362642 0.07
30   26  235  362642 0.06
31   30  223  362642 0.06
32   36  218  362642 0.06
33   42  213  362642 0.06
34   28  211  362642 0.06
35   31  210  362642 0.06
36   29  205  362642 0.06
37   37  203  362642 0.06
38   33  202  362642 0.06
39   39  200  362642 0.06
40   35  197  362642 0.05



junc_len[junc_len==1] -> junc_len1
strsplit(names(junc_len1), split="_") -> junc_len1_junc

te[junc_len1_junc[[1]],] -> test
data.frame(as.data.frame(test)[,-1], row.names = test[[1]]) -> test
rowSums(!is.na(test))
test[,colSums(is.na(test))!=2]
test[,colSums(!is.na(test))==2]


lapply(junc_len1_junc, function(x){
	te[x,] -> test
	data.frame(as.data.frame(test)[,-1], row.names = test[[1]]) -> test
	if(sum(colSums(!is.na(test))==2)>0){
		subset.data.frame(test,select = which(colSums(!is.na(test))==2)) -> res
		return(res)
	}
	}) -> test_result

sum(unlist(lapply(test_result, is.null)))
# 1239
length(test_result)
# [1] 1253

test_result2 <- test_result[!unlist(lapply(test_result, is.null))]
length(test_result2)
# [1] 14




junc_len[junc_len==2] -> junc_len2
strsplit(names(junc_len2), split="_") -> junc_len2_junc

lapply(junc_len2_junc, function(x){
	te[x,] -> test
	data.frame(as.data.frame(test)[,-1], row.names = test[[1]]) -> test
	if(sum(colSums(!is.na(test))==2)>0){
		subset.data.frame(test,select = which(colSums(!is.na(test))==2)) -> res
		return(res)
	}
	}) -> test2_result

sum(unlist(lapply(test2_result, is.null)))
# 1059
length(test2_result)
# [1] 1726

test2_result2 <- test2_result[!unlist(lapply(test2_result, is.null))]
length(test2_result2)
# [1] 667


lapply(lapply(test2_result2, colnames),function(x) grepl(x,pattern="UCB1\\.00248"))
sum(unlist(lapply(lapply(test2_result2, colnames),function(x) grepl(x,pattern="UCB1\\.00248"))))
# [1] 0

save(test_result2, test2_result2, file = "/mnt/data1/tangchao/junc_distance_test.RData")




lapply(test_result2, function(x) {t(x)}) -> test
lapply(test, as.data.frame) -> test
lapply(test, function(x) {x[,order(colSums(x))]}) -> test
for(i in 1:length(test)){
	#test[[i]] <- as.data.frame(test[[i]])
	test[[i]]$junc <- paste(colnames(test[[i]]), collapse = "_")
	test[[i]]$cell <- row.names(test[[i]])
	colnames(test[[i]])[1:2] <- c("minor","major")
}
do.call(rbind, test) -> junc_dis1_tab
table(junc_dis1_tab$minor)
#   1   2   3   4   5   6   7   8  14
# 113  38  18   7   8   3   3   1   1

lapply(test2_result2, function(x) {t(x)}) -> test
lapply(test, as.data.frame) -> test
lapply(test, function(x) {x[,order(colSums(x))]}) -> test
for(i in 1:length(test)){
	#test[[i]] <- as.data.frame(test[[i]])
	test[[i]]$junc <- paste(colnames(test[[i]]), collapse = "_")
	test[[i]]$cell <- row.names(test[[i]])
	colnames(test[[i]])[1:2] <- c("minor","major")
}
do.call(rbind, test) -> junc_dis2_tab
## in the subset of distance equal to 2, there some cells, in which the minor become major situation.




