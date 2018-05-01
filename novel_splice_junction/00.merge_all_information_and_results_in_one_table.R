#### I want to merge all information and results in one table, so that we can plot all figures from this immediately!
## Mon Mar 26 2018
## By: Tang Chao


######## 1: merge SJ reads from STAR result =====================================================================================================


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


save(te, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj_merge.RData")

data.frame(as.data.frame(te)[,-1], row.names = te[[1]], stringsAsFactors=F) -> reads

identical(row.names(overhang),row.names(reads))
# [1] TRUE


SJ_dic <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/SJ_dic.txt", header = F, stringsAsFactors = F)
table(SJ_dic$V4)/nrow(SJ_dic)
#          0          1
# 0.90389341 0.09610659
SJ_tu <- SJ_dic[SJ_dic$V1 %in% row.names(reads),]
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
sj$ReadSum <- rowSums(reads, na.rm = T)
sj$CellSum <- rowSums(!is.na(reads))

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj1.RData")


######## 2: merge SJ overhang from STAR result ==================================================================================================


for(i in 1:length(filenames)){
	if(i%%100 == 0) print(paste(i, "of", length(filenames), "--", date()))
	if(i == 1){
		tmp <- unique(fread(filenames[i],header=FALSE, select = c(1,2,3,7,9)))
    	tmp <- tmp[V1 %in% c(1:22,"X","Y") & V9>2 & V7>1,]
    	tmp$name <- paste(tmp[[1]],":",tmp[[2]],"-",tmp[[3]],sep="")
    	tmp <- tmp[,c("name","V9")]
    	setnames(tmp, "V9",tail(strsplit(filenames[i],"/")[[1]],n=1))
    	setkey(tmp,name)
    	te <- tmp
	}else{
		tmp <- unique(fread(filenames[i],header=FALSE, select = c(1,2,3,7,9)))
    	tmp <- tmp[V1 %in% c(1:22,"X","Y") & V9>2 & V7>1,]
    	tmp$name <- paste(tmp[[1]],":",tmp[[2]],"-",tmp[[3]],sep="")
    	tmp <- tmp[,c("name","V9")]
    	setnames(tmp, "V9",tail(strsplit(filenames[i],"/")[[1]],n=1))
    	setkey(tmp,name)
    	te <- merge(x = te, y = tmp, all=T, by="name")
	}
	gc()	
}

save(te, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/overhang_merge.RData")
data.frame(as.data.frame(te)[,-1], row.names = te[[1]], stringsAsFactors=F) -> overhang
identical(sj$sj, row.names(overhang))
# [1] TRUE

sj$overhang <- apply(overhang,1,function(x) max(x, na.rm = T))
save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj2.RData")
rm(list = ls())
gc()


######## 3: Find single cell specific and novel splice junction =================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj2.RData")
mrna <- read.table("/mnt/data1/projects/UCB/results_ucb/tangchao/GMAP/mRNA_intron.sj", stringsAsFactors = F)
est <- read.table("/mnt/data1/projects/UCB/results_ucb/tangchao/GMAP/est_intron.sj", stringsAsFactors = F)
bulk <- read.table("/mnt/data1/projects/UCB/data_ucb/TXT/bulk_AS_events_dictionary.txt",header=T, stringsAsFactors = F)

#### 2. overview of four methods' SJ 

mRNA_SJ <- as.character(mrna$V1)
EST_SJ <- as.character(est$V1)
bulk_SJ <- as.character(bulk$event_id)
sc_SJ <- as.character(sj$sj)

sc_sj <- sj[!sj$sj %in% unique(c(mRNA_SJ, EST_SJ, bulk_SJ)),]

table(sc_sj$annotation)
#      0      1
# 613535   7018

sj$SCSpe <- "N"
sj$novel <- "N"
sj[sj$sj %in% sc_sj$sj, ]$SCSpe <- "Y"
sj[sj$sj %in% sc_sj[sc_sj$annotation==0,]$sj, ]$novel <- "Y"

sj$SCSpe <- as.factor(sj$SCSpe)
sj$novel <- as.factor(sj$novel)

head(sj)
table(sj$SCSpe)
table(sj$novel)
dim(sj)

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj3.RData")
rm(list = ls())
gc()


######## 4: Call gene sequence from genome fastq file ===========================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj3.RData")
write.table(sj,"/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj3.txt", sep = "\t", row.names = F, quote = F)


#### shell
cd /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table
for f in `awk '{print $2":"$3-50"-"$3+49}' ./sj3.txt | sed -n '2,$p'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./start_site.fa
done

sed ':a;N;$!ba;s/\n/ /g' start_site.fa | sed 's/>/\n/g' | awk -F " " '{print $1"\t"$2$3}' > start_site.txt
sed -i '1d' start_site.txt


# right flank exons FASTA
for f in `awk '{print $2":"$4-49"-"$4+50}' ./sj3.txt | sed -n '2,$p'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./end_site.fa
done

sed ':a;N;$!ba;s/\n/ /g' end_site.fa | sed 's/>/\n/g' | awk -F " " '{print $1"\t"$2$3}' > end_site.txt
sed -i '1d' end_site.txt


paste start_site.txt end_site.txt | awk '{print $1"\t"$3"\t"$2$4}' > SS.fa.txt
####


#### R
fa <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/SS.fa.txt", sep = "\t", header = F, stringsAsFactors = F)

sj$fa <- as.character(fa$V3)

sj$ParsedFasta <- as.character(fa$V3)

library(Biostrings)
## Get the Reverse complementary chain

for (i in 1:nrow(sj)){
	if(sj[i, "strand"] == 2){
		sj[i, "ParsedFasta"] <- as.character(reverseComplement(DNAStringSet(sj[i, "fa"])))
	}
} 

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj4.RData")
rm(list = ls())
gc()


######## 5: MaxEntScan for Maximum Entropy of splicing sites ====================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj4.RData")

remain = grepl(sj$ParsedFasta, pattern = "N")

sj <- sj[!remain, ]

donor <- substr(with(sj,ParsedFasta),48,56)
acceptor <- substr(with(sj,ParsedFasta),131,153)
for_MaxEnt <- as.data.frame(cbind(donor,acceptor), stringsAsFactors=F)


write.table(for_MaxEnt$donor, "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/donor_for_MaxEnt.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(for_MaxEnt$acceptor, "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/acceptor_for_MaxEnt.txt", row.names = F, col.names = F, quote = F, sep = "\t")


#### shell

#### MaxEntScan for Maximum Entropy:


cd /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/MaxEntScan/
perl score5.pl /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/donor_for_MaxEnt.txt > /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/donor_MaxEnt.txt
perl score3.pl /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/acceptor_for_MaxEnt.txt > /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/acceptor_MaxEnt.txt

####


#### R 

sj$donor_MaxEnt <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/donor_MaxEnt.txt", header = F, stringsAsFactors = F)$V2
sj$acceptor_MaxEnt <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/acceptor_MaxEnt.txt", header = F, stringsAsFactors = F)$V2

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj5.RData")
rm(list = ls())
gc()


######## 6: Stop codon around splicing sites ====================================================================================================


library(Biostrings)
## Get the Reverse complementary chain

stopCoden = function(fasta, strand = NULL){
  library(Biostrings)
  remain = grepl(fasta, pattern = "N") == FALSE
  
  fasta <- fasta[remain]
  strand <- strand[remain]
  
  if(is.null(strand)){
    output <- data.frame(sum = NA, min = NA)
    
    for (i in (1:length(fasta))) {
      seq<-fasta[i]
      seq_orf1<-seq
      seq_orf2<-substring(seq,2,nchar(as.character(seq)))
      seq_orf3<-substring(seq,3,nchar(as.character(seq)))
      
      t_orf1<-DNAStringSet(seq_orf1)
      t_orf2<-DNAStringSet(seq_orf2)
      t_orf3<-DNAStringSet(seq_orf3)
      
      aa_orf1<-as.character(translate(t_orf1))
      aa_orf2<-as.character(translate(t_orf2))
      aa_orf3<-as.character(translate(t_orf3))
      
      freq1 <- cbind(aa_orf1,as.numeric(length(which(strsplit(aa_orf1, "")[[1]] == "*"))))
      freq2 <- cbind(aa_orf2,length(which(strsplit(aa_orf2, "")[[1]] == "*")))
      freq3 <- cbind(aa_orf3,length(which(strsplit(aa_orf3, "")[[1]] == "*")))
      
      freq<-rbind(freq1,freq2,freq3)
      
      ####reverse
      seq_rorf1<-reverseComplement(t_orf1)
      seq_rorf2<-substring(reverseComplement(t_orf1),2,nchar(as.character(seq)))
      seq_rorf3<-substring(reverseComplement(t_orf1),3,nchar(as.character(seq)))
      
      t_rorf1<-DNAStringSet(seq_rorf1)
      t_rorf2<-DNAStringSet(seq_rorf2)
      t_rorf3<-DNAStringSet(seq_rorf3)
      
      aa_rorf1<-as.character(translate(t_rorf1))
      aa_rorf2<-as.character(translate(t_rorf2))
      aa_rorf3<-as.character(translate(t_rorf3))
      
      rfreq1 <- cbind(aa_rorf1,length(which(strsplit(aa_rorf1, "")[[1]] == "*")))
      rfreq2 <- cbind(aa_rorf2,length(which(strsplit(aa_rorf2, "")[[1]] == "*")))
      rfreq3 <- cbind(aa_rorf3,length(which(strsplit(aa_rorf3, "")[[1]] == "*")))
      
      rfreq<-rbind(rfreq1,rfreq2,rfreq3)
      
      stop_count_toadd <- sum(as.numeric(rfreq[,2]))+sum(as.numeric(freq[,2]))
      stop_coden_min <- min(c(as.numeric(rfreq[,2]),as.numeric(freq[,2])))
      output[i,"sum"] <- stop_count_toadd
      output[i,"min"] <- stop_coden_min
    }
    return(output)
    
  }else{
    output <- data.frame(sum = NA, min = NA)
    for (i in (1:length(fasta))){
      if(strand[i] != 2){
        seq<-fasta[i]
      }else{
        seq<-as.character(reverseComplement(DNAStringSet(fasta[i])))
      }
      
      seq_orf1<-seq
      seq_orf2<-substring(seq,2,nchar(as.character(seq)))
      seq_orf3<-substring(seq,3,nchar(as.character(seq)))
      
      t_orf1<-DNAStringSet(seq_orf1)
      t_orf2<-DNAStringSet(seq_orf2)
      t_orf3<-DNAStringSet(seq_orf3)
      
      aa_orf1<-as.character(translate(t_orf1))
      aa_orf2<-as.character(translate(t_orf2))
      aa_orf3<-as.character(translate(t_orf3))
      
      freq1 <- cbind(aa_orf1,as.numeric(length(which(strsplit(aa_orf1, "")[[1]] == "*"))))
      freq2 <- cbind(aa_orf2,length(which(strsplit(aa_orf2, "")[[1]] == "*")))
      freq3 <- cbind(aa_orf3,length(which(strsplit(aa_orf3, "")[[1]] == "*")))
      
      freq<-rbind(freq1,freq2,freq3)
      
      stop_count_toadd <- sum(as.numeric(freq[,2]))
      stop_coden_min <- min(as.numeric(freq[,2]))
      output[i,"sum"] <- stop_count_toadd
      output[i,"min"] <- stop_coden_min
      
    }
    return(output)
  }
}


#### 1. Load data 


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj5.RData")

library(parallel)
donor_exon <- mclapply(substr(sj$ParsedFasta,1,50), function(x) stopCoden(x, strand = NULL), mc.cores=10)
donor_intron <- mclapply(substr(sj$ParsedFasta,51,100), function(x) stopCoden(x, strand = NULL), mc.cores=10)
acceptor_intron <- mclapply(substr(sj$ParsedFasta,101,150), function(x) stopCoden(x, strand = NULL), mc.cores=10)
acceptor_exon <- mclapply(substr(sj$ParsedFasta,151,200), function(x) stopCoden(x, strand = NULL), mc.cores=10)

donor_exon <- do.call(rbind,donor_exon)
donor_intron <- do.call(rbind,donor_intron)
acceptor_intron <- do.call(rbind,acceptor_intron)
acceptor_exon <- do.call(rbind,acceptor_exon)

colnames(donor_exon) <- c("donor_exon_sum_stop_codon", "donor_exon_min_stop_codon")
colnames(donor_intron) <- c("donor_intron_sum_stop_codon", "donor_intron_min_stop_codon")
colnames(acceptor_intron) <- c("acceptor_intron_sum_stop_codon", "acceptor_intron_min_stop_codon")
colnames(acceptor_exon) <- c("acceptor_exon_sum_stop_codon", "acceptor_exon_min_stop_codon")

cbind(donor_exon,donor_intron,acceptor_intron,acceptor_exon) -> stop_codon
sj <- cbind(sj, stop_codon)

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj6.RData")


######## 7: Is validated by EST or mRNA =========================================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj6.RData")
mrna <- read.table("/mnt/data1/projects/UCB/results_ucb/tangchao/GMAP/mRNA_introns_2017_12_29.txt", stringsAsFactors = F)
mrna_bed <- as.data.frame(do.call(rbind,strsplit(mrna$V2,split="[:\\..]")), stringsAsFactors = F)
mrna_bed <- data.frame(chr = mrna_bed$V1, start = as.numeric(mrna_bed$V2)+1, end = as.numeric(mrna_bed$V4)-1)
mrna_bed$name <- paste(mrna_bed$chr, ":", mrna_bed$start, "-", mrna_bed$end, sep = "")
mrna_bed <- mrna_bed[order(mrna_bed$chr, mrna_bed$start, mrna_bed$end),]
mrna_bed_list <- split(mrna_bed, mrna_bed$name)
mrna_junction <- as.data.frame(unlist(lapply(mrna_bed_list,nrow)))
colnames(mrna_junction) <- "n_mRNA"

est <- read.table("/mnt/data1/projects/UCB/results_ucb/tangchao/GMAP/est_introns_2017_12_29.txt", stringsAsFactors = F)
est_bed <- as.data.frame(do.call(rbind,strsplit(est$V2,split="[:\\..]")), stringsAsFactors = F)
est_bed <- data.frame(chr = est_bed$V1, start = as.numeric(est_bed$V2)+1, end = as.numeric(est_bed$V4)-1)
est_bed$name <- paste(est_bed$chr, ":", est_bed$start, "-", est_bed$end, sep = "")
est_bed <- est_bed[order(est_bed$chr, est_bed$start, est_bed$end),]
est_bed_list <- split(est_bed, est_bed$name)
est_junction <- as.data.frame(unlist(lapply(est_bed_list,nrow)))
colnames(est_junction) <- "n_EST"

library("data.table")
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames = list.files(path=mypath, pattern="tab", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp <- unique(fread(x,header=TRUE))
    setkey(tmp,event_id)
    return(tmp)})
  Reduce(function(x,y) {merge(x,y,all=T,by="event_id")}, datalist)
}
mypath<-file.path("/mnt/data1/projects/UCB/data_ucb/BP_ucb_bulk_SJ/")
system.time(te<-multmerge(mypath))
bulk_junction <- data.frame(n_bulk = rowSums(as.data.frame(te)[,-1], na.rm = T), row.names = substr(te[[1]],4,100))

save(mrna_junction, est_junction, bulk_junction, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/EST_mRNA_bulk_junction.RData")

mrna_junction$sj <- row.names(mrna_junction)
est_junction$sj <- row.names(est_junction)
bulk_junction$sj <- row.names(bulk_junction)

class(sj$sj)
sum(sj$sj %in% row.names(mrna_junction))
# [1] 72926
sum(sj$sj %in% row.names(est_junction))
# [1] 90372
sum(sj$sj %in% row.names(bulk_junction))
# [1] 243550

sj <- merge(x = sj, y = mrna_junction, by = "sj", all.x = T)
sj <- merge(x = sj, y = est_junction, by = "sj", all.x = T)
sj <- merge(x = sj, y = bulk_junction, by = "sj", all.x = T)

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj7.RData")
rm(list = ls())


######## 8: The closest splicing points =========================================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj7.RData")
sj$start <- as.numeric(sj$start)
sj$end <- as.numeric(sj$end)
sj <- sj[order(sj$chr, sj$start, sj$end),]

test <- c(diff(sj$start),0)
summary(test)
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
#-248905317          9        262         65       1552   30163750
sum(test<0)
# [1] 23

sj$closest_start <- c(diff(sj$start),0)
sj[sj$closest_start<0, ]$closest_start <- 0

sj$closest_end <- c(diff(sj$end),0)
sj[sj$closest_end<0, ]$closest_end <- 0

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj8.RData")


######## 9: GC content of sequence ==============================================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj8.RData")

## 
library(Biostrings)
letterFrequency(DNAString(sj[1,"fa"]), "GC", as.prob=TRUE)

letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],1,50)), "GC", as.prob=TRUE)
letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],51,100)), "GC", as.prob=TRUE)
letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],101,150)), "GC", as.prob=TRUE)
letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],151,200)), "GC", as.prob=TRUE)

sj$donor_exon_GC <- as.vector(unlist(mclapply(sj$ParsedFasta, function(x){
												letterFrequency(DNAString(substr(x,1,50)), "GC", as.prob=TRUE)
												}, mc.cores = 30)))

sj$donor_intron_GC <- as.vector(unlist(mclapply(sj$ParsedFasta, function(x){
												letterFrequency(DNAString(substr(x,51,100)), "GC", as.prob=TRUE)
												}, mc.cores = 30)))

sj$acceptor_exon_GC <- as.vector(unlist(mclapply(sj$ParsedFasta, function(x){
												letterFrequency(DNAString(substr(x,101,150)), "GC", as.prob=TRUE)
												}, mc.cores = 30)))

sj$acceptor_intron_GC <- as.vector(unlist(mclapply(sj$ParsedFasta, function(x){
												letterFrequency(DNAString(substr(x,151,200)), "GC", as.prob=TRUE)
												}, mc.cores = 30)))


######## 10: Entropy of sequence ================================================================================================================


## function:
MyEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log2(x))
  Reduce("+", -tmp)
}

round(apply(annotated_donor_percentage, 2, MyEntropy),2)


MyRelativeEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log(x/0.25))
  Reduce("+", tmp)
}


tab <- as.data.frame(matrix(NA,nrow=4,ncol=1),row.names=c("A", "C", "G", "T"))
tab[1,] <- letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],1,50)), "A", as.prob=TRUE)
tab[2,] <- letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],1,50)), "C", as.prob=TRUE)
tab[3,] <- letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],1,50)), "G", as.prob=TRUE)
tab[4,] <- letterFrequency(DNAString(substr(sj[1,"ParsedFasta"],1,50)), "T", as.prob=TRUE)


tab[1,] <- letterFrequency(DNAString(paste(substr(sj[1,"ParsedFasta"],1,50), substr(sj[1,"ParsedFasta"],151,200), sep = "")), "A", as.prob=TRUE)
tab[2,] <- letterFrequency(DNAString(paste(substr(sj[1,"ParsedFasta"],1,50), substr(sj[1,"ParsedFasta"],151,200), sep = "")), "C", as.prob=TRUE)
tab[3,] <- letterFrequency(DNAString(paste(substr(sj[1,"ParsedFasta"],1,50), substr(sj[1,"ParsedFasta"],151,200), sep = "")), "G", as.prob=TRUE)
tab[4,] <- letterFrequency(DNAString(paste(substr(sj[1,"ParsedFasta"],1,50), substr(sj[1,"ParsedFasta"],151,200), sep = "")), "T", as.prob=TRUE)

apply(tab, 2, MyRelativeEntropy)
apply(subset.data.frame(tab, tab>0), 2, MyEntropy)


exon_entropy <- mclapply(sj$ParsedFasta, function(x){
							tab <- as.data.frame(matrix(NA,nrow=4,ncol=1),row.names=c("A", "C", "G", "T"))
							tab[1,] <- letterFrequency(DNAString(paste(substr(x,1,50), substr(x,151,200), sep = "")), "A", as.prob=TRUE)
							tab[2,] <- letterFrequency(DNAString(paste(substr(x,1,50), substr(x,151,200), sep = "")), "C", as.prob=TRUE)
							tab[3,] <- letterFrequency(DNAString(paste(substr(x,1,50), substr(x,151,200), sep = "")), "G", as.prob=TRUE)
							tab[4,] <- letterFrequency(DNAString(paste(substr(x,1,50), substr(x,151,200), sep = "")), "T", as.prob=TRUE)
							cbind(apply(tab, 2, MyRelativeEntropy),apply(subset.data.frame(tab, tab>0), 2, MyEntropy))
							}, mc.cores = 30)

intron_entropy <- mclapply(sj$ParsedFasta, function(x){
							tab <- as.data.frame(matrix(NA,nrow=4,ncol=1),row.names=c("A", "C", "G", "T"))
							tab[1,] <- letterFrequency(DNAString(substr(x,51,150)), "A", as.prob=TRUE)
							tab[2,] <- letterFrequency(DNAString(substr(x,51,150)), "C", as.prob=TRUE)
							tab[3,] <- letterFrequency(DNAString(substr(x,51,150)), "G", as.prob=TRUE)
							tab[4,] <- letterFrequency(DNAString(substr(x,51,150)), "T", as.prob=TRUE)
							cbind(apply(tab, 2, MyRelativeEntropy),apply(subset.data.frame(tab, tab>0), 2, MyEntropy))
							}, mc.cores = 30)

exon_entropy <- do.call(rbind,exon_entropy)
intron_entropy <- do.call(rbind,intron_entropy)

# If MyRelativeEntropy == NaN, means that this sequence only have 3 (or less than 3) base type.

sj$exon_Shannon_entropy <- exon_entropy[,2]
sj$exon_Relative_entropy <- exon_entropy[,1]

sj$intron_Shannon_entropy <- intron_entropy[,2]
sj$intron_Relative_entropy <- intron_entropy[,1]

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
dim(sj)
# [1] 870046     37


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

write.table(sj, paste(to, "sj_tu.txt", sep = ""), sep = "\t", row.names = F, quote = F)



######## 11: Exon Length of SJ ==================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
dim(sj)

test <- sj[,1:5]
test$ub <- NA
test$db <- NA
for (i in 1:nrow(test)) {
  print(paste(i, "of", nrow(test)))
  up <- test[i, "start"]-test[test$chr == test[i, "chr"], ]$end
  if(sum(up>0)>0){
    test[i,"ub"] <- test[i, "start"] - min(up[up>0])
  }
  dw <- test[test$chr == test[i, "chr"], ]$start-test[i, "end"]
  if(sum(dw>0)>0){
    test[i,"db"] <- test[i, "end"] + min(dw[dw>0])
  }
}


test <- sj[,1:5]
test_start <- test[,-4]
test_start$type = "S"
test_end <- test[,-3]
test_end$type = "E"
colnames(test_start) <- c("sj", "chr", "pos", "strand", "type")
colnames(test_end) <- c("sj", "chr", "pos", "strand", "type")
test2 <- rbind(test_start,test_end)
test2 <- test2[order(test2$chr, test2$pos),]

test$ub <- NA
test$db <- NA

for(i in 1:nrow(test2)){
  print(paste(i, "of", nrow(test2)))
  if(test2[i, "type"] == "S"){
    for (j in i:1) {
      if(test2[j, "chr"] != test2[i, "chr"]){
        break()
      }
      if(test2[j, "type"]=="E" & test2[j, "pos"] < test2[i, "pos"]){
        test[test$sj == test2[i, "sj"],"ub"] <- test2[j, "pos"]
        break()
      }
    }
  }else{
    for(j in i:nrow(test2)){
      if(test2[j, "chr"] != test2[i, "chr"]){
        break()
      }
      if(test2[j, "type"]=="S" & test2[j, "pos"] > test2[i, "pos"]){
        test[test$sj == test2[i, "sj"],"db"] <- test2[j, "pos"]
        break()
      }
    }
  }
}

row.names(test) <- 1:nrow(test)

#head(test[test$chr==2,],40)
#c(which(!duplicated(test$chr)),
# which(!duplicated(test$chr))+1,
# which(!duplicated(test$chr))+2,
# which(!duplicated(test$chr))+3,
# which(!duplicated(test$chr))+4,
# which(!duplicated(test$chr))+5,
# which(!duplicated(test$chr))+6,
# which(!duplicated(test$chr))+7,
# which(!duplicated(test$chr))+8,
# which(!duplicated(test$chr))+9) -> n
#for(i in n){
# print(paste(i, "of", length(n)))
# if(!test[i, "ub"] %in% test[test$chr == test[i, "chr"],]$end){
#   test[i, "ub"] <- NA
# }
# if(!test[i, "db"] %in% test[test$chr == test[i, "chr"],]$start){
#   test[i, "db"] <- NA
# }
#}

test$ue <-  test$start - test$ub
test$de <- test$db - test$end
colnames(test)[6:9] <- c("upstream_boundary","downstream_boundary","upstream_exon_length","downstream_exon_length")

sj_11 <- merge(x = sj, y = test[,c(1,6:9)], by = "sj")
sj <- sj_11
save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj11.RData")

pdf("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/exon_length.pdf")
par(mfrow = c(3,1))
hist(with(sj,downstream_exon_length[classification=="Annotated"])[with(sj,downstream_exon_length[classification=="Annotated"])< 800],breaks = 100, main = "Annotated", xlab = "Length")
hist(with(sj,downstream_exon_length[classification=="Unannotated"])[with(sj,downstream_exon_length[classification=="Unannotated"])<800],breaks = 100, main = "Unannotated", xlab = "Length")
hist(with(sj,downstream_exon_length[classification=="Novel"])[with(sj,downstream_exon_length[classification=="Novel"])<800],breaks = 100, main = "Novel", xlab = "Length")
dev.off()



#===============================================================================================================================================#
#********************************************************************** Plot *******************************************************************#
#===============================================================================================================================================#


#### Novel splice junction result illustrate
## Wed Mar 28 2018
## By Tang Chao

#### Basic settings =============================================================================================================================
depar <- par()

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/")


col.a <- "#008080" # col.annotated
col.u <- "#EEE8AA" # col.unannotated
col.n <- "#CD853F" # col.novel

col.e <- "#00CCCC" # col.exon
col.i <- "#FF8C00" # col.intron


require(VennDiagram)

#### Load data ==================================================================================================================================
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/EST_mRNA_bulk_junction.RData")


#### 1. Basic statistic =========================================================================================================================


library(lattice)

#### Histogram of overhang ----
pdf(paste(fo, "Histogram of overhang.pdf", sep = ""))
histogram(~ overhang , data = sj, type = "density", 
          ylim = c(0,.18),
          panel = function(x, ...) {
            panel.histogram(x, ...)
            panel.densityplot(x, type = "density", adjust = 1.5)
          })
dev.off()



#### Histogram of GC content ----

pdf(paste(fo, "Hist of GC content.pdf", sep = ""))
par(mfrow = c(2,2), mar = c(5,4,1,2))
with(sj,hist(donor_exon_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(donor_intron_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(acceptor_exon_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(acceptor_intron_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC, adjust = 1.5), col = "#1E90FF"))
dev.off()


pdf(paste(fo, "Histogram of GC content.pdf", sep = ""))
GC_tab <- rbind(data.frame(GC = with(sj,donor_exon_GC), Host = "donor_exon"),
                data.frame(GC = with(sj,donor_intron_GC), Host = "donor_intron"),
                data.frame(GC = with(sj,acceptor_intron_GC), Host = "acceptor_intron"),
                data.frame(GC = with(sj,acceptor_exon_GC), Host = "acceptor_exon"))
GC_tab$Host <- factor(GC_tab$Host, levels = c("donor_exon", "donor_intron", "acceptor_exon", "acceptor_intron"))

histogram(~ GC|Host , data = GC_tab, xlab = "GC content", type = "density", 
          ylim = c(0,5), 
          panel = function(x, ...) {
            panel.histogram(x, ...)
            panel.densityplot(x, list(adjust = 1.5), list(type = "density"))
          }, layout = c(2, 2))

dev.off()


pdf(paste(fo, "trend of GC content in donor exon.pdf", sep = ""))

par(mfrow = c(4,2), mar = c(2,2,1,1))
with(sj,hist(donor_exon_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(donor_exon_GC[ReadSum>=10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[ReadSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10", bty = "n")

with(sj,hist(donor_exon_GC[CellSum >= 2], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[CellSum>=2], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 2", bty = "n")

with(sj,hist(donor_exon_GC[ReadSum >= 100], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[ReadSum>=100], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 100", bty = "n")

with(sj,hist(donor_exon_GC[CellSum >= 10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[CellSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 10", bty = "n")

with(sj,hist(donor_exon_GC[ReadSum >= 1000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[ReadSum>=1000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 1000", bty = "n")

with(sj,hist(donor_exon_GC[CellSum >= 50], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[CellSum>=50], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 50", bty = "n")

with(sj,hist(donor_exon_GC[ReadSum >= 10000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_exon_GC[ReadSum>=10000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10000", bty = "n")

dev.off()


pdf(paste(fo, "trend of GC content in donor intron.pdf", sep = ""))

par(mfrow = c(4,2), mar = c(2,2,1,1))
with(sj,hist(donor_intron_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(donor_intron_GC[ReadSum>=10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[ReadSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10", bty = "n")

with(sj,hist(donor_intron_GC[CellSum >= 2], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[CellSum>=2], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 2", bty = "n")

with(sj,hist(donor_intron_GC[ReadSum >= 100], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[ReadSum>=100], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 100", bty = "n")

with(sj,hist(donor_intron_GC[CellSum >= 10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[CellSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 10", bty = "n")

with(sj,hist(donor_intron_GC[ReadSum >= 1000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[ReadSum>=1000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 1000", bty = "n")

with(sj,hist(donor_intron_GC[CellSum >= 50], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[CellSum>=50], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 50", bty = "n")

with(sj,hist(donor_intron_GC[ReadSum >= 10000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(donor_intron_GC[ReadSum>=10000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10000", bty = "n")

dev.off()


pdf(paste(fo, "trend of GC content in acceptor intron.pdf", sep = ""))

par(mfrow = c(4,2), mar = c(2,2,1,1))
with(sj,hist(acceptor_intron_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(acceptor_intron_GC[ReadSum>=10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[ReadSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10", bty = "n")

with(sj,hist(acceptor_intron_GC[CellSum >= 2], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[CellSum>=2], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 2", bty = "n")

with(sj,hist(acceptor_intron_GC[ReadSum >= 100], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[ReadSum>=100], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 100", bty = "n")

with(sj,hist(acceptor_intron_GC[CellSum >= 10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[CellSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 10", bty = "n")

with(sj,hist(acceptor_intron_GC[ReadSum >= 1000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[ReadSum>=1000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 1000", bty = "n")

with(sj,hist(acceptor_intron_GC[CellSum >= 50], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[CellSum>=50], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 50", bty = "n")

with(sj,hist(acceptor_intron_GC[ReadSum >= 10000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_intron_GC[ReadSum>=10000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10000", bty = "n")

dev.off()


pdf(paste(fo, "trend of GC content in acceptor exon.pdf", sep = ""))

par(mfrow = c(4,2), mar = c(2,2,1,1))
with(sj,hist(acceptor_exon_GC, freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC, adjust = 1.5), col = "#1E90FF"))

with(sj,hist(acceptor_exon_GC[ReadSum>=10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[ReadSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10", bty = "n")

with(sj,hist(acceptor_exon_GC[CellSum >= 2], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[CellSum>=2], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 2", bty = "n")

with(sj,hist(acceptor_exon_GC[ReadSum >= 100], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[ReadSum>=100], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 100", bty = "n")

with(sj,hist(acceptor_exon_GC[CellSum >= 10], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[CellSum>=10], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 10", bty = "n")

with(sj,hist(acceptor_exon_GC[ReadSum >= 1000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[ReadSum>=1000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 1000", bty = "n")

with(sj,hist(acceptor_exon_GC[CellSum >= 50], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[CellSum>=50], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "CellSum >= 50", bty = "n")

with(sj,hist(acceptor_exon_GC[ReadSum >= 10000], freq = F, main = "", col = "#00FFFF"))
with(sj,lines(density(acceptor_exon_GC[ReadSum>=10000], adjust = 1.5), col = "#1E90FF"))
legend("topright", legend = "ReadSum >= 10000", bty = "n")

dev.off()


#### Histogram of sequence entropy ----


summary(sj$exon_Shannon_entropy)
hist(sj$exon_Shannon_entropy)
hist(sj$exon_Relative_entropy)

sum(sj$exon_Shannon_entropy < 1.8 | sj$intron_Shannon_entropy < 1.8)

pdf(paste(fo, "Histogram of sequence entropy.pdf", sep = ""))

entropy_tab <- rbind(data.frame(entropy = with(sj,exon_Shannon_entropy), region = "Exon", type = "Shannon Entropy"),
                     data.frame(entropy = with(sj,intron_Shannon_entropy), region = "Intron", type = "Shannon Entropy"),
                     data.frame(entropy = with(sj,exon_Relative_entropy), region = "Exon", type = "Relative Entropy"),
                     data.frame(entropy = with(sj,intron_Relative_entropy), region = "Intron", type = "Relative Entropy"))
# entropy_tab$Host <- factor(entropy_tab$Host, levels = c("donor_exon", "donor_intron", "acceptor_exon", "acceptor_intron"))

histogram(~ entropy|region*type , data = na.omit(entropy_tab), xlab = "Entropy", type = "density", 
          ylim = c(0,20), na.rm = T,
          panel = function(x, na.rm = T, ...) {
            panel.histogram(x, na.rm = T, ...)
            panel.densityplot(x, list(adjust = 1), list(type = "density"))
          }, layout = c(2, 2))

dev.off()
#

pdf(paste(fo, "Hist of sequence entropy.pdf", sep = ""))
par(mfrow = c(2,2), mar = c(5,4,1,2))
with(sj,hist(exon_Shannon_entropy, freq = F, main = "", ylim = c(0,15), breaks = 40, col = "#00FFFF"))
with(sj,lines(density(exon_Shannon_entropy), col = "#1E90FF"))

with(sj,hist(exon_Relative_entropy, freq = F, main = "", ylim = c(0,20), breaks = 40, col = "#00FFFF"))
with(sj,lines(density(exon_Relative_entropy[!is.na(exon_Relative_entropy)]), col = "#1E90FF"))

with(sj,hist(intron_Shannon_entropy, freq = F, main = "", ylim = c(0,15), breaks = 40, col = "#00FFFF"))
with(sj,lines(density(intron_Shannon_entropy), col = "#1E90FF"))

with(sj,hist(intron_Relative_entropy, freq = F, main = "", ylim = c(0,20), breaks = 40, col = "#00FFFF"))
with(sj,lines(density(intron_Relative_entropy[!is.na(intron_Relative_entropy)]), col = "#1E90FF"))
dev.off()


#### cutoff ----
sj_raw <- sj

sj <- sj[sj$exon_Shannon_entropy >= 1.8 & 
           !is.na(sj$exon_Relative_entropy) &
           sj$ReadSum >= 10 & 
           sj$overhang >= 10, ]
dim(sj)
# [1] 457282     37
sum(sj_raw$ReadSum>=10 & sj_raw$overhang>=10)
# [1] 477114
sum(sj_raw$ReadSum>=10 & sj_raw$overhang>=10 & sj_raw$CellSum>=2)
# [1] 247901


## Bar chart of annotation ----
pdf(paste(fo,"Bar chart of annotation.pdf", sep = ""))
x = barplot(table(sj$annotation)/nrow(sj)*100, names.arg = c("Unannotated","Annotated"), ylim = c(0,100), 
	main = "Annotation", ylab = "%", col = "white")
p <- round(table(sj$annotation)/nrow(sj)*100,2)
n <- table(sj$annotation)
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
dev.off()


## Boxplot of cells and reads ----
pdf(paste(fo,"Boxplot of cells and reads.pdf", sep = ""))
par(mfrow = c(1,2))
with(sj,boxplot(log10(CellSum[annotation==1]),log10(CellSum[annotation==0]),
                names = c("Annotated", "Unannotated"), main = "No. of cell sum", 
                ylab = "Log10(cells)"))

with(sj,boxplot(log10(ReadSum[annotation==1]),log10(ReadSum[annotation==0]),
                names = c("Annotated", "Unannotated"), main = "No. of reads sum", 
                ylab = "Log10(reads)"))
par(depar)
dev.off()


## Barchart of motif ----
pdf(paste(fo,"Barchart of motif.pdf", sep = ""), width = 10, height = 7)
par(mfrow = c(2,1))
par(mar = c(1,4,1,1))
x = with(sj, barplot(table(motif[annotation==1])/sum(table(motif[annotation==1]))*100, names.arg = F, ylim = c(0,60), 
	main = "SS motif", ylab = "%", col = "white"))
p <- round(with(sj, table(motif[annotation==1]))/sum(with(sj, table(motif[annotation==1])))*100,2)
n <- with(sj, table(motif[annotation==1]))
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
legend("topright", legend = "Annotated", bty = "n")

par(mar = c(3,4,1,1))
x = with(sj, barplot(table(motif[annotation==0])/sum(table(motif[annotation==0]))*100, ylim = c(0,42), ylab = "%", 
	col = "white", names.arg = c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT")))
p <- round(with(sj, table(motif[annotation==0]))/sum(with(sj, table(motif[annotation==0])))*100,2)
n <- with(sj, table(motif[annotation==0]))
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
legend("topright", legend = "Unannotated", bty = "n")
par(depar)
dev.off()

## Boxplot of junc distance
pdf(paste(fo,"Boxplot of junc distance.pdf", sep = ""))
d1 <- log10((as.numeric(sj[sj$annotation == 1, ]$end)+1) - (as.numeric(sj[sj$annotation == 1, ]$start)-1))
d2 <- log10((as.numeric(sj[sj$annotation == 0, ]$end)+1) - (as.numeric(sj[sj$annotation == 0, ]$start)-1))
boxplot(d1,d2, names = c("Annotated", "Unannotated"), main = "Distance", ylab = "log10(bp)")
dev.off()


#### Parse strand 


head(sj)
table(sj$strand)
# strand (0: undefined, 1: +, 2: -)
# intron motif: 0: non-canonical; 
# 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
table(sj[sj$strand==1,]$motif)
#      1      3      5 
# 170403  18151    510 
table(sj[sj$strand==2,]$motif)
#      2      4      6 
# 169028  18723    470 
table(sj[sj$strand==0,]$motif)
#     0 
# 79997 
sj -> sj_sp
sj_sp[sj_sp$motif==1,]$motif <- "GT/AG"
sj_sp[sj_sp$motif==2,]$motif <- "GT/AG"
sj_sp[sj_sp$motif==3,]$motif <- "GC/AG"
sj_sp[sj_sp$motif==4,]$motif <- "GC/AG"
sj_sp[sj_sp$motif==5,]$motif <- "AT/AC"
sj_sp[sj_sp$motif==6,]$motif <- "AT/AC"
sj_sp[sj_sp$motif==0,]$motif <- "non-canonical"

sj_sp$motif <- as.factor(sj_sp$motif)


## Barchart of motif Parse strand 
pdf(paste(fo,"Barchart of motif parsed.pdf", sep = ""), width = 10, height = 7)
par(mfrow = c(2,1))
par(mar = c(1,4,1,1))
x = with(sj_sp, barplot(table(motif[annotation==1])/sum(table(motif[annotation==1]))*100, names.arg = F, ylim = c(0,119), main = "SS motif", ylab = "%", col = "white"))
p <- round(with(sj_sp, table(motif[annotation==1]))/sum(with(sj_sp, table(motif[annotation==1])))*100,2)
n <- with(sj_sp, table(motif[annotation==1]))
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
legend("topright", legend = "Annotated", bty = "n")

par(mar = c(3,4,1,1))
x = with(sj_sp, barplot(table(motif[annotation==0])/sum(table(motif[annotation==0]))*100, ylim = c(0,100), ylab = "%", col = "white"))
p <- round(with(sj_sp, table(motif[annotation==0]))/sum(with(sj_sp, table(motif[annotation==0])))*100,2)
n <- with(sj_sp, table(motif[annotation==0]))
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
legend("topright", legend = "Unannotated", bty = "n")
par(depar)
dev.off()


#### 2. Basic statistic of novel splicing junction ==============================================================================================


mRNA_SJ <- row.names(mrna_junction)
EST_SJ <- row.names(est_junction)
bulk_SJ <- row.names(bulk_junction)
sc_SJ <- as.character(sj$sj)

pdf(paste(fo,"SJ_numble_of_different_splicing_methods.pdf",sep = ""), width = 10, height = 7)
x = barplot(c(length(sc_SJ),length(bulk_SJ),length(EST_SJ),length(mRNA_SJ)), 
            names.arg=c("SC","Bulk","EST","mRNA"), 
            main="SJ numble of different splicing methods", 
            ylim = c(0, 900000), col = "white")
text(x, c(length(sc_SJ),length(bulk_SJ),length(EST_SJ),length(mRNA_SJ)), labels = c(length(sc_SJ),length(bulk_SJ),length(EST_SJ),length(mRNA_SJ)), pos = 3)
dev.off()


library(VennDiagram)
venn_plot <- venn.diagram(list(mRNA = mRNA_SJ, EST = EST_SJ, Bulk = bulk_SJ, SC = sc_SJ), fill = rainbow(4), filename = NULL, main.cex = 2, main = "Venn plot of splicing junction in different sequencing methods", height = 6, width = 10)

pdf(paste(fo,"venn_plot_of_sc_bulk_mrna_EST_SJ.pdf",sep = ""), width = 10, height = 7)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


pdf(paste(fo, "Annotation of SC specific junction.pdf", sep = ""))
x = barplot(table(with(sj, annotation[SCSpe == "Y"])), main = "Annotation of SC specific junction",
            ylim = c(0,table(with(sj, annotation[SCSpe == "Y"]))[1]*1.2), col = "white")

p <- round(table(with(sj, annotation[SCSpe == "Y"]))/sum(table(with(sj, annotation[SCSpe == "Y"])))*100,2)
n <- table(with(sj, annotation[SCSpe == "Y"]))
text(x = x, y = n, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
dev.off()


pdf(paste(fo, "No. and percentage of SC SJ.pdf", sep = ""))
x = barplot(c(nrow(sj[sj$annotation==1, ]),nrow(sj[sj$SCSpe=="N" & sj$annotation==0, ]), nrow(sj[sj$SCSpe=="Y" & sj$annotation==0, ])),
            ylim = c(0, nrow(sj[sj$SCSpe=="Y" & sj$annotation==0, ])*1.2), col = c(col.a, col.u, col.n),
            main = "No. and percentage of SC SJ", names.arg = c("Annotated","Unannotated","Novel"))
n = c(nrow(sj[sj$annotation==1, ]),nrow(sj[sj$SCSpe=="N" & sj$annotation==0, ]), nrow(sj[sj$SCSpe=="Y" & sj$annotation==0, ]))
p = round(n/sum(n)*100,2)
text(x = x, y = n, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
dev.off()


dim(sj)
sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)

## Boxplot of cells and reads
pdf(paste(fo,"Boxplot of cells and reads of novel splice junction.pdf", sep = ""))
par(mfrow = c(1,2))
boxplot(log10(CellSum) ~ classification, data = sj, main = "No. of Cells", 
        ylab = "Log10(cells)", col = c(col.a, col.u, col.n), pch = 20)
boxplot(log10(ReadSum) ~ classification, data = sj, main = "No. of Reads", 
        ylab = "Log10(Reads)", col = c(col.a, col.u, col.n), pch = 20)
par(depar)
dev.off()

library(ggplot2)
## Vioplot of cells and reads
pdf(paste(fo,"Vioplot of cells and reads of novel splice junction.pdf", sep = ""))
ggplot(data = sj, aes(x = classification, y = log10(ReadSum), fill = classification))+
  geom_violin()+
  scale_fill_manual(values=c(col.a, col.u, col.n))+
  theme(legend.position="none")
       
ggplot(data = sj, aes(x = classification, y = log10(CellSum), fill = classification))+
  geom_violin()+
  scale_fill_manual(values=c(col.a, col.u, col.n))+
  theme(legend.position="none")

dev.off()


pdf(paste(fo,"Boxplot of novel junc distance.pdf", sep = ""))
d1 <- log10((as.numeric(sj[sj$classification == "Annotated", ]$end)+1) - (as.numeric(sj[sj$classification == "Annotated", ]$start)-1))
d2 <- log10((as.numeric(sj[sj$classification == "Unannotated", ]$end)+1) - (as.numeric(sj[sj$classification == "Unannotated", ]$start)-1))
d3 <- log10((as.numeric(sj[sj$classification == "Novel", ]$end)+1) - (as.numeric(sj[sj$classification == "Novel", ]$start)-1))

boxplot(d1,d2,d3, names = c("Annotated","Unannotated","Novel"), 
        main = "Distance", ylab = "log10(bp)", col = c(col.a, col.u, col.n))
dev.off()



#### 3. PhastCons conservation scores of annotated and novel SS =================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
sj_raw <- sj

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


write.table(sj, "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj_for_PhastCons.txt", sep = "\t", row.names = F, quote = F)


#### shell
cd /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/

## Prepare bed file:
awk '$7==1' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > Annotated_sj_start_site.bed
awk '$7==1' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > Annotated_sj_end_site.bed

awk '$7==0&&$11=="N"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > Uuannotated_sj_start_site.bed
awk '$7==0&&$11=="N"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > Unannotated_sj_end_site.bed

awk '$7==0&&$11=="Y"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > Novel_sj_start_site.bed
awk '$7==0&&$11=="Y"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > Novel_sj_end_site.bed

## Extract PhastCons conservation scores from bigwig file:
# Annotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_start_site_matrix.txt

# Annotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_end_site_matrix.txt


# Unannotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Uuannotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Unannotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Uuannotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Unannotated_sj_start_site_matrix.txt

# Unannotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Unannotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Unannotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Unannotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Unannotated_sj_end_site_matrix.txt


# Novel_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_start_site_matrix.txt

# Novel_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_end_site_matrix.txt


depar <- par()
require(ggplot2)
library(reshape2)


AEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
AEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
ASM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
ASA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)

UEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
UEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
USM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
USA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)

NEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
NEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
NSM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
NSA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_start_site_agg.txt", stringsAsFactors = F, header = T)


## PhastCons conservation scores agg:
pdf(paste(fo, "PhastCons conservation scores agg.pdf",sep = ""), width = 10, height = 7)
#### Plot in one figure
rbind(ASA, USA, NSA) -> test
test$Type <- factor(rep(c("Annotated","Unannotated","Novel"), c(nrow(ASA),nrow(USA),nrow(NSA))),
                    levels = c("Annotated","Unannotated","Novel"))

rbind(AEA, UEA, NEA) -> test2
test2$Type <- factor(rep(c("Annotated","Unannotated","Novel"), c(nrow(AEA),nrow(UEA),nrow(NEA))),
                     levels = c("Annotated","Unannotated","Novel"))

rbind(test2,test) -> test3
test3$Role <- factor(rep(c("Acceptor","Donor"),c(nrow(test2),nrow(test))),
                     levels = c("Acceptor","Donor"))

ggplot(test3, aes(x = Position, y = Mean_1, colour = Type)) + 
  geom_path(size = 1.5)+
  facet_grid(. ~ Role)+
  scale_color_manual(values = c(col.a, col.u, col.n))

dev.off()


## PhastCons conservation scores matrix:
library(reshape2)
p1 <- ggplot(melt(ASM), aes(x = variable, y = value, fill = col.a))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank())+
  scale_fill_manual(values = col.a)

p2 <- ggplot(melt(AEM), aes(x = variable, y = value, fill = col.a))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())+
  scale_fill_manual(values = col.a)

p3 <- ggplot(melt(USM), aes(x = variable, y = value, fill = col.u))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank())+
  scale_fill_manual(values = col.u)

p4 <- ggplot(melt(UEM), aes(x = variable, y = value, fill = col.u))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())+
  scale_fill_manual(values = col.u)

p5 <- ggplot(melt(NSM), aes(x = variable, y = value, fill = col.n))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.y = element_blank())+
  scale_fill_manual(values = col.n)

p6 <- ggplot(melt(NEM), aes(x = variable, y = value, fill = col.n))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none")+
  scale_fill_manual(values = col.n)

library(cowplot)
pdf(paste(fo, "PhastCons conservation scores matrix.pdf",sep = ""), width = 10, height = 7)
plot_grid(p2, p1, p4, p3, p6, p5, ncol = 2)
dev.off()


#### 4. Seqlogo of splice sites =================================================================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
sj <- sj[sj$overhang>=10 & sj$ReadSum >= 10,]
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

## plot seqLogo plot:
seqLogoFromMatrix <- function(tableResult, plot = TRUE, ic.scale=FALSE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15){
  ## 
  N = ncol(tableResult)
  
  if("seqLogo" %in% installed.packages()){
    library(seqLogo)
  }else{
    source("https://bioconductor.org/biocLite.R")
    biocLite("seqLogo")
    library(seqLogo)
  }
  ## 
  if(is.matrix(tableResult)){
    data_tu <- matrix(rep(NA,4*N),nrow = 4)
    rownames(data_tu) <- c("A", "C", "G", "T")
    
    for(j in 1:N){
      for(i in c("A", "C", "G", "T")){
        data_tu[i,j] <- sum(tableResult[,j] == i)
      }
    }
    data_tu_p <- data_tu/colSums(data_tu)
    if(plot == FALSE){
      return(data_tu_p)
    }else{
      pwm <- makePWM(data_tu_p)
      
      seqLogo(pwm, ic.scale = ic.scale, xaxis = xaxis, yaxis = yaxis, xfontsize = xfontsize, yfontsize = yfontsize)
    }
  }else{
    stop("The input data of tableResult must be a matrix")
  }
}


pdf(paste(fo,"Seqlogo.pdf",sep=""), width = 10, height = 5)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 48, 56), split="")), xaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 0 & SCSpe=="N"]), 48, 56), split="")), xaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 48, 56), split="")), xaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 145, 153), split="")), xaxis = F, yaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 0 & SCSpe=="N"]), 145, 153), split="")), xaxis = F, yaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 145, 153), split="")), xaxis = F, yaxis = F)

dev.off()


## function:
MyEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log2(x))
  Reduce("+", -tmp)
}

round(apply(annotated_donor_percentage, 2, MyEntropy),2)


MyRelativeEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log(x/0.25))
  Reduce("+", tmp)
}


annotated_donor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 48, 56), split="")), plot = F)
unannotated_donor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 0 & SCSpe=="N"]), 48, 56), split="")), plot = F)
nover_donor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 48, 56), split="")), plot = F)
annotated_acceptor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 145, 153), split="")), plot = F)
unannotated_acceptor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 0 & SCSpe=="N"]), 145, 153), split="")), plot = F)
nover_acceptor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 145, 153), split="")), plot = F)


round(apply(annotated_donor, 2, MyEntropy),2)
# [1] 1.88 1.50 0.96 0.03 0.10 1.27 1.38 1.14 1.81
round(apply(unannotated_donor, 2, MyEntropy),2)
# [1] 1.92 1.58 1.20 0.37 0.53 1.52 1.61 1.40 1.84
round(apply(nover_donor, 2, MyEntropy),2)
# [1] 1.92 1.69 1.56 1.16 1.52 1.97 1.95 1.86 1.99
round(apply(annotated_acceptor, 2, MyEntropy),2)
# [1] 1.58 1.55 1.98 1.25 0.01 0.03 1.78 1.94 1.99
round(apply(unannotated_acceptor, 2, MyEntropy),2)
# [1] 1.75 1.74 1.99 1.59 0.40 0.36 1.85 1.96 1.99
round(apply(nover_acceptor, 2, MyEntropy),2)
# [1] 1.98 1.96 2.00 1.87 1.23 1.17 1.69 1.91 1.99


round(apply(annotated_donor, 2, MyRelativeEntropy),2)
# [1] 0.08 0.35 0.72 1.37 1.32 0.51 0.43 0.60 0.13
round(apply(unannotated_donor, 2, MyRelativeEntropy),2)
# [1] 0.05 0.29 0.55 1.13 1.02 0.33 0.27 0.42 0.11
round(apply(nover_donor, 2, MyRelativeEntropy),2)
# [1] 0.05 0.21 0.30 0.58 0.33 0.02 0.04 0.10 0.01
round(apply(annotated_acceptor, 2, MyRelativeEntropy),2)
# [1] 0.29 0.31 0.01 0.52 1.38 1.37 0.15 0.04 0.00
round(apply(unannotated_acceptor, 2, MyRelativeEntropy),2)
# [1] 0.17 0.18 0.01 0.29 1.11 1.13 0.10 0.03 0.01
round(apply(nover_acceptor, 2, MyRelativeEntropy),2)
# [1] 0.01 0.03 0.00 0.09 0.54 0.57 0.22 0.06 0.01


#### 5. MaxEntScan for Maximum Entropy of splicing sites ========================================================================================


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


library(MASS)
library(ggplot2)
library(viridis)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

annotated <- sj[sj$annotation == 1, c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(annotated) <- c("Donor","Acceptor")
annotated$density <- get_density(annotated$Donor, annotated$Acceptor)

unannotated <- sj[sj$annotation == 0 & sj$SCSpe == "N", c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(unannotated) <- c("Donor","Acceptor")
unannotated$density <- get_density(unannotated$Donor, unannotated$Acceptor)

novel <- sj[sj$novel == "Y", c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(novel) <- c("Donor","Acceptor")
novel$density <- get_density(novel$Donor, novel$Acceptor)


p1 <- ggplot(annotated, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")
p2 <- ggplot(unannotated, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")
p3 <- ggplot(novel, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")
library(cowplot)
pdf(paste(fo, "Maximum Entropy plot of novel and other sj.pdf", sep = ""), width = 12, height = 4)
plot_grid(p1,p2,p3, ncol = 3)
dev.off()



#### 6. GC Content ==============================================================================================================================


p1 <- ggplot(sj, aes(x = classification, y = donor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor exon")+
  labs(y = "GC Content")

p2 <- ggplot(sj, aes(x = classification, y = acceptor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n))+
  theme(legend.position = "none",
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  ggtitle("Acceptor exon")


p3 <- ggplot(sj, aes(x = classification, y = donor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor intron")+
  labs(y = "GC Content")


p4 <- ggplot(sj, aes(x = classification, y = acceptor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n))+
  theme(legend.position = "none",
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  ggtitle("Acceptor intron")

pdf(paste(fo,"Vioplot of GC content.pdf", sep = ""))
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()



#### 7. Stop Codon ==============================================================================================================================


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


pdf(paste(fo,"Boxplot of stop codons.pdf", sep = ""))
par(mfrow = c(1,2), las = 3)
### Donor
boxplot(donor_exon_sum_stop_codon ~ classification, data = sj,
        boxwex = 0.25, at = 1:3 - 0.2, ylim = c(0,8.3),
        col = c(col.a, col.u, col.n), outline = F,
        main = "Donor", names = c("Exon", "Exon", "Exon"),
        ylab = "Sum of Stop Codon")
boxplot(donor_intron_sum_stop_codon ~ classification, data = sj, add = TRUE,
        boxwex = 0.25, at = 1:3 + 0.2, outline = F, names = c("Intron", "Intron", "Intron"),
        col = c(col.a, col.u, col.n))
### Acceptor
boxplot(acceptor_exon_sum_stop_codon ~ classification, data = sj,
        boxwex = 0.25, at = 1:3 - 0.2, ylim = c(0,8.3),
        col = c(col.a, col.u, col.n), outline = F,
        main = "Acceptor", names = c("Exon", "Exon", "Exon"),
        ylab = "Sum of Stop Codon")
boxplot(acceptor_intron_sum_stop_codon ~ classification, data = sj, add = TRUE,
        boxwex = 0.25, at = 1:3 + 0.2, outline = F, names = c("Intron", "Intron", "Intron"),
        col = c(col.a, col.u, col.n))

dev.off()



#### 8. Heat map ================================================================================================================================


depar <- par()

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/")


col.a <- "#008080" # col.annotated
col.u <- "#EEE8AA" # col.unannotated
col.n <- "#CD853F" # col.novel

col.e <- "#00CCCC" # col.exon
col.i <- "#FF8C00" # col.intron


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

library(pheatmap)
pdf(paste(fo,"Heat map of SJ.pdf", sep = ""))
pheatmap(data[,row.names(annotation3)],annotation_col = annotation3, cluster_cols = F,
         show_colnames = F, show_rownames = F)
dev.off()


#### 9. Motif enrichment ========================================================================================================================


depar <- par()

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/")


col.a <- "#008080" # col.annotated
col.u <- "#EEE8AA" # col.unannotated
col.n <- "#CD853F" # col.novel

col.e <- "#00CCCC" # col.exon
col.i <- "#FF8C00" # col.intron


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


## We can use this function to read a meme format ifle(download from http://meme-suite.org/db/motifs) and output their PWM format file
## We just need to input our meme file path and file name.
## Usage: pwms <- ReadMEME2PWM("/Users/tangchao/Downloads/motif_databases/RNA/Ray2013_rbp_Homo_sapiens.meme")
ReadMEME2PWM <- function(path){
  test <- read.table(path, fill = T)
  remain <- c(grep(as.character(test$V1), pattern =  "MOTIF"), grep(as.character(test$V1), pattern = "[[:digit:]]"))
  ## We just need the first column is MOTIF and numeric elements
  test2 <- test[sort(remain),1:4]
  test3 <- split.data.frame(test2, rep(1:length(grep(test2[,1], pattern = "MOTIF")), diff(c(grep(test2[,1], pattern = "MOTIF"),nrow(test2)+1))))
  ## Split data.frame by MOTIF
  pwms <- test3
  for(i in 1:length(pwms)){
    names(pwms)[i] <- paste(as.character(pwms[[i]][1,2]), as.character(pwms[[i]][1,3]), sep = "_")
    tmp <- as.matrix(data.frame(pwms[[i]][-1,]))
    colnames(tmp) <- c("A","C","G","T")
    row.names(tmp) <- 1:nrow(tmp)
    mode(tmp) <- "numeric"
    pwms[[i]] <- t(log2(tmp/0.25))
  }
  return(pwms)
}


#### **** enrichTest **** ============================================================================================================================================================
##########################################################
## Over-Representation of a Motif in a Set of Sequences ##
##########################################################
## Author: Thomas Girke
## Last update: http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/enrichMotif.R
## Utility: calculates enrichment/depletion p-values of motif matches in sequence 
## sets relative to frequency in search space (e.g. all promoters of genome)

## Variables for: phyper(k, D, n-D, N, lower.tail=FALSE)
  ## k: number of sample sequences with at least one motif match (number of good items in sequence sample)  
  ## D: number of database sequence with at least one motif match (number of good items in sequence database)
  ## n: number of sequences in database minus D (number of bad items in sequence database)
  ## N: number of sample sequences (number of drawn items) 

## Input data sets/arguments
  ## pwm: position weight matrix
  ## cutoff: minimum score cutoff for matching PWM
  ## set: DNAStringSet containing sample sequences
  ## db: DNAStringSet containing sequences of search space

## Motif enrichment/depletion function
enrichTest <- function(pwm=pwm, revcomp=TRUE, cutoff="90%", set=set, db=db, occurrence=1, lower.tail=FALSE) { 
# pwm: query position weight matrix; 
# revcomp: search only one or both strands; 
# cutoff: min.score for matching; 
# set: test set; 
# db: background sequence data set; 
# occurrence: min number of matches in each sequence; 
# low.tail: see phyper fct.
  if(length(db) < length(set)) stop("The size of 'db' needs to be >= 'set'.")
  ## Obtain values for required variables and calculate hypergeometric p-values
  if(revcomp==TRUE) {
    krclog <- sapply(set, function(x) countPWM(reverseComplement(pwm), x, min.score=cutoff)) >= occurrence
    k <- sum(klog | krclog)

    Drclog <- sapply(db, function(x) countPWM(reverseComplement(pwm), x, min.score=cutoff)) >= occurrence
    D <- sum(Dlog | Drclog)
  }else{
  	klog <- sapply(set, function(x) countPWM(pwm, x, min.score=cutoff)) >= occurrence
  	k <- sum(klog)

  	Dlog <- sapply(db, function(x) countPWM(pwm, x, min.score=cutoff)) >= occurrence
  	D <- sum(Dlog)
  }

  n <- length(db) - D
  N <- length(set)
  pval <- phyper(k-1, D, n, N, lower.tail=lower.tail)
  ## Bonferroni p-value adjustment
  if(k==0) {
    adj_pval <- pval
  } else {
    adj_pval <- pval * k
  }
  adj_pval[adj_pval>1] <- 1
  enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}


pwms <- ReadMEME2PWM("/mnt/data1/projects/UCB/results_ucb/tangchao/motif/Ray2013_rbp_Homo_sapiens.meme")


#### **** MotifSearch **** ===========================================================================================================================================================
## We can use this function to search the given motif PWM in target sequence
## Usage: motif_search_intron <- mclapply(set_intron, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=10)
## In this example, set_intron is a list of FASTA
MotifSearch <- function(fasta = x, pwms = pwms){
  motif_count<-NULL
  for (seed in names(pwms)) {
    #print(seed)
    cmd<-paste("pwms$'",seed,"'",sep="")
    pwm<-eval(parse(text=cmd))
    #print(pwmb)
    motif_count_to_add<-countPWM(pwm, DNAString(fasta), min.score="90%")
    names(motif_count_to_add)<-seed
    motif_count<-c(motif_count,motif_count_to_add)
  }
  return(motif_count)
}


Novel_donor_exon <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,1,50))
names(Novel_donor_exon) <- sj[sj$classification=="Novel",]$sj

Novel_donor_intron <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,51,100))
names(Novel_donor_intron) <- sj[sj$classification=="Novel",]$sj

Novel_acceptor_intron <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,101,150))
names(Novel_acceptor_intron) <- sj[sj$classification=="Novel",]$sj

Novel_acceptor_exon <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,151,200))
names(Novel_acceptor_exon) <- sj[sj$classification=="Novel",]$sj


library(parallel)
library(BSgenome)
Novel_donor_exon_motif_search<-mclapply(Novel_donor_exon, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
Novel_donor_intron_motif_search<-mclapply(Novel_donor_intron, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
Novel_acceptor_intron_motif_search<-mclapply(Novel_acceptor_intron, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
Novel_acceptor_exon_motif_search<-mclapply(Novel_acceptor_exon, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
## 
save(Novel_donor_exon_motif_search, Novel_donor_intron_motif_search, Novel_acceptor_intron_motif_search, Novel_acceptor_exon_motif_search,
	file = paste(ro, "motif_search.RData", sep = ""))



donor_exon <- do.call(rbind,Novel_donor_exon_motif_search)
donor_exon_sd <- apply(donor_exon,1,sd)
donor_exon_sd_motif <- apply(donor_exon,2,sd)

donor_intron <- do.call(rbind,Novel_donor_intron_motif_search)
donor_intron_sd <- apply(donor_intron,1,sd)
donor_intron_sd_motif <- apply(donor_intron,2,sd)

acceptor_intron <- do.call(rbind,Novel_acceptor_intron_motif_search)
acceptor_intron_sd <- apply(acceptor_intron,1,sd)
acceptor_intron_sd_motif <- apply(acceptor_intron,2,sd)

acceptor_exon <- do.call(rbind,Novel_acceptor_exon_motif_search)
acceptor_exon_sd <- apply(acceptor_exon,1,sd)
acceptor_exon_sd_motif <- apply(acceptor_exon,2,sd)

library(pheatmap)
plot_list=list()
plot_list[[1]] <- pheatmap(donor_exon[donor_exon_sd >= quantile(donor_exon_sd, probs = .95),], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F)[[4]]
plot_list[[2]] <- pheatmap(donor_intron[donor_intron_sd >= quantile(donor_intron_sd, probs = .95),], show_rownames = F, show_colnames = F,cluster_rows = F, cluster_cols = F)[[4]]
plot_list[[3]] <- pheatmap(acceptor_intron[acceptor_intron_sd >= quantile(acceptor_intron_sd, probs = .95),], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F)[[4]]
plot_list[[4]] <- pheatmap(acceptor_exon[acceptor_exon_sd >= quantile(acceptor_exon_sd, probs = .95),], show_rownames = F, show_colnames = F,cluster_rows = F, cluster_cols = F)[[4]]
library(gridExtra)
g<-do.call(grid.arrange, plot_list)
#g2<-do.call(function(x) {grid.arrange(x, ncol = 1)}, plot_list)
library(ggplot2)
ggsave(paste(fo, "MEME_motif_count.pdf", sep = ""),g)

pdf(paste(fo, "MEME_motif_counts.pdf", sep = ""))
pheatmap(donor_exon[donor_exon_sd >= quantile(donor_exon_sd, probs = .95),], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F)
pheatmap(donor_intron[donor_intron_sd >= quantile(donor_intron_sd, probs = .95),], show_rownames = F, show_colnames = F,cluster_rows = F, cluster_cols = F)
pheatmap(acceptor_intron[acceptor_intron_sd >= quantile(acceptor_intron_sd, probs = .95),], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F)
pheatmap(acceptor_exon[acceptor_exon_sd >= quantile(acceptor_exon_sd, probs = .95),], show_rownames = F, show_colnames = F,cluster_rows = F, cluster_cols = F)
dev.off()


#### enrichment

library(parallel)
library(Biostrings)
db <- as.list(substr(sj$ParsedFasta,1,50))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,1,50))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment

enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 10)

depletion_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = TRUE)
						  }, mc.cores = 10)

donor_exon_enrich <- do.call(rbind,enrich_matrix)
donor_exon_depletion <- do.call(rbind,depletion_matrix)

save(donor_exon_enrich, donor_exon_depletion, file = paste(ro,"donor_exon_motif_enrich.RData",sep = ""))


## motif enrichment and depletion of novel SJ donor intron

db <- as.list(substr(sj$ParsedFasta,51,100))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,51,100))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment
enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 10)

depletion_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = TRUE)
						  }, mc.cores = 10)

donor_intron_enrich <- do.call(rbind,enrich_matrix)
donor_intron_depletion <- do.call(rbind,depletion_matrix)

save(donor_intron_enrich, donor_intron_depletion, file = paste(ro,"donor_intron_motif_enrich.RData",sep = ""))


## motif enrichment and depletion of novel SJ accptor intron

db <- as.list(substr(sj$ParsedFasta,101,150))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,101,150))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment
enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 10)

depletion_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = TRUE)
						  }, mc.cores = 10)

accptor_intron_enrich <- do.call(rbind,enrich_matrix)
accptor_intron_depletion <- do.call(rbind,depletion_matrix)

save(accptor_intron_enrich, accptor_intron_depletion, file = paste(ro,"accptor_intron_motif_enrich.RData",sep = ""))


## motif enrichment and depletion of novel SJ accptor exon

db <- as.list(substr(sj$ParsedFasta,151,200))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,151,200))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment
enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 10)

depletion_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1,lower.tail = TRUE)
						  }, mc.cores = 10)

accptor_exon_enrich <- do.call(rbind,enrich_matrix)
accptor_exon_depletion <- do.call(rbind,depletion_matrix)

save(accptor_intron_enrich, accptor_intron_depletion, file = paste(ro,"accptor_intron_motif_enrich.RData",sep = ""))


enrich <- data.frame(DonorExon = donor_exon_enrich[,"adj_pval"],
           DonorIntron = donor_intron_enrich[,"adj_pval"],
           AcceptorIntron = accptor_intron_enrich[,"adj_pval"],
           AcceptorExon = accptor_exon_enrich[,"adj_pval"])

-log10(enrich) -> logenrich
logenrich[,1][is.infinite(logenrich[,1])] <- 300
logenrich[,2][is.infinite(logenrich[,2])] <- 300
logenrich[,3][is.infinite(logenrich[,3])] <- 300
logenrich[,4][is.infinite(logenrich[,4])] <- 300

pdf(paste(fo,"MEME motif enrichment.pdf", sep = ""))
pheatmap(logenrich, show_rownames = F)
dev.off()


#### RBPDB ======================================================================================================================================

pfm_list <- list.files("/mnt/data5/BGI/UCB/tangchao/data/RBPDB/matrices_human/PFMDir", pattern = "pfm", full.names = T)
matrix_list <- read.table("/mnt/data5/BGI/UCB/tangchao/data/RBPDB/PFMDir/matrix_list.txt", row.names = 1, stringsAsFactors = F)

pwms <- as.list(1:length(pfm_list))
names(pwms) <- matrix_list[do.call(rbind,strsplit(pfm_list, split = "[/\\.]"))[, 11],]$V3

library(Biostrings)

for(i in 1:length(pfm_list)){
  print(i)
  tmp <- as.matrix(read.table(pfm_list[i], row.names = c("A","C","G","T")))
  if(mean(round(colSums(tmp)))>1){
  	pwms[[i]] <- log2(apply(tmp,2,function(x) x/sum(x))/0.25)
  }else{
  	pwms[[i]] <- log2(tmp/.25)
  }  
}


Novel_donor_exon <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,1,50))
names(Novel_donor_exon) <- sj[sj$classification=="Novel",]$sj

Novel_donor_intron <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,51,100))
names(Novel_donor_intron) <- sj[sj$classification=="Novel",]$sj

Novel_acceptor_intron <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,101,150))
names(Novel_acceptor_intron) <- sj[sj$classification=="Novel",]$sj

Novel_acceptor_exon <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,151,200))
names(Novel_acceptor_exon) <- sj[sj$classification=="Novel",]$sj


library(parallel)
library(BSgenome)
Novel_donor_exon_motif_search<-mclapply(Novel_donor_exon, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
Novel_donor_intron_motif_search<-mclapply(Novel_donor_intron, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
Novel_acceptor_intron_motif_search<-mclapply(Novel_acceptor_intron, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
Novel_acceptor_exon_motif_search<-mclapply(Novel_acceptor_exon, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=30)
## 
save(Novel_donor_exon_motif_search, Novel_donor_intron_motif_search, Novel_acceptor_intron_motif_search, Novel_acceptor_exon_motif_search,
	file = paste(ro, "motif_search_RBPDB.RData", sep = ""))


#### enrichment

library(parallel)
library(Biostrings)
db <- as.list(substr(sj$ParsedFasta,1,50))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,1,50))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment

enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff="90%", set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 10)


donor_exon_enrich <- do.call(rbind,enrich_matrix)

save(donor_exon_enrich, file = paste(ro,"donor_exon_motif_enrich_RBPDB_PFM_90%.RData",sep = ""))


## motif enrichment and depletion of novel SJ donor intron

db <- as.list(substr(sj$ParsedFasta,51,100))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,51,100))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment
enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff="90%", set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 10)


donor_intron_enrich <- do.call(rbind,enrich_matrix)

save(donor_intron_enrich, file = paste(ro,"donor_intron_motif_enrich_RBPDB_PFM_90%.RData",sep = ""))


## motif enrichment and depletion of novel SJ accptor intron

db <- as.list(substr(sj$ParsedFasta,101,150))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,101,150))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment
enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff="90%", set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 20)


accptor_intron_enrich <- do.call(rbind,enrich_matrix)

save(accptor_intron_enrich, file = paste(ro,"accptor_intron_motif_enrich_RBPDB_PFM_90%.RData",sep = ""))


## motif enrichment and depletion of novel SJ accptor exon

db <- as.list(substr(sj$ParsedFasta,151,200))
names(db) <- sj$sj
db <- lapply(db, DNAString)

set <- as.list(substr(sj[sj$classification=="Novel",]$ParsedFasta,151,200))
names(set) <- sj[sj$classification=="Novel",]$sj
set <- lapply(set, DNAString)


## enrichment
enrich_matrix <- mclapply(pwms, function(x) {
							enrichTest(pwm=x, revcomp=FALSE, cutoff="90%", set=set, db=db, occurrence=1,lower.tail = FALSE)
						  }, mc.cores = 20)


accptor_exon_enrich <- do.call(rbind,enrich_matrix)

save(accptor_intron_enrich, file = paste(ro,"accptor_intron_motif_enrich_PFM_90%.RData",sep = ""))


enrich <- data.frame(DonorExon = donor_exon_enrich[,"adj_pval"],
           DonorIntron = donor_intron_enrich[,"adj_pval"],
           AcceptorIntron = accptor_intron_enrich[,"adj_pval"],
           AcceptorExon = accptor_exon_enrich[,"adj_pval"])

as.matrix(-log10(enrich)) -> logenrich

row.names(logenrich) <- paste(row.names(accptor_exon_enrich),row.names(matrix_list[do.call(rbind,strsplit(pfm_list, split = "[/\\.]"))[, 11],]), sep = "|")
logenrich <- logenrich[rowSums(logenrich)>0,]

logenrich[,1][is.infinite(logenrich[,1])] <- 300
logenrich[,2][is.infinite(logenrich[,2])] <- 300
logenrich[,3][is.infinite(logenrich[,3])] <- 300
logenrich[,4][is.infinite(logenrich[,4])] <- 300

logenrich[,1][logenrich[,1]>200] <- 200
logenrich[,2][logenrich[,2]>200] <- 200
logenrich[,3][logenrich[,3]>200] <- 200
logenrich[,4][logenrich[,4]>200] <- 200

donor_exon_enrich[donor_exon_enrich[,"adj_pval"]<1,]
donor_intron_enrich[donor_intron_enrich[,"adj_pval"]<1,]
accptor_intron_enrich[accptor_intron_enrich[,"adj_pval"]<1,]
accptor_exon_enrich[accptor_exon_enrich[,"adj_pval"]<1,]

library(pheatmap)
pdf(paste(fo,"RBPDB motif PFM enrichment_2.pdf", sep = ""), width = 7, height = 10)
pheatmap(logenrich, show_rownames = T)
dev.off()











