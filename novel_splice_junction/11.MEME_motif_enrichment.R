#### MEME motif 
## Mon Mar 26 2018
## By: Tang Chao

load("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/sj_parsed_fastq.RData")


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
enrichTest <- function(pwm=pwm, revcomp=TRUE, cutoff=0.6, set=set, db=db, occurrence=1, lower.tail=FALSE) { 
# pwm: query position weight matrix; 
# revcomp: search only one or both strands; 
# cutoff: min.score for matching; 
# set: test set; 
# db: background sequence data set; 
# occurrence: min number of matches in each sequence; 
# low.tail: see phyper fct.
  if(length(db) < length(set)) stop("The size of 'db' needs to be >= 'set'.")
  ## Obtain values for required variables and calculate hypergeometric p-values
  klog <- sapply(set, function(x) countPWM(pwm, x, min.score=cutoff)) >= occurrence
  k <- sum(klog)
  if(revcomp==TRUE) {
    krclog <- sapply(set, function(x) countPWM(reverseComplement(pwm), x, min.score=cutoff)) >= occurrence
    k <- sum(klog | krclog)
  }
  Dlog <- sapply(db, function(x) countPWM(pwm, x, min.score=cutoff)) >= occurrence
  D <- sum(Dlog)
  if(revcomp==TRUE) {
    Drclog <- sapply(db, function(x) countPWM(reverseComplement(pwm), x, min.score=cutoff)) >= occurrence
    D <- sum(Dlog | Drclog)
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


## motif count of novel SJ donor exon

head(sj)

set <- as.list(substr(sj[sj$annotation==1,]$ParsedFasta,1,50))
names(set) <- sj[sj$annotation==1,]$sj

library(BSgenome)
motif_search<-mclapply(set, function(x) {
	#to_search<-x
	motif_count<-NULL
	for (seed in names(pwms)) {
		#print(seed)
		cmd<-paste("pwms$'",seed,"'",sep="")
		pwm<-eval(parse(text=cmd))
		#print(pwmb)
		motif_count_to_add<-countPWM(pwm, x, min.score="90%")
		names(motif_count_to_add)<-seed
		motif_count<-c(motif_count,motif_count_to_add)
	}
	return(motif_count)
},mc.cores=1)


motif_search_bind <- do.call(rbind,motif_search)
motif_search_bind_sd <- apply(motif_search_bind,1,sd)
motif_search_bind_sd_motif <- apply(motif_search_bind,2,sd)


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


library(parallel)
motif_search<-mclapply(set, function(x) MotifSearch(fasta = x, pwms = pwms), mc.cores=10)
## This way is better than before one.


motif_search_bind <- do.call(rbind,motif_search)
motif_search_bind_sd <- apply(motif_search_bind,1,sd)
motif_search_bind_sd_motif <- apply(motif_search_bind,2,sd)


library(pheatmap)
pheatmap(log10(as.matrix(motif_search_intron_bind)+1),show_colnames = F, show_rownames = F)
pheatmap(as.matrix(motif_search_intron_bind[motif_search_intron_bindd_sd<1,motif_search_intron_bind_sd_motif<=.4]),show_colnames = T, show_rownames = F)


## motif seqlogo ##
seqLogo((2^(pwms$RNCMPT00004_BRUNOL4)*.25), xaxis = F)
grid.text("RNCMPT00004_BRUNOL4", x = .6, y = .9)


## enrichment
enrich_matrix <- as.list(rep(NA,length(pwms)))
names(enrich_matrix) <- names(pwms)

for(i in 1:length(pwms)){
	seed <- names(pwms)[i]
	cmd<-paste("pwms$'",seed,"'",sep="")
	pwm<-eval(parse(text=cmd))

	enrich<-enrichTest(pwm=pwm, revcomp=FALSE, cutoff=0.9, set=set, db=db, occurrence=1)	
	enrich_matrix[[i]]<-enrich
}

enrich_matrix_bind <- do.call(rbind,enrich_matrix)




