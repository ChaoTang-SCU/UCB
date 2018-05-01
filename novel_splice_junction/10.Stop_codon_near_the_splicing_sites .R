#### Stop codon near the splicing sites 
## Mon Mar 26 2018
## By: Tang Chao


#### 0. Basic settings ======================================================================


depar <- par()
setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/stopcodon")


## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/stopcodon/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/stopcodon/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/stopcodon/")


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



#### 1. Load data ===========================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/sj_parsed_fastq.RData")


stopcodon <- stopCoden(fasta = substr(sj$ParsedFasta,1,50), strand = rep(1,length(sj$ParsedFasta)))

remain = grepl(sj$ParsedFasta, pattern = "N") == FALSE
sum(!remain)





