#### Call genome DNA sequence around SS
## Sat Mar 25 2018
## By Tang Chao

#!/bin/sh

#  samtools.sh
#  We can get the sequence of a given region of the human genome(FASTA)
#
#  Created by TangChao on Tue Jan  9 15:06:12 CST 2018.
#  
# left flank exons FASTA
for f in `awk '{print $2":"$3-50"-"$3+49}' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | sed -n '2,$p'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> /mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/fasta/start_site.fa
done

sed ':a;N;$!ba;s/\n/ /g' start_site.fa | sed 's/>/\n/g' | awk -F " " '{print $1"\t"$2$3}' > start_site.txt
sed -i '1d' start_site.txt



# right flank exons FASTA
for f in `awk '{print $2":"$4-49"-"$4+50}' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | sed -n '2,$p'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> /mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/fasta/end_site.fa
done

sed ':a;N;$!ba;s/\n/ /g' end_site.fa | sed 's/>/\n/g' | awk -F " " '{print $1"\t"$2$3}' > end_site.txt
sed -i '1d' end_site.txt



paste start_site.txt end_site.txt | awk '{print $1"\t"$3"\t"$2$4}' > SS.fa.txt




#### Merge with sj information


## R
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.RData")
fa <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/fasta/SS.fa.txt", sep = "\t", header = F, stringsAsFactors = F)

sj$fa <- as.character(fa$V3)

sj$ParsedFasta <- as.character(fa$V3)

library(Biostrings)
## Get the Reverse complementary chain

for (i in 1:nrow(sj)){
	if(sj[i, "strand"] == 2){
		sj[i, "ParsedFasta"] <- as.character(reverseComplement(DNAStringSet(sj[i, "fa"])))
	}
} 

save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/sj_parsed_fastq.RData")
write.table(sj, "/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/sj_parsed_fastq.txt", sep = "\t", row.names = F, quote = F)



fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/seqlogo/")


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

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 48, 56), split="")))
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 48, 56), split="")))
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 145, 153), split="")))
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 145, 153), split="")))

dev.off()

colSums(test)



test <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 48, 56), split="")), plot=FALSE)











