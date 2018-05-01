#### MaxEntScan for Maximum Entropy of splicing sites
## Mon Mar 26 2018
## By: Tang Chao

## Prepare FASTA
## R
setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy")

load("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/sj_parsed_fastq.RData")
with(sj,ParsedFasta[annotation==1])
head(substr(with(sj,ParsedFasta[novel=="N" & strand!=0]),48,56))
donor <- substr(with(sj,ParsedFasta[novel=="N" & strand!=0]),48,56)
acceptor <- substr(with(sj,ParsedFasta[novel=="N" & strand!=0]),131,153)
other_for_MaxEnt <- as.data.frame(cbind(donor,acceptor), stringsAsFactors=F)

donor <- substr(with(sj,ParsedFasta[novel=="Y" & strand!=0]),48,56)
acceptor <- substr(with(sj,ParsedFasta[novel=="Y" & strand!=0]),131,153)
novel_for_MaxEnt <- as.data.frame(cbind(donor,acceptor), stringsAsFactors=F)

write.table(other_for_MaxEnt$donor, "other_donor.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(other_for_MaxEnt$acceptor, "other_acceptor.txt", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(novel_for_MaxEnt$donor, "novel_donor.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(novel_for_MaxEnt$acceptor, "novel_acceptor.txt", row.names = F, col.names = F, quote = F, sep = "\t")




#### MaxEntScan for Maximum Entropy:



cd /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/MaxEntScan/
perl score3.pl /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/novel_acceptor.txt > /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/novel_acceptor_MaxEnt.txt

perl score3.pl /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/other_acceptor.txt > /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/other_acceptor_MaxEnt.txt

perl score5.pl /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/novel_donor.txt > /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/novel_donor_MaxEnt.txt

perl score5.pl /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/other_donor.txt > /mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/other_donor_MaxEnt.txt





#### plot:

fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy/Figure/")

setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Maximum_Entropy")

NDME <- read.table("novel_donor_MaxEnt.txt", sep = "\t", header = F)
NAME <- read.table("novel_acceptor_MaxEnt.txt", sep = "\t", header = F)

ODME <- read.table("other_donor_MaxEnt.txt", sep = "\t", header = F)
OAME <- read.table("other_acceptor_MaxEnt.txt", sep = "\t", header = F)


other <- subset.data.frame(ODME, select = 2)
colnames(other) <- "donor"
other$acceptor <- OAME$V2

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

other$density <- get_density(other$acceptor, other$donor)


novel <- subset.data.frame(NDME, select = 2)
colnames(novel) <- "donor"
novel$acceptor <- NAME$V2

novel$density <- get_density(novel$acceptor, novel$donor)


pdf(paste(fo, "Maximum Entropy plot of novel and other sj.pdf"))
ggplot(other, aes(x = acceptor, y = donor,color = density))+
  geom_point(size = .6)+
  scale_color_viridis()

ggplot(novel, aes(x = acceptor, y = donor,color = density))+
  geom_point(size = .6)+
  scale_color_viridis()

dev.off()






