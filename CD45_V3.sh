

## The stringtie can't achieve the goal that assembly CDC45 transcript.
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/UCB3.01421/UCB3.01421_Aligned.sortedByCoord.out.bam \
				   -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed \
				   > UCB3.01421.bam

bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/UCB3.01421/UCB3.01421_Aligned.sortedByCoord.out.bam \
				   -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed \
				   > UCB3.01421.bam



for f in UCB1.00003 UCB1.00006 UCB1.00008 UCB1.00016 UCB4.01529 UCB4.01516 UCB4.01522
do
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam \
				   -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed \
				   > $f.bam
done

for f in UCB1.00003 UCB1.00006 UCB1.00008 UCB1.00016 UCB4.01529 UCB4.01516 UCB4.01522
do
samtools index $f.bam
done

for f in UCB1.00003 UCB1.00006 UCB1.00008 UCB1.00016 UCB4.01529 UCB4.01516 UCB4.01522
do
/mnt/data4/software/IGVTools/igvtools count -z 5 -w 25 $f.bam $f.bam.tdf /mnt/data1/reference/ensembl/human/GRCh38.91/GRCh38.chrom.sizes
done


/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G RA.gtf \
												 -e \
												 -B \
												 -o ./stout_RA/UCB3.01421.gtf \
												 /mnt/data5/BGI/UCB/tangchao/data/STAR/UCB3.01421/UCB3.01421_Aligned.sortedByCoord.out.bam


samtools fastq -n -t -c 9 UCB3.01421.bam > UCB3.01421.fq.gz


/usr/local/bin/hisat2-2.1.0/hisat2 \
				-x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran \
				--dta \
				-U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/UCB3.01421.fq \
				-S UCB3.01421.sam

samtools view -bS UCB3.01421.sam | samtools sort - file_sorted

/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G custom.gff -e file_sorted.bam | less


bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/UCB3.01421/UCB3.01421_Aligned.sortedByCoord.out.bam \
				-b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed |\
				samtools fastq -n -t - | \
				hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran \
				--dta \
				-U - |\
				samtools view -bS - |\
				samtools sort - file_sorted



bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/UCB3.01421/UCB3.01421_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t - | hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U - | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS - | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - file_sorted2

bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/UCB3.01421/UCB3.01421_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.hisat.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.hisat.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.hisat.file_sorted2
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/custom.gff -e /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.hisat.file_sorted2 > test.txt

/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/bedtools/test.hisat.file_sorted2 -o test


cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `ls `
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/tmp/$f*
done



#R
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

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/result")

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_merge[,row.names(tsne)] -> cov_tu

library(ggplot2)
t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey',guide = F)+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/CD45_RABC_v3_2.pdf", height = 12,width = 6)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, ncol = 2)
dev.off()


















## So, I have devolope ourself methold:
# Filstly, we call the expression of E3~E7 using bedtools like IR

# E3
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/Exon37.bed -b /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/genomecov/UCB1.01421_genomecov -wa -wb > UCB1.01421.Exon37.cov

Rscript --vanilla Exon_expression.Rscript UCB1.01421.Exon37.cov UCB1.01421.Exon37.exp.txt


#### For loop for all cells
cd /mnt/data5/BGI/UCB/Star
for f in $(ls )
do
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/Exon37.bed -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_cov/${f}_cov
done

for f in $(cat /mnt/data1/projects/UCB/results_ucb/tangchao/IR/cellInfo.txt)
do
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/Exon37.bed -b /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_cov/${f}_cov
done

cd /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_cov
for f in $(ls )
do
Rscript --vanilla /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v3/Exon_expression.Rscript $f /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_exp/${f:0:10}.txt
done


# 
library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "median",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_exp/")
system.time(median_merge<-multmerge(path))
median_merge <- as.data.frame(median_merge)
colnames(median_merge) <- substr(colnames(median_merge),1,10)

colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB3", replacement="UCB4")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB1", replacement="UCB3")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB5", replacement="UCB1")

row.names(median_merge) <- paste("E", 3:7, sep = "")
median_merge <- t(median_merge)

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

library(ggplot2)
f1 <- ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(median_merge[row.names(tsne),1]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(colnames(median_merge)[1])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

f2 <- ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(median_merge[row.names(tsne),2]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(colnames(median_merge)[2])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

f3 <- ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(median_merge[row.names(tsne),3]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(colnames(median_merge)[3])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

f4 <- ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(median_merge[row.names(tsne),4]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(colnames(median_merge)[4])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

f5 <- ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(median_merge[row.names(tsne),5]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(colnames(median_merge)[5])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Exon_Exp.pdf", height = 5,width = 10)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, ncol = 3)
dev.off()



load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
# 1	198638671	198757283

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

loci <- do.call(rbind,strsplit(row.names(psi_sj_same_start_table), split="[:-]"))
loci <- data.frame(loci, stringsAsFactors = F)
loci$X2 <- as.numeric(loci$X2)
loci$X3 <- as.numeric(loci$X3)

sum(loci$X1 == 1 & loci$X2 >= 198638671 & loci$X3 <= 198757283)
# 150
psi <- psi_sj_same_start_table[loci[,1] == 1 & loci[,2] >= 198638671 & loci[,3] <= 198757283, row.names(tsne)]
#psi[is.na(psi)] <- -0.5
dim(psi)
# [1]  150 3039

library(ggplot2)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Intron_centric_PSI.pdf", height = 5,width = 10)
for (i in 1:nrow(psi)){
	print(ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  		geom_point(aes(colour = as.numeric(psi[i,row.names(tsne)])), size = .6)+
  		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
  		ggtitle(row.names(psi)[i])+
  		theme(panel.background = element_blank()))
}#statements

dev.off()


library(ggplot2)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Intron_centric_PSI2.pdf", height = 5,width = 10)
for (i in 1:nrow(psi)){
	tmp_data <- data.frame(tsne, PSI = as.numeric(psi[i,row.names(tsne)]))
	tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
	print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
  		geom_point(aes(colour = PSI), size = 2)+
  		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
  		ggtitle(row.names(psi)[i])+
  		theme(panel.background = element_blank())
)
}#statements
dev.off()




pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/CD45_Intron_centric_PSI3.pdf", height = 5,width = 10)
for (i in 1:nrow(psi)){
	tmp_data <- data.frame(tsne, PSI = as.numeric(psi[i,row.names(tsne)]))
	tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
	tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
	tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
	print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
  		geom_point(aes(colour = PSI), size = 2)+
  		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
  		ggtitle(row.names(psi)[i])+
  		theme(panel.background = element_blank())
)
}#statements
dev.off()




	



