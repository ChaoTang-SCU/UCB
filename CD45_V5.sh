
cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '1,400p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '401,800p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '801,1200p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '1201,1600p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '1601,2000p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '2001,2400p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '2401,2800p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done

cd /mnt/data5/BGI/UCB/tangchao/data/STAR
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt | sed -n '2801,3039p'`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f
bedtools intersect -a /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result/$f/$f.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
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

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result")

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_merge[cov_merge<10] <- 0

cov_merge[,row.names(tsne)] -> cov_tu

library(ggplot2)
t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45_RABC_v5.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()



#### Only plot the major isoform
cov_merge[,row.names(tsne)] -> cov_tu

for(i in 1:ncol(cov_tu)){
	cov_tu[,i][-which.max(cov_tu[,i])] <- 0
}


t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45_RABC_v5_2.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()




cov_merge[,row.names(tsne)] -> cov_tu

table(apply(cov_tu,2,function(x) sum(x>0)))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/Bar chart of cell isoform numbers.pdf")
x = barplot(table(apply(cov_tu,2,function(x) sum(x>0))), ylim = c(0,1200),
	ylab = "No. Cells", col = "white")
n <- table(apply(cov_tu,2,function(x) sum(x>0)))
text(x = x, y = n, labels = n, cex = .8, pos = 3)
dev.off()


#### Only plot the cells have only one isoform
for(i in 1:ncol(cov_tu)){
	if(sum(cov_tu[,i] > 0 )>1){
		cov_tu[,i] <- 0
	}	
}


t(t(cov_tu)[order(t(cov_tu)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab1)[1])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = 1)+
  scale_colour_gradient(high = 'Grey',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

t(t(cov_tu)[order(t(cov_tu)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = 1)+
  scale_colour_gradient(high = 'red',low = 'Grey')+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.title = element_blank(),
  		axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45_RABC_v5_3.pdf", height = 6,width = 16)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, nrow = 2)
dev.off()









cd /mnt/nfs_nas/EGAD00001002671/EGAD00001002671_bam_decrypted
for f in `ls *mRNA*bam | head -n 10`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result_bulk/Tcell/${f:0:15}
bedtools intersect -a /mnt/nfs_nas/EGAD00001002671/EGAD00001002671_bam_decrypted/$f -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45_e3_7_hg19.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result_bulk/Tcell/${f:0:15}/${f:0:15}.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/$f*
done



cd /mnt/nfs_nas/EGAD00001002674/EGAD00001002674_bam_decrypted
for f in `ls *bam | head -n 10`
do 
mkdir /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result_bulk/Monocytes/${f:0:15}
bedtools intersect -a /mnt/nfs_nas/EGAD00001002674/EGAD00001002674_bam_decrypted/$f -b /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/CD45_e3_7_hg19.bed > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.bam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools fastq -n -t /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.bam > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.fq
hisat2 -x /mnt/data4/reference/GRCh38/grch38_snp_tran/genome_snp_tran --dta -U /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.fq > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sam
/mnt/data7/tangchao/biosoft/samtools-1.8/samtools view -bS /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sam | /mnt/data7/tangchao/biosoft/samtools-1.8/samtools sort - > /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sorted.bam
/home/zhaoyuancun/bin/stringtie-1.3.4c/stringtie -G /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/custom.gff -e -B /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}.sorted.bam -o /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result_bulk/Monocytes/${f:0:15}/${f:0:15}.txt
rm /mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/tmp/${f:0:15}*
done




















