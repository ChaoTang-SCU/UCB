#### Use R package BCRANK to find annotated regions' motifs
#### Then find the motifs in novel regions

#### Firstly, we should extract the FASTA sequence of the target area:
#### shell
cd /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table

#### Donor 

for f in `awk -F "\t" '$38=="Annotated"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$3-300"-"$3-1} else {print $2":"$4+1"-"$4+300}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Annotated_Donor_exon300.fa
done


for f in `awk -F "\t" '$38=="Annotated"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$3"-"$3+299} else {print $2":"$4-299"-"$4}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Annotated_Donor_intron300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$3-300"-"$3-1} else {print $2":"$4+1"-"$4+300}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Novel_Donor_exon300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$3"-"$3+299} else {print $2":"$4-299"-"$4}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Novel_Donor_intron300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0&&$9>1' sj_tu.txt | awk '{if($5==1) {print $2":"$3-300"-"$3-1} else {print $2":"$4+1"-"$4+300}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./HFNovel_Donor_exon300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0&&$9>1' sj_tu.txt | awk '{if($5==1) {print $2":"$3"-"$3+299} else {print $2":"$4-299"-"$4}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./HFNovel_Donor_intron300.fa
done




#### Acceptor

for f in `awk -F "\t" '$38=="Annotated"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$4-299"-"$4} else {print $2":"$3"-"$3+299}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Annotated_Acceptor_intron300.fa
done


for f in `awk -F "\t" '$38=="Annotated"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$4+1"-"$4+300} else {print $2":"$3-300"-"$3-1}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Annotated_Acceptor_exon300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$4-299"-"$4} else {print $2":"$3"-"$3+299}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Novel_Acceptor_intron300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0' sj_tu.txt | awk '{if($5==1) {print $2":"$4+1"-"$4+300} else {print $2":"$3-300"-"$3-1}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./Novel_Acceptor_exon300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0&&$9>1' sj_tu.txt | awk '{if($5==1) {print $2":"$4-299"-"$4} else {print $2":"$3"-"$3+299}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./HFNovel_Acceptor_intron300.fa
done


for f in `awk -F "\t" '$38=="Novel"&&$5!=0&&$9>1' sj_tu.txt | awk '{if($5==1) {print $2":"$4+1"-"$4+300} else {print $2":"$3-300"-"$3-1}}'`
do
	samtools faidx /mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa $f >> ./HFNovel_Acceptor_exon300.fa
done














library(BCRANK)
Annotated_Donor_exon_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_Donor_exon300.fa")
set.seed(0)
Annotated_Donor_exon_BCRANKout <- bcrank(Annotated_Donor_exon_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(Annotated_Donor_exon_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/Annotated_Donor_exon_BCRANKout.RData")


HFNovel_Donor_exon_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_Donor_exon300.fa")
set.seed(0)
HFNovel_Donor_exon_BCRANKout <- bcrank(HFNovel_Donor_exon_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(HFNovel_Donor_exon_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/HFNovel_Donor_exon_BCRANKout.RData")


Annotated_Donor_intron_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_Donor_intron300.fa")
set.seed(0)
Annotated_Donor_intron_BCRANKout <- bcrank(Annotated_Donor_intron_fastaFile, restarts=20, use.P1=TRUE, use.P2=TRUE)
save(Annotated_Donor_intron_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/Annotated_Donor_intron_BCRANKout20.RData")


HFNovel_Donor_intron_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_Donor_intron300.fa")
set.seed(0)
HFNovel_Donor_intron_BCRANKout <- bcrank(HFNovel_Donor_intron_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(HFNovel_Donor_intron_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/HFNovel_Donor_intron_BCRANKout.RData")


library(BCRANK)
HFNovel_Acceptor_intron_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_Acceptor_intron300.fa")
set.seed(0)
HFNovel_Acceptor_intron_BCRANKout <- bcrank(HFNovel_Acceptor_intron_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(HFNovel_Acceptor_intron_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/HFNovel_Acceptor_intron_BCRANKout.RData")



library(BCRANK)
HFNovel_Acceptor_exon_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_Acceptor_exon300.fa")
set.seed(0)
HFNovel_Acceptor_exon_BCRANKout <- bcrank(HFNovel_Acceptor_exon_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(HFNovel_Acceptor_exon_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/HFNovel_Acceptor_exon_BCRANKout.RData")


library(BCRANK)
Annotated_Acceptor_intron_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_Acceptor_intron300.fa")
set.seed(0)
Annotated_Acceptor_intron_BCRANKout <- bcrank(Annotated_Acceptor_intron_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(Annotated_Acceptor_intron_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/Annotated_Acceptor_intron_BCRANKout.RData")


library(BCRANK)
Annotated_Acceptor_exon_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_Acceptor_exon300.fa")
set.seed(0)
Annotated_Acceptor_exon_BCRANKout <- bcrank(Annotated_Acceptor_exon_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(Annotated_Acceptor_exon_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/Annotated_Acceptor_exon_BCRANKout.RData")



library(BCRANK)
Annotated_Donor_intron_fastaFile <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_Donor_intron300.fa")
set.seed(0)
Annotated_Donor_intron_BCRANKout <- bcrank(Annotated_Donor_intron_fastaFile, restarts=100, use.P1=TRUE, use.P2=TRUE)
save(Annotated_Donor_intron_BCRANKout, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/Annotated_Donor_intron_BCRANKout.RData")







load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/Annotated_Donor_exon_BCRANKout.RData")
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/BCRANK/HFNovel_Donor_exon_BCRANKout.RData")
library(BCRANK)



data(Annotated_Donor_exon_BCRANKout)
topMotif <- toptable(Annotated_Donor_exon_BCRANKout, 1)
weightMatrix <- pwm(topMotif, normalize = FALSE)
weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
library(seqLogo)
seqLogo(weightMatrixNormalized)
plot(topMotif)


