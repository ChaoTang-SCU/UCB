
#### GMAP for introns:
## Sat Dec 29 2017
## By Tang Chao

#### GMAP

#### UCSC

#### gmap build
gmap_build -D /mnt/data1/reference/GMAPDB/UCSC/human -d GRCh38 -C /mnt/data1/reference/UCSC/human/GRCh38/hg38.fa.gz -g

#### mrna
## splicesite
gmap -D /mnt/data1/reference/GMAPDB/UCSC/human/ -d GRCh38 -B 5 -t 8 -w 300000 /mnt/data1/reference/UCSC/human/GRCh38/mrna.fa -n 1 -z sense_filter -f 6  > mRNA_splicesite_2017_12_29.txt
## introns
gmap -D /mnt/data1/reference/GMAPDB/UCSC/human/ -d GRCh38 -B 5 -t 8 -w 300000 /mnt/data1/reference/UCSC/human/GRCh38/mrna.fa -n 1 -z sense_filter -f introns  > mRNA_introns_2017_12_29.txt



#### EST
## splicesite
gmap -D /mnt/data1/reference/GMAPDB/UCSC/human/ -d GRCh38 -B 5 -t 8 -w 300000 /mnt/data1/reference/UCSC/human/GRCh38/est.fa -n 1 -z sense_filter -f 6  > est_splicesite_2017_12_29.txt
## introns
gmap -D /mnt/data1/reference/GMAPDB/UCSC/human/ -d GRCh38 -B 5 -t 8 -w 300000 /mnt/data1/reference/UCSC/human/GRCh38/est.fa -n 1 -z sense_filter -f introns  > est_introns_2017_12_29.txt


awk -F " " '{print $2}' mRNA_introns_2017_12_29.txt | sort | uniq > temp
awk -F "[\.:]" '{print $1"\t"$2"\t"$4}' temp > mRNA_intron.bed
rm temp
awk '{print $1":"$2+1"-"$3-1}' mRNA_intron.bed > mRNA_intron.sj

awk -F " " '{print $2}' est_introns_2017_12_29.txt | sort | uniq > temp
awk -F "[\.:]" '{print $1"\t"$2"\t"$4}' temp > est_intron.bed
rm temp
awk '{print $1":"$2+1"-"$3-1}' est_intron.bed > est_intron.sj