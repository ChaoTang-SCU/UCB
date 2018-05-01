#### Call PhastCons conservation scores of annotated and novel SS
## Sat Mar 25 2018
## By Tang Chao

cd /mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons

## Prepare bed file:
awk '$7==1' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' | sed -n '2,$p' > Annotated_sj_start_site.bed
awk '$7==1' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' | sed -n '2,$p' > Annotated_sj_end_site.bed

awk '$11=="Y"' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' | sed -n '2,$p' > Novel_sj_start_site.bed
awk '$11=="Y"' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' | sed -n '2,$p' > Novel_sj_end_site.bed

## Extract PhastCons conservation scores from bigwig file:
# Annotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_start_site_matrix.txt

# Annotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_end_site_matrix.txt

# Novel_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_start_site_matrix.txt

# Novel_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_end_site_matrix.txt




#### Cutoff: 10 reads


cd /mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffReads10

## Prepare bed file:
awk '$7==1&&$8>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' | sed -n '2,$p' > Annotated_sj_start_site.bed
awk '$7==1&&$8>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' | sed -n '2,$p' > Annotated_sj_end_site.bed

awk '$11=="Y"&&$8>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' | sed -n '2,$p' > Novel_sj_start_site.bed
awk '$11=="Y"&&$8>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' | sed -n '2,$p' > Novel_sj_end_site.bed

## Extract PhastCons conservation scores from bigwig file:
# Annotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_start_site_matrix.txt

# Annotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_end_site_matrix.txt

# Novel_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_start_site_matrix.txt

# Novel_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_end_site_matrix.txt




#### Cutoff: 10 cells


cd /mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffCells10

## Prepare bed file:
awk '$7==1&&$9>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' | sed -n '2,$p' > Annotated_sj_start_site.bed
awk '$7==1&&$9>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' | sed -n '2,$p' > Annotated_sj_end_site.bed

awk '$11=="Y"&&$9>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' | sed -n '2,$p' > Novel_sj_start_site.bed
awk '$11=="Y"&&$9>=10' /mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/NSS.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' | sed -n '2,$p' > Novel_sj_end_site.bed

## Extract PhastCons conservation scores from bigwig file:
# Annotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_start_site_matrix.txt

# Annotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_end_site_matrix.txt

# Novel_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_start_site_matrix.txt

# Novel_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_end_site_matrix.txt


