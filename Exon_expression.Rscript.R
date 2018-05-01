## We can use this script to calculate the median(M), mad(M), coverage(C) of given region(s) after bedtools intersect
args <- commandArgs(T)
## args1: The path and name of bedtools intersect result
## args2: The path and name of results

## Usage: Rscript --vanilla

library(data.table)

MMCLGenomeRegion = function(intron){
  print(intron)
  tmp<-fread(intron,header=FALSE)
  setnames(tmp,c("chr","start","end","chr.b","start.b","end.b","count"))
  setkey(tmp,"chr","start","end")
  tmp_sortby<-tmp[,,by="chr,start,end"]
  class(tmp_sortby$count)<-"double"

  #print(paste("step 1 of 4: calculate median"))
  tmp_median <- tmp_sortby[,median(rep(count,(end.b-start.b))),by="chr,start,end"]

  #print(paste("step 2 of 4: calculate mad"))
  tmp_mad <- tmp_sortby[,round(mad(rep(count,(end.b-start.b))),1),by="chr,start,end"]
  #tmp_median_mad <- tmp_sortby[,cbind(median(rep(count,(1+end.b-start.b))),mad(rep(count,(1+end.b-start.b)))),by="chr,start,end"]
  #coverage

  #print(paste("step 3 of 4: calculate coverage"))
  tmp_cov <- tmp_sortby[,round(length(rep(count,(end.b-start.b))[which(rep(count,(1+end.b-start.b))>0)])/length(rep(count,(1+end.b-start.b)))*100,1),by="chr,start,end"]

  #print(paste("step 4 of 4: calculate intron retention level"))

  output<-cbind(tmp_median, tmp_mad$V1, tmp_cov$V1)
  setnames(output, 4:6, c("median","mad","coverage"))

  return(output)
  #example, add "" to the myfilename
  #ts<-medianGenomeRegion("myfilename")
}

te <- MMCLGenomeRegion(args[1])

write.table(te, args[2], sep="\t", quote=F, row.names=F)
