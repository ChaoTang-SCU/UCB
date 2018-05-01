gtf <- read.table("/mnt/data5/BGI/UCB/tangchao/Homo_sapiens.GRCh38.87.transccript.exon.gtf", header = FALSE, sep = "\t", stringsAsFactors = F)

test <- as.data.frame(do.call(rbind, strsplit(gtf$V9, split = "[; ]")))


valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}



gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})), stringsAsFactors=F)

gtf2 <- gtf

gtf2$V9 <- gsub(pattern = "transcript_id", x = gtf2$V9, replacement = "Transcript_ID")

transcript_id <- paste("transcript_id ",'"',  gtf_infor$transcript_id, "_", gtf_infor$gene_name, "_", gtf_infor$gene_id,'"',  "; ", sep = "")

head(paste(gtf2$V9, transcript_id))

gtf$V9 <- paste(gtf2$V9, transcript_id)

write.table(gtf, "/mnt/data5/BGI/UCB/tangchao/MY.GRCh38.87.for.ggsashimi.gtf", col.names = F, row.names = F, quote = F, sep = "\t")


