
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
dim(sj)

test <- sj[,1:5]
test$ub <- NA
test$db <- NA
for (i in 1:nrow(test)) {
	print(paste(i, "of", nrow(test)))
	up <- test[i, "start"]-test[test$chr == test[i, "chr"], ]$end
	if(sum(up>0)>0){
		test[i,"ub"] <- test[i, "start"] - min(up[up>0])
	}
	dw <- test[test$chr == test[i, "chr"], ]$start-test[i, "end"]
	if(sum(dw>0)>0){
		test[i,"db"] <- test[i, "end"] + min(dw[dw>0])
	}
}


test <- sj[,1:5]
test_start <- test[,-4]
test_start$type = "S"
test_end <- test[,-3]
test_end$type = "E"
colnames(test_start) <- c("sj", "chr", "pos", "strand", "type")
colnames(test_end) <- c("sj", "chr", "pos", "strand", "type")
test2 <- rbind(test_start,test_end)
test2 <- test2[order(test2$chr, test2$pos),]

test$ub <- NA
test$db <- NA

for(i in 1:nrow(test2)){
  print(paste(i, "of", nrow(test2)))
  if(test2[i, "type"] == "S"){
    for (j in i:1) {
      if(test2[j, "chr"] != test2[i, "chr"]){
        break()
      }
      if(test2[j, "type"]=="E" & test2[j, "pos"] < test2[i, "pos"]){
        test[test$sj == test2[i, "sj"],"ub"] <- test2[j, "pos"]
        break()
      }
    }
  }else{
    for(j in i:nrow(test2)){
      if(test2[j, "chr"] != test2[i, "chr"]){
        break()
      }
      if(test2[j, "type"]=="S" & test2[j, "pos"] > test2[i, "pos"]){
        test[test$sj == test2[i, "sj"],"db"] <- test2[j, "pos"]
        break()
      }
    }
  }
}

row.names(test) <- 1:nrow(test)

#head(test[test$chr==2,],40)
#c(which(!duplicated(test$chr)),
#	which(!duplicated(test$chr))+1,
#	which(!duplicated(test$chr))+2,
#	which(!duplicated(test$chr))+3,
#	which(!duplicated(test$chr))+4,
#	which(!duplicated(test$chr))+5,
#	which(!duplicated(test$chr))+6,
#	which(!duplicated(test$chr))+7,
#	which(!duplicated(test$chr))+8,
#	which(!duplicated(test$chr))+9) -> n
#for(i in n){
#	print(paste(i, "of", length(n)))
#	if(!test[i, "ub"] %in% test[test$chr == test[i, "chr"],]$end){
#		test[i, "ub"] <- NA
#	}
#	if(!test[i, "db"] %in% test[test$chr == test[i, "chr"],]$start){
#		test[i, "db"] <- NA
#	}
#}

test$ue <-  test$start - test$ub
test$de <- test$db - test$end
colnames(test)[6:9] <- c("upstream_boundary","downstream_boundary","upstream_exon_length","downstream_exon_length")

sj_11 <- merge(x = sj, y = test[,c(1,6:9)], by = "sj")
sj <- sj_11
save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj11.RData")


sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)
  Annotated Unannotated       Novel
     176529       79987      613530


summary(with(sj,upstream_exon_length[classification=="Annotated"]))
summary(with(sj,upstream_exon_length[classification=="Unannotated"]))
summary(with(sj,upstream_exon_length[classification=="Novel"]))


summary(with(sj,downstream_exon_length[classification=="Annotated"]))
summary(with(sj,downstream_exon_length[classification=="Unannotated"]))
summary(with(sj,downstream_exon_length[classification=="Novel"]))


par(mfrow = c(3,1))
hist(with(sj,downstream_exon_length[classification=="Annotated"])[with(sj,downstream_exon_length[classification=="Annotated"])< 800],breaks = 100, main = "Annotated", xlab = "Length")
hist(with(sj,downstream_exon_length[classification=="Unannotated"])[with(sj,downstream_exon_length[classification=="Unannotated"])<800],breaks = 100, main = "Unannotated", xlab = "Length")
hist(with(sj,downstream_exon_length[classification=="Novel"])[with(sj,downstream_exon_length[classification=="Novel"])<800],breaks = 100, main = "Novel", xlab = "Length")





#### FASTA of cell type specific novel SJ

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/cell_type_specific_unannotated_SJ.RData")
load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj11.RData")

novel_sj_maxreads <- merge(x = novel_sj_maxreads, y = sj[,c("sj","chr","start","end","upstream_boundary","downstream_boundary","upstream_exon_length","downstream_exon_length")], by.x = "SJ", by.y = "sj", all.x = T)

write.table(novel_sj_maxreads, "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/cell_type_specific_unannotated_SJ.txt", sep = "\t", row.names=F, quote=F)




















