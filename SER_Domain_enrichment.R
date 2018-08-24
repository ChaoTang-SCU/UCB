
#### CDD ====================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/CDD.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/CDD.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$CDD == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/CDD.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)

#### Gene3D =================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Gene3D.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/Gene3D.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$Gene3D.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/gene3D.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### HAMAP ==================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/HAMAP.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/HAMAP.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$HAMAP.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/HAMAP.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### hmmpanther =============================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/hmmpanther.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/hmmpanther.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$hmmpanther.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/hmmpanther.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### Interpro ===============================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Interpro.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/Interpro.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$Interpro.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/Interpro.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### pfam ===================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/pfam.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/pfam.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$Pfam.domain.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/pfam.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### PIRST ==================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/PIRST.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/PIRST.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$PIRSF.domain.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/PIRST.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### Prints =================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Prints.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/Prints.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$Prints.domain.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/Prints.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### PROSITE_pattern ========================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/PROSITE.pattern.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/PROSITE.pattern.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$PROSITE.patterns.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/PROSITE.pattern.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### PROSITE.profile ========================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/PROSITE.profile.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/PROSITE.profile.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$PROSITE.profiles.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/PROSITE.profile.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### SFLD ===================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/SFLD.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/SFLD.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$SFLD == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/SFLD.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### SMART ==================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/SMART.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/SMART.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$SMART.domains.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/SMART.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### Superfamily ============================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Superfamily.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/Superfamily.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$Superfamily.domains.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/Superfamily.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)


#### TIGRFAM ================================

cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/TIGRFAM.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/TIGRFAM.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$Superfamily.domains.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

pdf("/Users/tangchao/Desktop/Domain/Result/TIGRFAM.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)







#### Cell specific SE ===========================================================================================================================
load("~/Desktop/cell_sep_SE.RData")
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))

Bcell_sep_range <- mapply(function(x) {paste(x[1], ":", as.numeric(x[3])+1, "-", as.numeric(x[4])-1, sep = "")}, strsplit(sj_Bcell, "[:-]"))
Tcell_sep_range <- mapply(function(x) {paste(x[1], ":", as.numeric(x[3])+1, "-", as.numeric(x[4])-1, sep = "")}, strsplit(sj_Tcell, "[:-]"))
MK_sep_range <- mapply(function(x) {paste(x[1], ":", as.numeric(x[3])+1, "-", as.numeric(x[4])-1, sep = "")}, strsplit(sj_MK, "[:-]"))
NK_sep_range <- mapply(function(x) {paste(x[1], ":", as.numeric(x[3])+1, "-", as.numeric(x[4])-1, sep = "")}, strsplit(sj_NK, "[:-]"))
Mono_sep_range <- mapply(function(x) {paste(x[1], ":", as.numeric(x[3])+1, "-", as.numeric(x[4])-1, sep = "")}, strsplit(sj_Mono, "[:-]"))


#### TIGRFAM
cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/TIGRFAM.res.txt", sep = "\t", header = T, stringsAsFactors = F)
length(unique(cdd$Domain))
dim(cdd)
cdd <- cdd[cdd$CDSLength >= cdd$CodingLength,]

allcdd <- read.table("~/Desktop/Domain/TIGRFAM.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)

cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
sum(cdd$ASrange %in% Bcell_sep_range)

cdd <- cdd[cdd$ASrange %in% Bcell_sep_range,]
length(unique(cdd$Domain))

tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
for(i in unique(cdd$Domain)){
  test <- cdd[cdd$Domain == i, ]
  test2 <- allcdd[allcdd$TIGRFAM.domain.ID == i, ]
  tab[i, 1] <- length(unique(test$HostGene))
  tab[i, 2] <- length(unique(cdd$HostGene))
  tab[i, 3] <- length(unique(test2$Gene.stable.ID))
  tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
}

enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)

#pdf("/Users/tangchao/Desktop/Domain/Result/TIGRFAM.pdf")
par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))
#dev.off()

length(enrich_p2);sum(enrich_p2<0.05);sum(enrich_p2<0.01);sum(enrich_p3<0.05);sum(enrich_p3<0.01)
sum(p.adjust(enrich_p2, method = "BH")<0.05);sum(p.adjust(enrich_p2, method = "BH")<0.01);sum(p.adjust(enrich_p3, method = "BH")<0.05);sum(p.adjust(enrich_p3, method = "BH")<0.01)

length(delete_p2);sum(delete_p2<0.05);sum(delete_p2<0.01);sum(delete_p3<0.05);sum(delete_p3<0.01)
sum(p.adjust(delete_p2, method = "BH")<0.05);sum(p.adjust(delete_p2, method = "BH")<0.01);sum(p.adjust(delete_p3, method = "BH")<0.05);sum(p.adjust(delete_p3, method = "BH")<0.01)

for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/TIGRFAM.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  allcdd <- read.table("~/Desktop/Domain/TIGRFAM.txt", sep = "\t", header = T, stringsAsFactors = F)
  allcdd <- na.omit(allcdd)
  
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]

  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$TIGRFAM.domain.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}


#### SMART
allcdd <- read.table("~/Desktop/Domain/SMART.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/SMART.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$SMART.domains.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}


#### PROSITE.profile
allcdd <- read.table("~/Desktop/Domain/PROSITE.profile.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/PROSITE.profile.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$PROSITE.profiles.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}


#### PROSITE.pattern
allcdd <- read.table("~/Desktop/Domain/PROSITE.pattern.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/PROSITE.pattern.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$PROSITE.patterns.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}


#### Prints
allcdd <- read.table("~/Desktop/Domain/Prints.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Prints.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$Prints.domain.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}


#### pfam
allcdd <- read.table("~/Desktop/Domain/pfam.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/pfam.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$Pfam.domain.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}


#### Gene3D
allcdd <- read.table("~/Desktop/Domain/Gene3D.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Gene3D.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$Gene3D.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  delete_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "less")$p.value)
  
  print(length(enrich_p2)); print(sum(enrich_p2<0.05)); print(sum(enrich_p2<0.01)); print(sum(enrich_p3<0.05)); print(sum(enrich_p3<0.01))
  print(sum(p.adjust(enrich_p2, method = "BH")<0.05)); print(sum(p.adjust(enrich_p2, method = "BH")<0.01)); print(sum(p.adjust(enrich_p3, method = "BH")<0.05)); print(sum(p.adjust(enrich_p3, method = "BH")<0.01))
  
  print(length(delete_p2)); print(sum(delete_p2<0.05)); print(sum(delete_p2<0.01)); print(sum(delete_p3<0.05)); print(sum(delete_p3<0.01))
  print(sum(p.adjust(delete_p2, method = "BH")<0.05)); print(sum(p.adjust(delete_p2, method = "BH")<0.01)); print(sum(p.adjust(delete_p3, method = "BH")<0.05)); print(sum(p.adjust(delete_p3, method = "BH")<0.01))
}

par(mfrow = c(2,2))
hist(enrich_p2, main = "prop.test", xlab = "P value")
hist(enrich_p3, main = "fisher.test", xlab = "P value")
hist(p.adjust(enrich_p2, method = "BH"), main = "prop.test",  xlab = "Adjusted P value", xlim = c(0,1))
hist(p.adjust(enrich_p3, method = "BH"), main = "fisher.test",  xlab = "Adjusted P value", xlim = c(0,1))






#### Cell specific heatmap ======================================================================================================================

#### 
#### Gene3D
allcdd <- read.table("~/Desktop/Domain/Gene3D.txt", sep = "\t", header = T, stringsAsFactors = F)
allcdd <- na.omit(allcdd)
for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
  print(j)
  cdd <- read.table("/Users/tangchao/Desktop/Domain/Merge/Gene3D.res.txt", sep = "\t", header = T, stringsAsFactors = F)
  cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
  cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
  
  cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
  tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
  for(i in unique(cdd$Domain)){
    test <- cdd[cdd$Domain == i, ]
    test2 <- allcdd[allcdd$Gene3D.ID == i, ]
    tab[i, 1] <- length(unique(test$HostGene))
    tab[i, 2] <- length(unique(cdd$HostGene))
    tab[i, 3] <- length(unique(test2$Gene.stable.ID))
    tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
  }
  enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
  assign(paste("enrich_p",j, sep = "_"), data.frame(ID = names(enrich_p2), p = enrich_p2, row.names = 1:length(enrich_p2)))
  assign(paste("enrich_f",j, sep = "_"), data.frame(ID = names(enrich_p3), p = enrich_p3, row.names = 1:length(enrich_p3)))
}
my_list_p <- list(enrich_p_Bcell_sep_range, enrich_p_Tcell_sep_range, enrich_p_MK_sep_range, enrich_p_NK_sep_range, enrich_p_Mono_sep_range)
my_list_f <- list(enrich_f_Bcell_sep_range, enrich_f_Tcell_sep_range, enrich_f_MK_sep_range, enrich_f_NK_sep_range, enrich_f_Mono_sep_range)

pro_p <- data.frame(Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_p)[,-1], row.names = Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_p)[,1])
colnames(pro_p) <- c("B cell", "T cell", "MK", "NK", "Monocytes")
fisher_p <- data.frame(Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_f)[,-1], row.names = Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_f)[,1])
colnames(fisher_p) <- c("B cell", "T cell", "MK", "NK", "Monocytes")

pro_p <- pro_p[rowSums(pro_p<0.05, na.rm = T) > 0, ]
pro_p <- -log10(pro_p)
pro_p[is.na(pro_p)] <- 0
library(pheatmap)
library(RColorBrewer)
#breaksList = seq(-1,1,0.01)
#colors = rev(colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(199))
#pheatmap(pro_p, show_rownames = F, breaks = breaksList, color = colors)
pheatmap(pro_p, show_rownames = F)


#pro_p[pro_p>0.05] <- 0.05
#pro_p[pro_p<0] <- -0.05
#breaksList = seq(-0.05,0.05,0.001)
#colors = rev(colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(100))
#pheatmap(pro_p, show_rownames = F, breaks = breaksList, color = colors)




pheatmap(pro_p, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 2.4)(100))

fisher_p <- fisher_p[rowSums(fisher_p<0.05, na.rm = T) > 0, ]
fisher_p <- -log10(fisher_p)
fisher_p[is.na(fisher_p)] <- 0
pheatmap(fisher_p, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = .7)(100))




file.path("/Users/tangchao/Desktop/Domain/Merge/")
doman_info <- list.files(file.path("/Users/tangchao/Desktop/Domain"), pattern = "txt", full.names = T)
doman_res <- list.files(file.path("/Users/tangchao/Desktop/Domain/Merge"), pattern = "res.txt", full.names = T)

all_pro_p.list <- list()
all_fisher_p.list <- list()

for(p in 1:length(doman_info)){
  allcdd <- read.table(doman_info[p], sep = "\t", header = T, stringsAsFactors = F)
  allcdd <- na.omit(allcdd)
  
  for(j in c("Bcell_sep_range", "Tcell_sep_range", "MK_sep_range" ,"NK_sep_range" ,"Mono_sep_range")){
    print(j)
    cdd <- read.table(doman_res[p], sep = "\t", header = T, stringsAsFactors = F)
    cdd <- cdd[cdd$CodingExonLength/cdd$ASlength > 0.8, ]
    cdd$ASrange <- paste(cdd$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(cdd$ASexon,",")),sep = ":")
    
    cdd <- cdd[cdd$ASrange %in% eval(parse(text = j)),]
    if(nrow(cdd)==0){
      #  assign(paste("enrich_p",j, sep = "_"),NULL)
      #  assign(paste("enrich_f",j, sep = "_"),NULL)
      next()
    }
    tab <- data.frame(matrix(NA, nrow = length(unique(cdd$Domain)), ncol = 4), row.names = unique(cdd$Domain))
    for(i in unique(cdd$Domain)){
      test <- cdd[cdd$Domain == i, ]
      test2 <- allcdd[allcdd[,9] == i, ]
      tab[i, 1] <- length(unique(test$HostGene))
      tab[i, 2] <- length(unique(cdd$HostGene))
      tab[i, 3] <- length(unique(test2$Gene.stable.ID))
      tab[i, 4] <- length(unique(allcdd$Gene.stable.ID))
    }
    enrich_p2 <- apply(tab, 1, function(x)prop.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
    enrich_p3 <- apply(tab, 1, function(x)fisher.test(matrix(c(x[1], x[3], x[2] - x[1], x[4] - x[3]), ncol = 2), alternative = "greater")$p.value)
    assign(paste("enrich_p",j, sep = "_"), data.frame(ID = names(enrich_p2), p = enrich_p2, row.names = 1:length(enrich_p2)))
    assign(paste("enrich_f",j, sep = "_"), data.frame(ID = names(enrich_p3), p = enrich_p3, row.names = 1:length(enrich_p3)))
  }
  my_list_p <- list(enrich_p_Bcell_sep_range, enrich_p_Tcell_sep_range, enrich_p_MK_sep_range, enrich_p_NK_sep_range, enrich_p_Mono_sep_range)
  my_list_f <- list(enrich_f_Bcell_sep_range, enrich_f_Tcell_sep_range, enrich_f_MK_sep_range, enrich_f_NK_sep_range, enrich_f_Mono_sep_range)
  #rm(list = c("enrich_p_Bcell_sep_range", "enrich_p_Tcell_sep_range", "enrich_p_MK_sep_range", "enrich_p_NK_sep_range", "enrich_p_Mono_sep_range"))
  
  #my_list_p[!mapply(is.null, my_list_p)]
  
  pro_p <- data.frame(Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_p)[,-1], row.names = Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_p)[,1])
  colnames(pro_p) <- c("B cell", "T cell", "MK", "NK", "Monocytes")
  fisher_p <- data.frame(Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_f)[,-1], row.names = Reduce(function(x,y) {merge(x,y,all=T,by="ID")}, my_list_f)[,1])
  colnames(fisher_p) <- c("B cell", "T cell", "MK", "NK", "Monocytes")
  all_pro_p.list[[p]] <- pro_p
  all_fisher_p.list[[p]] <- fisher_p
}

for(i in 1:length(all_pro_p.list)){
  all_pro_p.list[[i]] <- all_pro_p.list[[i]][row.names(all_pro_p.list[[i]]) %in% unique(read.table(doman_res[i], sep = "\t", header = T, stringsAsFactors = F)$Domain), ]
  all_fisher_p.list[[i]] <- all_fisher_p.list[[i]][row.names(all_fisher_p.list[[i]]) %in% unique(read.table(doman_res[i], sep = "\t", header = T, stringsAsFactors = F)$Domain), ]
}

all_pro_p.list_tu <- all_pro_p.list[mapply(FUN = nrow, all_pro_p.list)>1]
all_fisher_p.list_tu <- all_fisher_p.list[mapply(FUN = nrow, all_fisher_p.list)>1]

mapply(function(x)x[length(x)-1], strsplit(doman_info, "[/\\.]"))[mapply(FUN = nrow, all_pro_p.list)>1]
all_database <- mapply(function(x)x[length(x)-1], strsplit(doman_info, "[/\\.]"))[mapply(FUN = nrow, all_fisher_p.list)>1]

all_pro_p.list_tu2 <- lapply(all_pro_p.list_tu, function(x) x[rowSums(x<0.05, na.rm = T)>0,])
all_fisher_p.list_tu2 <- lapply(all_fisher_p.list_tu, function(x) x[rowSums(x<0.05, na.rm = T)>0,])


all_pro_p <- do.call(rbind, all_pro_p.list_tu2)
all_fisher_p <- do.call(rbind, all_fisher_p.list_tu2)

database_ano <- data.frame(Database = rep(all_database, mapply(nrow, all_pro_p.list_tu2)), row.names = row.names(all_pro_p))

all_pro_p <- -log10(all_pro_p)
all_pro_p[is.na(all_pro_p)] <- 0
library(pheatmap)
library(RColorBrewer)
pheatmap(all_pro_p, show_rownames = F, annotation_row = database_ano,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 2.4)(100))



database_ano <- data.frame(Database = rep(all_database, mapply(nrow, all_fisher_p.list_tu2)), row.names = row.names(all_fisher_p))

all_fisher_p <- -log10(all_fisher_p)
all_fisher_p[is.na(all_fisher_p)] <- 0
library(pheatmap)
library(RColorBrewer)
pheatmap(all_fisher_p, show_rownames = F, annotation_row = database_ano,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 1)(100))


row.names(database_ano)[database_ano$Database != "hmmpanther"]

pheatmap(all_fisher_p[row.names(database_ano)[database_ano$Database != "hmmpanther"],], 
         show_rownames = F, annotation_row = subset.data.frame(database_ano, row.names(database_ano) != "hmmpanther"),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 1)(100))

row.names(database_ano)[database_ano$Database %in% c("hmmpanther", "CDD")]

pheatmap(all_fisher_p[row.names(database_ano)[!database_ano$Database %in% c("hmmpanther", "CDD")],], 
         show_rownames = F, annotation_row = subset.data.frame(database_ano, !row.names(database_ano) %in% c("hmmpanther", "CDD")),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = 1)(100))



do.call(rbind, all_fisher_p.list_tu) -> test
test$ID <- row.names(test)
melt(test, id.vars="ID") -> test2

test2$padjust <- NA
test2[!is.na(test2$value), ]$padjust <- p.adjust(test2[!is.na(test2$value), ]$value, method = "BH")

dcast(test2[,-3], ID ~ variable, value.var="padjust") -> test3

data.frame(test3[,-1], row.names = test3[,1]) -> test4

do.call(rbind, all_fisher_p.list_tu) -> test
dim(test)
dim(test4)
test4$database <- rep(all_database, mapply(nrow, all_fisher_p.list_tu))
test4$database <- factor(test4$database)
test5 <- split(x = test4, f = test4$database)
lapply(test5, function(x) x[,-6]) -> test6
mapply(nrow, test6)

test7 <- lapply(test6, function(x) x[rowSums(x<0.05, na.rm = T)>0,])
mapply(nrow, test7)


test8 <- do.call(rbind, test7)

database_ano <- data.frame(Database = rep(all_database, mapply(nrow, test7)), row.names = row.names(test8))

test8 <- -log10(test8)
test8[is.na(test8)] <- 0
library(pheatmap)
library(RColorBrewer)
pdf("/Users/tangchao/Desktop/Domain/Result/fisher_FDR_Heatmap.pdf")
pheatmap(test8, show_rownames = F, annotation_row = database_ano,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias = .7)(100))
dev.off()




doman_res <- list.files(file.path("/home/zhouran/data/DataBase/Ensembl_bioMart/result"), pattern = "res.txt", full.names = T)

res_list <- lapply(doman_res, function(x) read.table(x, sep = "\t", header = T, stringsAsFactors = F))
res_tab <- do.call(rbind, res_list)
res_tab <- res_tab[res_tab$CodingExonLength/res_tab$ASlength > 0.8, ]

load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
Doman_SE <- unique(paste(res_tab$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(res_tab$ASexon,",")),sep = ":"))

SE_psi[SE_psi$SE_type == "exon-skipping_exactly" & SE_psi$ASrange %in% Doman_SE, grepl("UCB", colnames(SE_psi))] -> Domain_psi
dim(Domain_psi)
#[1] 6959 3574
SE_psi[SE_psi$SE_type == "exon-skipping_exactly" & !SE_psi$ASrange %in% Doman_SE, grepl("UCB", colnames(SE_psi))] -> No_Domain_psi
dim(No_Domain_psi)
#[1] c 3574

library(reshape2)
na.omit(melt(Domain_psi)) -> Domain_psi
na.omit(melt(No_Domain_psi)) -> No_Domain_psi

Domain_psi$variable <- "Domain"
No_Domain_psi$variable <- "No Domain"

test <- rbind(Domain_psi, No_Domain_psi)
colnames(test) <- c("Type", "PSI")
pdf("~/PSI_vioplot.pdf")
library(ggplot2)
ggplot(data = test, aes(x = Type, y = PSI))+
  geom_violin()
dev.off()


res_lis2t <- lapply(res_list, function(y) {unique(paste(y$Chromosome, mapply(function(x) {paste(x, collapse = "-")}, strsplit(y$ASexon,",")),sep = ":"))})

all_database <- mapply(function(x)x[length(x)-2], strsplit(doman_res, "[/\\.]"))

my_list <- list()
for(i in 1:length(res_lis2t)){
  my_list[[i]] <- data.frame(PSI = na.omit(melt(SE_psi[SE_psi$SE_type == "exon-skipping_exactly" & SE_psi$ASrange %in% res_lis2t[[i]], grepl("UCB", colnames(SE_psi))]))$value,
                             DataBase = all_database[i])
}

do.call(rbind, my_list) -> all_tab
colnames(No_Domain_psi) <- c("DataBase", "PSI")

rbind(all_tab, No_Domain_psi[,2:1]) -> all_tabs
all_tabs$DataBase <- factor(all_tabs$DataBase)


library(ggplot2)
library(gridExtra)



pdf("/mnt/data5/BGI/UCB/tangchao/Domain/result/annotated/PSI_vioplot_all_database.pdf")
library(ggplot2)
ggplot(data = all_tabs, aes(x = DataBase, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

library(ggridges)
ggplot(all_tabs, aes(x = PSI, y = DataBase)) +
  geom_density_ridges2() +
  theme_ridges(center_axis_labels = TRUE)
dev.off()



load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_sep_SE.RData")

sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
unlist(lapply(list(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell),length))


Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

dim(SE_psi[SE_psi$loci %in% sj_Bcell, row.names(Cell_type[Cell_type$CellType=="B Cells",])])
[1] 385 607

Bcell_sep <- melt(SE_psi[SE_psi$loci %in% sj_Bcell,row.names(Cell_type)])
Bcell_sep <- merge(x=Bcell_sep, y = Cell_type, by.x = "variable", by.y=0, all.x=T)
colnames(Bcell_sep)[1] <-"UCB"
colnames(Bcell_sep)[2] <-"PSI"

pdf("~/PSI_vioplot_of_Bcell_sep_SE.pdf")
library(ggplot2)
ggplot(data = Bcell_sep, aes(x = CellType, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
dev.off()


CD4_sep <- melt(SE_psi[SE_psi$loci %in% sj_CD4,row.names(Cell_type)])
CD4_sep <- merge(x=CD4_sep, y = Cell_type, by.x = "variable", by.y=0, all.x=T)
colnames(CD4_sep)[1] <-"UCB"
colnames(CD4_sep)[2] <-"PSI"

pdf("~/PSI_vioplot_of_CD4_sep_SE.pdf")
library(ggplot2)
ggplot(data = CD4_sep, aes(x = CellType, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
dev.off()


CD8_sep <- melt(SE_psi[SE_psi$loci %in% sj_CD8,row.names(Cell_type)])
CD8_sep <- merge(x=CD8_sep, y = Cell_type, by.x = "variable", by.y=0, all.x=T)
colnames(CD8_sep)[1] <-"UCB"
colnames(CD8_sep)[2] <-"PSI"

pdf("~/PSI_vioplot_of_CD8_sep_SE.pdf")
library(ggplot2)
ggplot(data = CD8_sep, aes(x = CellType, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
dev.off()


MK_sep <- melt(SE_psi[SE_psi$loci %in% sj_MK,row.names(Cell_type)])
MK_sep <- merge(x=MK_sep, y = Cell_type, by.x = "variable", by.y=0, all.x=T)
colnames(MK_sep)[1] <-"UCB"
colnames(MK_sep)[2] <-"PSI"

pdf("~/PSI_vioplot_of_MK_sep_SE.pdf")
library(ggplot2)
ggplot(data = MK_sep, aes(x = CellType, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
dev.off()


Mono_sep <- melt(SE_psi[SE_psi$loci %in% sj_Mono, row.names(Cell_type)])
Mono_sep <- merge(x=Mono_sep, y = Cell_type, by.x = "variable", by.y=0, all.x=T)
colnames(Mono_sep)[1] <-"UCB"
colnames(Mono_sep)[2] <-"PSI"

pdf("~/PSI_vioplot_of_Mono_sep_SE.pdf")
library(ggplot2)
ggplot(data = Mono_sep, aes(x = CellType, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
dev.off()


NK_sep <- melt(SE_psi[SE_psi$loci %in% sj_NK, row.names(Cell_type)])
NK_sep <- merge(x=NK_sep, y = Cell_type, by.x = "variable", by.y=0, all.x=T)
colnames(NK_sep)[1] <-"UCB"
colnames(NK_sep)[2] <-"PSI"

pdf("~/PSI_vioplot_of_NK_sep_SE.pdf")
library(ggplot2)
ggplot(data = NK_sep, aes(x = CellType, y = PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
dev.off()




SE_psi_cell_sep <- SE_psi[SE_psi$ASrange %in% Doman_SE & SE_psi$loci %in% c(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell), row.names(Cell_type)]
SE_psi_noncell_sep <- SE_psi[SE_psi$ASrange %in% Doman_SE & !SE_psi$loci %in% c(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK, sj_Tcell), row.names(Cell_type)]

cel_sep <- melt(SE_psi_cell_sep)
cel_sep$Category <- "Cell-Type-Specific"

noncel_sep <- melt(SE_psi_noncell_sep)
noncel_sep$Category <- "Common"

cel_sep_Category <- rbind(cel_sep, noncel_sep)
cel_sep_Category <- na.omit(cel_sep_Category)
colnames(cel_sep_Category)[1:2] <- c("Cell", "PSI")

pdf("/mnt/data5/BGI/UCB/tangchao/Domain/result/annotated/PSI_vioplot_of_cell_spe_and_common.pdf")

ggplot(data = cel_sep_Category, aes(x = Category, y = PSI))+
  geom_violin()

library(ggridges)
ggplot(cel_sep_Category, aes(x = PSI, y = Category)) +
  geom_density_ridges2() +
  theme_ridges(center_axis_labels = TRUE)

dev.off()





#### Unannotated SE


doman_res <- list.files(file.path("/mnt/data5/BGI/UCB/tangchao/Domain/result/unannotated"), pattern = "res.txt", full.names = T)

res_list <- lapply(doman_res, function(x) read.delim(x, sep = "\t", header = F, stringsAsFactors = F))
for (i in 1:length(res_list)) {
  colnames(res_list[[i]]) <- as.character(res_list[[i]][1,])
}
for (i in 1:length(res_list)) {
  res_list[[i]] <- res_list[[i]][-1,]
}

#res_tab <- do.call(rbind, res_list)
#res_tab <- res_tab[res_tab$CodingExonLength/res_tab$ASlength > 0.8, ]

load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
table(SE_psi$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                 10372                  3398                  3900

Doman_SE <- mapply(function(x){x[1]}, strsplit(unique(do.call(c, mapply(function(x){x$ASinfo}, res_list))),";"))

SE_psi[SE_psi$SE_type == "novel_exon-skipping" & SE_psi$ASrange %in% Doman_SE, grepl("UCB", colnames(SE_psi))] -> Domain_psi
dim(Domain_psi)
#[1] 473 3574
SE_psi[SE_psi$SE_type == "novel_exon-skipping" & !SE_psi$ASrange %in% Doman_SE, grepl("UCB", colnames(SE_psi))] -> No_Domain_psi
dim(No_Domain_psi)
#[1] 2925 3574

library(reshape2)
na.omit(melt(Domain_psi)) -> Domain_psi
na.omit(melt(No_Domain_psi)) -> No_Domain_psi

Domain_psi$variable <- "Domain"
No_Domain_psi$variable <- "No Domain"

test <- rbind(Domain_psi, No_Domain_psi)
colnames(test) <- c("Type", "PSI")
pdf("~/Unannotated_SE_PSI_vioplot.pdf")
library(ggplot2)
ggplot(data = test, aes(x = Type, y = PSI))+
  geom_violin()

ggplot(data = test[test$PSI>0,], aes(x = Type, y = PSI))+
  geom_violin()

ggplot(data = test[test$PSI>0&test$PSI<1,], aes(x = Type, y = PSI))+
  geom_violin()


par(mfrow = c(3,2))
hist(Domain_psi$value, main = "Domain", xlab = "PSI")
hist(No_Domain_psi$value, main = "No Domain", xlab = "PSI")
hist(Domain_psi$value[Domain_psi$value>0], main = "Domain", xlab = "PSI", breaks = 10)
hist(No_Domain_psi$value[No_Domain_psi$value>0], main = "No Domain", xlab = "PSI", breaks = 10)
hist(Domain_psi$value[Domain_psi$value>0&Domain_psi$value<1], main = "Domain", xlab = "PSI")
hist(No_Domain_psi$value[No_Domain_psi$value>0&No_Domain_psi$value<1], main = "No Domain", xlab = "PSI")


par(mfrow = c(3,2))
barplot(as.numeric(table(cut(Domain_psi$value, seq(0,1,0.1), include.lowest = T))/nrow(Domain_psi)*100), ylim = c(0,100))
barplot(as.numeric(table(cut(No_Domain_psi$value, seq(0,1,0.1), include.lowest = T))/nrow(No_Domain_psi)*100), ylim = c(0,100))
barplot(as.numeric(table(cut(Domain_psi$value, seq(0,1,0.1), include.lowest = F))/nrow(Domain_psi[Domain_psi$value>0,])*100), ylim = c(0,60))
barplot(as.numeric(table(cut(No_Domain_psi$value, seq(0,1,0.1), include.lowest = F))/nrow(No_Domain_psi[No_Domain_psi$value>0,])*100), ylim = c(0,60))
barplot(as.numeric(table(cut(Domain_psi[Domain_psi$value<1,]$value, seq(0,1,0.1), include.lowest = F))/nrow(Domain_psi[Domain_psi$value<1 & Domain_psi$value>0,])*100), ylim = c(0,60))
barplot(as.numeric(table(cut(No_Domain_psi[No_Domain_psi$value<1,]$value, seq(0,1,0.1), include.lowest = F))/nrow(No_Domain_psi[No_Domain_psi$value<1&No_Domain_psi$value>0,])*100), ylim = c(0,60))

dev.off()















