load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
SJ1 = "1:198692374-198696711"
SJ2 = "1:198692374-198699563"
SJ3 = "1:198696910-198699563"
SJ4 = "1:198696910-198702386"
SJ5 = "1:198699705-198702386"
SJ6 = "1:198699705-198703297"
SJ7 = "1:198702531-198703297"
SJ8 = "1:198692374-198703297"
SJ9 = "1:198696910-198703297"
SJ10= "1:198692374-198702386"

te[.(c(SJ1,SJ2,SJ3,SJ4,SJ5,SJ6,SJ7,SJ8,SJ9,SJ10)),1:5]
as.data.frame(te[.(c(SJ1,SJ2,SJ3,SJ4,SJ5,SJ6,SJ7,SJ8,SJ9,SJ10))])[,-1] -> CD45_SJ

CD45_SJ[is.na(CD45_SJ)] <- 0
CD45_SJ[CD45_SJ<2] <- 0


load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)


#### Use the juction reads information also =======================================================================================


CD45_R_r_tab <- CD45_SJ[1:9,]
row.names(CD45_R_r_tab) <- c("CD45-","RO","RAC","RBC","RAB","RA","RC","RB","RABC")
CD45_R_r_tab[CD45_R_r_tab>0] <- 0
for (i in 1:ncol(CD45_SJ)) {
  if(sum(CD45_SJ[,i]) == 0){
    CD45_R_r_tab["CD45-",i] <- 1
  }
  if(CD45_SJ[8,i] > 0){
    CD45_R_r_tab["RO",i] <- CD45_SJ[8,i]
  }
  if(CD45_SJ[4,i] > 0){
    CD45_R_r_tab["RAC",i] <- CD45_SJ[4,i]
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0){
    CD45_R_r_tab["RBC",i] <- CD45_SJ[2,i]
  }
  if(CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] > 0){
    CD45_R_r_tab["RAB",i] <- CD45_SJ[6,i]
  }
  if(CD45_SJ[9,i] > 0){
    CD45_R_r_tab["RA",i] <- CD45_SJ[9,i]
  }
  if(CD45_SJ[10,i] > 0){
    CD45_R_r_tab["RC",i] <- CD45_SJ[10,i]
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[6,i] > 0){
    CD45_R_r_tab["RB",i] <- mean(CD45_SJ[2,i], CD45_SJ[6,i])
  }
  if(sum(CD45_SJ[,i]) > 0 & CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0 & CD45_SJ[8,i] == 0 & CD45_SJ[9,i] == 0 & CD45_SJ[10,i] == 0){
    CD45_R_r_tab["RABC",i] <- mean(CD45_SJ[c(1,3,5,7),i])
}
}


colnames(CD45_R_r_tab) <- substr(colnames(CD45_R_r_tab),1,10)
CD45_R_r_tab[,row.names(tsne)] -> CD45_R_r_tab

CD45_R_r_tab[,1:5]
CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3]
CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2][,1:10]

apply(CD45_R_r_tab, 2, function(x) round(x/sum(x),2)) -> tang

tang[,colSums(tang>0)==2] -> test

hist(as.numeric(apply(test,2,function(x){ 1-max(x/sum(x)) })), breaks = 20, main = "", xlab = "Fraction of minor subset")
plot(density(as.numeric(apply(test,2,function(x){ 1-max(x/sum(x)) })), adjust = 0.2 ), main = "", xlab = "Fraction of minor subset")
abline(v = 0.05, col = 2)


tang[tang < 0.05] <- 0
tang[,colSums(tang>0)==2][,1:10]
tang[tang>0] <- 1


dim(CD45_R_r_tab)
dim(tang)
identical(colnames(CD45_R_r_tab),colnames(tang))

CD45_R_r_tab -> chao
chao[tang==0] <- 0
apply(chao,1,summary)


write.table(chao,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_RABC.txt", sep = "\t")
write.table(CD45_SJ,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_SJ.txt", sep = "\t")




t(t(chao)[order(t(chao)[,1]),]) -> chao1
f1 <- ggplot(data = tsne[colnames(chao1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao1[1,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle("R-")+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,2]),]) -> chao2
f2 <- ggplot(data = tsne[colnames(chao2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao2[2,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao2)[2])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,3]),]) -> chao3
f3 <- ggplot(data = tsne[colnames(chao3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao3[3,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao3)[3])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,4]),]) -> chao4
f4 <- ggplot(data = tsne[colnames(chao4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao4[4,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao4)[4])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,5]),]) -> chao5
f5 <- ggplot(data = tsne[colnames(chao5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao5[5,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao5)[5])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,6]),]) -> chao6
f6 <- ggplot(data = tsne[colnames(chao6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao6[6,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao6)[6])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,7]),]) -> chao7
f7 <- ggplot(data = tsne[colnames(chao7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao7[7,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao7)[7])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,8]),]) -> chao8
f8 <- ggplot(data = tsne[colnames(chao8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao8[8,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao8)[8])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,9]),]) -> chao9
f9 <- ggplot(data = tsne[colnames(chao9),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao9[9,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao9)[9])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())



pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_RABC.pdf")
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, f9, ncol = 3)
dev.off()






#### 10 reads cutoff ==================================================================================================================


CD45_SJ[is.na(CD45_SJ)] <- 0
CD45_SJ[CD45_SJ<10] <- 0


CD45_R_r_tab <- CD45_SJ[1:9,]
row.names(CD45_R_r_tab) <- c("CD45-","RO","RAC","RBC","RAB","RA","RC","RB","RABC")
CD45_R_r_tab[CD45_R_r_tab>0] <- 0
for (i in 1:ncol(CD45_SJ)) {
  if(sum(CD45_SJ[,i]) == 0){
    CD45_R_r_tab["CD45-",i] <- 1
  }
  if(CD45_SJ[8,i] > 0){
    CD45_R_r_tab["RO",i] <- CD45_SJ[8,i]
  }
  if(CD45_SJ[4,i] > 0){
    CD45_R_r_tab["RAC",i] <- CD45_SJ[4,i]
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0){
    CD45_R_r_tab["RBC",i] <- CD45_SJ[2,i]
  }
  if(CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] > 0){
    CD45_R_r_tab["RAB",i] <- CD45_SJ[6,i]
  }
  if(CD45_SJ[9,i] > 0){
    CD45_R_r_tab["RA",i] <- CD45_SJ[9,i]
  }
  if(CD45_SJ[10,i] > 0){
    CD45_R_r_tab["RC",i] <- CD45_SJ[10,i]
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[6,i] > 0){
    CD45_R_r_tab["RB",i] <- mean(CD45_SJ[2,i], CD45_SJ[6,i])
  }
  if(sum(CD45_SJ[,i]) > 0 & CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0 & CD45_SJ[8,i] == 0 & CD45_SJ[9,i] == 0 & CD45_SJ[10,i] == 0){
    CD45_R_r_tab["RABC",i] <- mean(CD45_SJ[c(1,3,5,7),i])
}
}


colnames(CD45_R_r_tab) <- substr(colnames(CD45_R_r_tab),1,10)
CD45_R_r_tab[,row.names(tsne)] -> CD45_R_r_tab

CD45_R_r_tab[,1:5]
CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3]
CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2][,1:10]

apply(CD45_R_r_tab, 2, function(x) round(x/sum(x),2)) -> tang

tang[,colSums(tang>0)==2] -> test


hist(as.numeric(apply(test,2,function(x){ 1-max(x/sum(x)) })), breaks = 20, main = "", xlab = "Fraction of minor subset")
plot(density(as.numeric(apply(test,2,function(x){ 1-max(x/sum(x)) })), adjust = 0.2 ), main = "", xlab = "Fraction of minor subset")
abline(v = 0.05, col = 2)


tang[tang < 0.05] <- 0
tang[,colSums(tang>0)==2][,1:10]
tang[tang>0] <- 1


dim(CD45_R_r_tab)
dim(tang)
identical(colnames(CD45_R_r_tab),colnames(tang))

CD45_R_r_tab -> chao
chao[tang==0] <- 0
apply(chao,1,summary)


write.table(chao,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_RABC_10reads.txt", sep = "\t")
write.table(CD45_SJ,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_SJ_10reads.txt", sep = "\t")


t(t(chao)[order(t(chao)[,1]),]) -> chao1
f1 <- ggplot(data = tsne[colnames(chao1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao1[1,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle("R-")+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,2]),]) -> chao2
f2 <- ggplot(data = tsne[colnames(chao2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao2[2,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao2)[2])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,3]),]) -> chao3
f3 <- ggplot(data = tsne[colnames(chao3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao3[3,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao3)[3])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,4]),]) -> chao4
f4 <- ggplot(data = tsne[colnames(chao4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao4[4,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao4)[4])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,5]),]) -> chao5
f5 <- ggplot(data = tsne[colnames(chao5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao5[5,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao5)[5])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,6]),]) -> chao6
f6 <- ggplot(data = tsne[colnames(chao6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao6[6,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao6)[6])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,7]),]) -> chao7
f7 <- ggplot(data = tsne[colnames(chao7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao7[7,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao7)[7])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,8]),]) -> chao8
f8 <- ggplot(data = tsne[colnames(chao8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao8[8,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao8)[8])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,9]),]) -> chao9
f9 <- ggplot(data = tsne[colnames(chao9),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao9[9,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao9)[9])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())



pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_RABC_10reads.pdf")
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, f9, ncol = 3)
dev.off()






