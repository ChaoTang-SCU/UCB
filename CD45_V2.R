
setwd("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2")
load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
SJ34 = "1:198692374-198696711"
SJ35 = "1:198692374-198699563"
SJ45 = "1:198696910-198699563"
SJ46 = "1:198696910-198702386"
SJ56 = "1:198699705-198702386"
SJ57 = "1:198699705-198703297"
SJ67 = "1:198702531-198703297"
SJ47 = "1:198696910-198703297"
SJ36 = "1:198692374-198702386"
SJ37 = "1:198692374-198703297"

as.data.frame(te[.(c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37))])[,-1] -> CD45_SJ
CD45_SJ[is.na(CD45_SJ)] <- 0
CD45_SJ[CD45_SJ<10] <- 0
colnames(CD45_SJ) <- substr(colnames(CD45_SJ),1,10)
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

sum(colSums(CD45_SJ[1:9,])>=10)
# [1] 2348
CD45_SJ[1:9,colSums(CD45_SJ[1:9,])<10] <- 0
## sum of junction reads (except the cells only have SJ37) must more than 10.
test <- as.data.frame(apply(CD45_SJ[1:9,],2, function(x) x/sum(x) <= 0.05))
test[is.na(test)] <- TRUE
test["SJ37",] <- FALSE
CD45_SJ[as.matrix(test)] <- 0



CD45_R_r_tab <- CD45_SJ[1:8,]
row.names(CD45_R_r_tab) <- c("RO","RAC","RBC","RAB","RA","RC","RB","RABC")
CD45_R_r_tab[CD45_R_r_tab>0] <- 0

for (i in 1:ncol(CD45_SJ)) {
  if(CD45_SJ["SJ37",i] > 0){
    CD45_R_r_tab["RO",i] <- CD45_SJ["SJ37",i]
  }
  if(CD45_SJ["SJ46",i] > 0){
    CD45_R_r_tab["RAC",i] <- CD45_SJ["SJ46",i]
  }
  if(CD45_SJ["SJ35",i] > 0 & CD45_SJ["SJ56",i] > 0 & CD45_SJ["SJ67",i] > 0){
  	if(CD45_SJ["SJ34",i] > 0 & CD45_SJ["SJ45",i] > 0){
  		CD45_R_r_tab["RBC",i] <- CD45_SJ["SJ35",i]
  	}else{
  		CD45_R_r_tab["RBC",i] <- mean(c(CD45_SJ["SJ56",i], CD45_SJ["SJ67",i]))
  	}#yin wei dang RBC yu RB tong shi cun zai de shi hou, yong J35 ji suan PSI RBC hui bu dui.
  }
  if(CD45_SJ["SJ45",i] > 0 & CD45_SJ["SJ57",i] > 0 & CD45_SJ["SJ34",i] > 0){
  	if(CD45_SJ["SJ56",i] > 0 & CD45_SJ["SJ67",i] > 0){
  		CD45_R_r_tab["RAB",i] <- CD45_SJ["SJ57",i]
  	}else{
  		CD45_R_r_tab["RAB",i] <- mean(c(CD45_SJ["SJ34",i], CD45_SJ["SJ45",i]))
  	}    
  }
  if(CD45_SJ["SJ47",i] > 0 & CD45_SJ["SJ34",i] > 0){
    CD45_R_r_tab["RA",i] <- CD45_SJ["SJ47",i]
  }
  if(CD45_SJ["SJ36",i] > 0 & CD45_SJ["SJ67",i] > 0){
    CD45_R_r_tab["RC",i] <- CD45_SJ["SJ36",i]
  }
  if(CD45_SJ["SJ35",i] > 0 & CD45_SJ["SJ57",i] > 0){
  	if(sum(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),i])>0){
  		CD45_R_r_tab["RB",i] <- min(c(CD45_SJ["SJ35",i], CD45_SJ["SJ57",i]))
  		}else{
  			CD45_R_r_tab["RB",i] <- mean(c(CD45_SJ["SJ35",i], CD45_SJ["SJ57",i]))
  		}
  }# dang RB yu RAB huozhe RBC tongshi cunzai de shihou, xuyao jianqu RAB huozhe RBC de reads.
  # zheli women congshi zhizhong dou meiyou kaolv tongshi cunzai 4 zhong isoform de qingkuang, shi yinwei women houmian yao jiang 4 zhong isoform de qingkuang tiaochulai chongxin dingyi.
  if(CD45_SJ["SJ34",i] > 0 & CD45_SJ["SJ45",i] > 0 & CD45_SJ["SJ56",i] > 0 & CD45_SJ["SJ67",i] > 0){
  	if(CD45_SJ["SJ35",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ34",i], CD45_SJ["SJ45",i]))
}
if(CD45_SJ["SJ36",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ34",i], CD45_SJ["SJ45",i],CD45_SJ["SJ56",i]))
}
if(CD45_SJ["SJ46",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ45",i], CD45_SJ["SJ56",i]))
}
if(CD45_SJ["SJ47",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ45",i], CD45_SJ["SJ56",i], CD45_SJ["SJ67",i]))
}
if(CD45_SJ["SJ57",i] > 0){
	CD45_R_r_tab["RABC",i] <- mean(c(CD45_SJ["SJ56",i], CD45_SJ["SJ67",i]))
}
if(CD45_SJ["SJ35",i] == 0 & CD45_SJ["SJ36",i] == 0 & CD45_SJ["SJ46",i] == 0 & CD45_SJ["SJ47",i] == 0 & CD45_SJ["SJ57",i] == 0){
    CD45_R_r_tab["RABC",i] <- (CD45_SJ["SJ34",i] + CD45_SJ["SJ45",i] + CD45_SJ["SJ56",i] + CD45_SJ["SJ67",i])/4
}
}
}





CD45_R_r_tab[,1:5]
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==0])
# [1]    8 1486
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==1])
# [1]    8 1480
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2])
# [1]    8 452
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3])
# [1]    8 22
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==4])
# [1]    8 129
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==5])
# [1]    8 5

dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==0])
# [1]   8 1486
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==1])
# [1]   8 1480
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==2])
# [1]   8 452
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==3])
# [1]   8 22
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==4])
# [1]   8 129
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==5])
# [1]   8 5




CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2][,1:5]
CD45_SJ[,colSums(CD45_R_r_tab!=0)==2][,1:5]

CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3][,1:5]
CD45_SJ[,colSums(CD45_R_r_tab!=0)==3][,1:5]

CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==4][,1:5]
CD45_SJ[,colSums(CD45_R_r_tab!=0)==4][,1:5]

recalculete <- colnames(CD45_R_r_tab)[c(colSums(CD45_R_r_tab!=0)==4 | colSums(CD45_R_r_tab!=0)==5)]
dim(CD45_R_r_tab[, recalculete])
# [1]   8 134
dim(CD45_SJ[, recalculete])
# [1]  10 134
rowSums(CD45_SJ[, recalculete])
#  SJ34  SJ35  SJ45  SJ46  SJ56  SJ57  SJ67  SJ47  SJ36  SJ37
# 80656 90594 72565     0 73570 78443 75796     0     0  1524

for(i in 1:length(recalculete)){
	# -+-, -+-
	if(CD45_SJ["SJ35",recalculete[i]] > CD45_SJ["SJ34",recalculete[i]] & CD45_SJ["SJ35",recalculete[i]] > CD45_SJ["SJ45",recalculete[i]]&
		CD45_SJ["SJ57",recalculete[i]] > CD45_SJ["SJ56",recalculete[i]] & CD45_SJ["SJ57",recalculete[i]] > CD45_SJ["SJ67",recalculete[i]]){
		CD45_R_r_tab["RB",recalculete[i]] <- mean(c(CD45_SJ["SJ35",recalculete[i]], CD45_SJ["SJ57",recalculete[i]]))
		CD45_R_r_tab["RAB",recalculete[i]] <- 0
		CD45_R_r_tab["RBC",recalculete[i]] <- 0
		CD45_R_r_tab["RABC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/4
	}
	# +-+, -+-
	if(CD45_SJ["SJ35",recalculete[i]] < CD45_SJ["SJ34",recalculete[i]] & CD45_SJ["SJ35",recalculete[i]] < CD45_SJ["SJ45",recalculete[i]]&
		CD45_SJ["SJ57",recalculete[i]] > CD45_SJ["SJ56",recalculete[i]] & CD45_SJ["SJ57",recalculete[i]] > CD45_SJ["SJ67",recalculete[i]]){
		CD45_R_r_tab["RB",recalculete[i]] <- 0
		CD45_R_r_tab["RAB",recalculete[i]] <- (CD45_SJ["SJ35",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/3
		CD45_R_r_tab["RBC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ57",recalculete[i]])/3
		CD45_R_r_tab["RABC",recalculete[i]] <- 0
	}
	# -+-, +-+
	if(CD45_SJ["SJ35",recalculete[i]] > CD45_SJ["SJ34",recalculete[i]] & CD45_SJ["SJ35",recalculete[i]] > CD45_SJ["SJ45",recalculete[i]]&
		CD45_SJ["SJ57",recalculete[i]] < CD45_SJ["SJ56",recalculete[i]] & CD45_SJ["SJ57",recalculete[i]] < CD45_SJ["SJ67",recalculete[i]]){
		CD45_R_r_tab["RB",recalculete[i]] <- 0
		CD45_R_r_tab["RAB",recalculete[i]] <- (CD45_SJ["SJ35",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/3
		CD45_R_r_tab["RBC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ57",recalculete[i]])/3
		CD45_R_r_tab["RABC",recalculete[i]] <- 0
	}
	# +-+, +-+
	if(CD45_SJ["SJ35",recalculete[i]] < CD45_SJ["SJ34",recalculete[i]] & CD45_SJ["SJ35",recalculete[i]] < CD45_SJ["SJ45",recalculete[i]]&
		CD45_SJ["SJ57",recalculete[i]] < CD45_SJ["SJ56",recalculete[i]] & CD45_SJ["SJ57",recalculete[i]] < CD45_SJ["SJ67",recalculete[i]]){
		CD45_R_r_tab["RB",recalculete[i]] <- mean(c(CD45_SJ["SJ35",recalculete[i]], CD45_SJ["SJ57",recalculete[i]]))
		CD45_R_r_tab["RAB",recalculete[i]] <- 0
		CD45_R_r_tab["RBC",recalculete[i]] <- 0
		CD45_R_r_tab["RABC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/4
	}
}

dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==4])
[1]  8 10
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==5])
[1] 8 0

CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==4]
     UCB3.00062 UCB3.00158 UCB3.00217 UCB3.00489 UCB3.00860 UCB4.00503
RO          0.0          0          0        0.0        0.0          0
RAC         0.0          0          0        0.0        0.0          0
RBC        82.0        232         21      154.0      180.0        130
RAB        94.0        222         17      106.0      420.0        110
RA          0.0          0          0        0.0        0.0          0
RC          0.0          0          0        0.0        0.0          0
RB        176.0        454         38      260.0      600.0        240
RABC       97.5        227         16      158.5      401.5         99
     UCB4.00730 UCB4.00797 UCB4.00965 UCB4.00983
RO            0        0.0        0.0        0.0
RAC           0        0.0        0.0        0.0
RBC         154      468.0       53.0      520.0
RAB          75      489.0       52.0      184.0
RA            0        0.0        0.0        0.0
RC            0        0.0        0.0        0.0
RB          229      957.0      105.0      704.0
RABC         72      193.5       20.5      415.5
CD45_SJ[,colSums(CD45_R_r_tab!=0)==4]
     UCB3.00062 UCB3.00158 UCB3.00217 UCB3.00489 UCB3.00860 UCB4.00503
SJ34         99        652         27        200        772        135
SJ35         82        232         21        154        180        130
SJ45         88        356         17        176        746         99
SJ46          0          0          0          0          0          0
SJ56         86        242         20         92        458         97
SJ57         94        222         17        106        420        110
SJ67        109        212         12        225        345        101
SJ47          0          0          0          0          0          0
SJ36          0          0          0          0          0          0
SJ37          0          0          0          0          0          0
     UCB4.00730 UCB4.00797 UCB4.00965 UCB4.00983
SJ34        117        498         59        598
SJ35        154        468         53        520
SJ45         56        443         36        350
SJ46          0          0          0          0
SJ56         77        209         21        484
SJ57         75        489         52        184
SJ67         67        178         20        347
SJ47          0          0          0          0
SJ36          0          0          0          0
SJ37          0          0          0          0


for(i in 1:length(recalculete)){
	# -+-, -+-
	if(CD45_SJ["SJ35",recalculete[i]] > mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] > mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- mean(c(CD45_SJ["SJ35",recalculete[i]], CD45_SJ["SJ57",recalculete[i]]))
		CD45_R_r_tab["RAB",recalculete[i]] <- 0
		CD45_R_r_tab["RBC",recalculete[i]] <- 0
		CD45_R_r_tab["RABC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/4
	}
	# +-+, -+-
	if(CD45_SJ["SJ35",recalculete[i]] < mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] > mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- 0
		CD45_R_r_tab["RAB",recalculete[i]] <- (CD45_SJ["SJ35",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/3
		CD45_R_r_tab["RBC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ57",recalculete[i]])/3
		CD45_R_r_tab["RABC",recalculete[i]] <- 0
	}
	# -+-, +-+
	if(CD45_SJ["SJ35",recalculete[i]] > mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] < mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- 0
		CD45_R_r_tab["RAB",recalculete[i]] <- (CD45_SJ["SJ35",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/3
		CD45_R_r_tab["RBC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ57",recalculete[i]])/3
		CD45_R_r_tab["RABC",recalculete[i]] <- 0
	}
	# +-+, +-+
	if(CD45_SJ["SJ35",recalculete[i]] < mean(c(CD45_SJ["SJ34",recalculete[i]], CD45_SJ["SJ45",recalculete[i]])) &
		CD45_SJ["SJ57",recalculete[i]] < mean(c(CD45_SJ["SJ56",recalculete[i]], CD45_SJ["SJ67",recalculete[i]]))){
		CD45_R_r_tab["RB",recalculete[i]] <- mean(c(CD45_SJ["SJ35",recalculete[i]], CD45_SJ["SJ57",recalculete[i]]))
		CD45_R_r_tab["RAB",recalculete[i]] <- 0
		CD45_R_r_tab["RBC",recalculete[i]] <- 0
		CD45_R_r_tab["RABC",recalculete[i]] <- (CD45_SJ["SJ34",recalculete[i]] + CD45_SJ["SJ45",recalculete[i]] + CD45_SJ["SJ56",recalculete[i]] + CD45_SJ["SJ67",recalculete[i]])/4
	}
}


dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==4])
# [1]   8 0
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab!=0)==5])
# [1]   8 0


rowSums(CD45_R_r_tab>0)
  RO  RAC  RBC  RAB   RA   RC   RB RABC
 160    1  788  263    2    3  758  748




dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==0])
# [1]    8 1486
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==1])
# [1]    8 1480
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2])
# [1]    8 581
dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3])
# [1]    8 27


write.table(CD45_R_r_tab,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2/CD45_RABC.txt", sep = "\t")
write.table(CD45_SJ,"/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2/CD45_SJ.txt", sep = "\t")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

library(ggplot2)
t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,1]),]) -> CD45_R_r_tab1
f1 <- ggplot(data = tsne[colnames(CD45_R_r_tab1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab1[1,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab2)[1])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,2]),]) -> CD45_R_r_tab2
f2 <- ggplot(data = tsne[colnames(CD45_R_r_tab2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab2[2,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab2)[2])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,3]),]) -> CD45_R_r_tab3
f3 <- ggplot(data = tsne[colnames(CD45_R_r_tab3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab3[3,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab3)[3])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,4]),]) -> CD45_R_r_tab4
f4 <- ggplot(data = tsne[colnames(CD45_R_r_tab4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab4[4,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab4)[4])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,5]),]) -> CD45_R_r_tab5
f5 <- ggplot(data = tsne[colnames(CD45_R_r_tab5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab5[5,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab5)[5])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,6]),]) -> CD45_R_r_tab6
f6 <- ggplot(data = tsne[colnames(CD45_R_r_tab6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab6[6,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab6)[6])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,7]),]) -> CD45_R_r_tab7
f7 <- ggplot(data = tsne[colnames(CD45_R_r_tab7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab7[7,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab7)[7])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(CD45_R_r_tab)[order(t(CD45_R_r_tab)[,8]),]) -> CD45_R_r_tab8
f8 <- ggplot(data = tsne[colnames(CD45_R_r_tab8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(CD45_R_r_tab8[8,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(CD45_R_r_tab8)[8])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2/CD45_RABC_v2.pdf", height = 5,width = 10)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, ncol = 4)
dev.off()



dim(CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==1])

CD45_R_r_tab -> chao
chao[,colSums(CD45_R_r_tab>0)>1] <- 0
dim(chao[,colSums(chao>0)==1])
dim(chao[,colSums(chao>0)>1])
#[1] 8 0
chao <- chao[rowSums(chao)>0,]
dim(chao)


library(ggplot2)
t(t(chao)[order(t(chao)[,1]),]) -> chao1
f1 <- ggplot(data = tsne[colnames(chao1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao1[1,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao2)[1])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,2]),]) -> chao2
f2 <- ggplot(data = tsne[colnames(chao2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao2[2,]+1))), size = .6)+
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
  geom_point(aes(colour = log10(as.numeric(chao3[3,]+1))), size = .6)+
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
  geom_point(aes(colour = log10(as.numeric(chao4[4,]+1))), size = .6)+
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
  geom_point(aes(colour = log10(as.numeric(chao5[5,]+1))), size = .6)+
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
  geom_point(aes(colour = log10(as.numeric(chao6[6,]+1))), size = .6)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao6)[6])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v2/CD45_RABC_v2_only_one_type.pdf", height = 4,width = 6)
library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, ncol = 3)
dev.off()


rowSums(chao>0)





