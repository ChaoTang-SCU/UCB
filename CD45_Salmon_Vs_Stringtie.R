st1 <- read.csv("/Users/tangchao/CloudStation/tangchao/project/UCB/Experiment_Validate/CD45_1/one_isoform_todo.csv", row.names = 1)
Cells <- c("UCB1.00091","UCB4.01155","UCB4.00206","UCB3.01077","UCB1.00252","UCB3.01523","UCB3.00718","UCB3.00730","UCB4.01381","UCB3.02063","UCB1.00070","UCB3.00594","UCB3.00309","UCB3.00538","UCB3.01856","UCB1.00039","UCB1.00046","UCB4.01527")

st1 <- st1[, colnames(st1) %in% Cells]
apply(st1,2,function(x)x/sum(x))
library(reshape2)
melt(apply(st1,2,function(x)x/sum(x)))

st2 <- read.csv("/Users/tangchao/CloudStation/tangchao/project/UCB/Experiment_Validate/CD45_2/two_isoform_todo.csv", row.names = 1)
st2 <- st2[, colnames(st2) %in% Cells]
melt(apply(st2,2,function(x) x/sum(x)))

st3 <- read.csv("/Users/tangchao/CloudStation/tangchao/project/UCB/Experiment_Validate/CD45_3/three_isoform_todo.csv", row.names = 1)
st3 <- st3[, colnames(st3) %in% Cells]
st4 <- read.csv("/Users/tangchao/CloudStation/tangchao/project/UCB/Experiment_Validate/CD45_4/four_isoform_todo.csv", row.names = 1)
st4 <- st4[, colnames(st4) %in% Cells]
st5 <- read.csv("/Users/tangchao/CloudStation/tangchao/project/UCB/Experiment_Validate/CD45_5/five_isoform_todo.csv", row.names = 1)
st5 <- st5[, colnames(st5) %in% Cells]

st <- cbind(st1,st2,st3,st4,st5)
melt(apply(st,2,function(x) x/sum(x)))

load(file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")
sm <- cov_tu[,c("UCB1.00091","UCB4.01155","UCB4.00206","UCB3.01077","UCB1.00252","UCB3.01523","UCB3.00718","UCB3.00730","UCB4.01381","UCB3.02063","UCB1.00070","UCB3.00594","UCB3.00309","UCB3.00538","UCB3.01856","UCB1.00039","UCB1.00046","UCB4.01527")]
melt(apply(sm,2,function(x) x/sum(x)))

load("~/Desktop/sm.RData")

st <- st[, colnames(sm)]

test <- cbind(melt(apply(st,2,function(x) x/sum(x))), melt(apply(sm,2,function(x) x/sum(x))))

test <- test[, c(2,1,3,6)]
colnames(test) <- c("Cell", "Isoform", "Stringtie", "Salmon")
test$Stringtie <- round(test$Stringtie*100,2)
test$Salmon <- round(test$Salmon*100,2)
write.csv(test, "/Users/tangchao/CloudStation/tangchao/project/UCB/CD45_sort_out/Isoform_ratio.csv", row.names = F, quote = F)


test <- read.csv("/Users/tangchao/CloudStation/tangchao/project/UCB/CD45_sort_out/Isoform_ratio.csv")
cor(test$Stringtie, test$Salmon)
cor.test(test$Stringtie, test$Salmon)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/CD45_sort_out/correlation_between_stringtie_and_salmon.pdf")
plot(test$Stringtie, test$Salmon, xlab = "Stringtie", ylab = "Salmon", pch = 19)
abline(lm(test$Salmon~test$Stringtie), col = 2)
lab = paste("cor = ", round(cor(test$Stringtie, test$Salmon),2), "\n", "p = ", signif(cor.test(test$Stringtie, test$Salmon)$p.value,2), "\n", "R2 = ", "0.94", sep = "")
legend("topleft", legend = lab, bty = "n")
dev.off()



library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
library(ggplot2)
ggplot(data = test, aes(x = Stringtie, y = Salmon, colour = Isoform))+
  geom_point()+
  geom_abline(intercept = coef(lm(test$Salmon~test$Stringtie))[1], slope = coef(lm(test$Salmon~test$Stringtie))[2])+
  annotation_custom(grob = grob)

test$Isoform <- factor(test$Isoform)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/CD45_sort_out/correlation_between_stringtie_and_salmon.pdf")
ggplot(data = test, aes(x = Stringtie, y = Salmon, colour = Isoform, shape = Isoform))+
  geom_point()+
  geom_abline(intercept = coef(lm(test$Salmon~test$Stringtie))[1], slope = coef(lm(test$Salmon~test$Stringtie))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(test$Isoform))
dev.off()



library(gdata)
test <- read.xls("/Users/tangchao/CloudStation/LabShare/Projects/UCB/Sup_materials/CD45 transcripts order.xlsx", sheet = 2)
test <- na.omit(test)


pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/CD45_sort_out/correlation_between_RT_PCR_and_salmon.pdf")
summary(lm(test$Salmon~test$RT.PCR))
lab = paste("cor = ", round(cor(test$RT.PCR, test$Salmon),2), "\n", "p = ", signif(cor.test(test$RT.PCR, test$Salmon)$p.value,2), "\n", "R2 = ", "0.73", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
ggplot(data = test, aes(x = RT.PCR, y = Salmon, colour = Isoform, shape = Isoform))+
  geom_point()+
  geom_abline(intercept = coef(lm(test$Salmon~test$RT.PCR))[1], slope = coef(lm(test$Salmon~test$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(test$Isoform))
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/CD45_sort_out/correlation_between_RT_PCR_and_Stringtie.pdf")
summary(lm(test$Stringtie~test$RT.PCR))
lab = paste("cor = ", round(cor(test$RT.PCR, test$Stringtie),2), "\n", "p = ", signif(cor.test(test$RT.PCR, test$Stringtie)$p.value,2), "\n", "R2 = ", "0.67", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
ggplot(data = test, aes(x = RT.PCR, y = Stringtie, colour = Isoform, shape = Isoform))+
  geom_point()+
  geom_abline(intercept = coef(lm(test$Stringtie~test$RT.PCR))[1], slope = coef(lm(test$Stringtie~test$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(test$Isoform))
dev.off()









