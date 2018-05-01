#### Calculate Shannon entropy and relative entropy of position
## Sat Mar 26 2018
## By Tang Chao

#### 0. Basic settings ======================================================================


depar <- par()
setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/")


## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/entropy/")


## function:
MyEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log2(x))
  Reduce("+", -tmp)
}

round(apply(annotated_donor_percentage, 2, MyEntropy),2)


MyRelativeEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log(x/0.25))
  Reduce("+", tmp)
}


seqLogoFromMatrix <- function(tableResult, plot = TRUE, ic.scale=FALSE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15){
  ## 
  N = ncol(tableResult)
  
  if("seqLogo" %in% installed.packages()){
    library(seqLogo)
  }else{
    source("https://bioconductor.org/biocLite.R")
    biocLite("seqLogo")
    library(seqLogo)
  }
  ## 
  if(is.matrix(tableResult)){
    data_tu <- matrix(rep(NA,4*N),nrow = 4)
    rownames(data_tu) <- c("A", "C", "G", "T")
    
    for(j in 1:N){
      for(i in c("A", "C", "G", "T")){
        data_tu[i,j] <- sum(tableResult[,j] == i)
      }
    }
    data_tu_p <- data_tu/colSums(data_tu)
    if(plot == FALSE){
      return(data_tu_p)
    }else{
      pwm <- makePWM(data_tu_p)
      
      seqLogo(pwm, ic.scale = ic.scale, xaxis = xaxis, yaxis = yaxis, xfontsize = xfontsize, yfontsize = yfontsize)
    }
  }else{
    stop("The input data of tableResult must be a matrix")
  }
}


#### 1. Load data ===========================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/sj_parsed_fastq.RData")


other_donor <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 48, 56), split="")), plot = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 48, 56), split="")), plot = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[annotation == 1]), 145, 153), split="")), plot = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj, ParsedFasta[novel == "Y"]), 145, 153), split="")), plot = F)


round(apply(other_donor, 2, MyEntropy),2)
round(apply(other_donor, 2, MyRelativeEntropy),2)













