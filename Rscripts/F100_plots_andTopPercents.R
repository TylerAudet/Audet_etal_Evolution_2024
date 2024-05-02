library(dplyr)
require(data.table)
library(ggplot2)
library(cowplot)

# not provided but generated from grenedalf via the F100 data from Turner et al. 2011
F100_fst <- read.csv("../data/F100_filteredSync_10000bp_windowsfst.csv", header = T)
F100_pi <- read.csv("../data/F100_50000_pidiversity.csv", header = T)

L1_S1 <- F100_fst[,c(1,2,3,15)]
L1_S2 <- F100_fst[,c(1,2,3,16)]
L2_S1 <- F100_fst[,c(1,2,3,17)]
L2_S2 <- F100_fst[,c(1,2,3,18)]

fst_95th <- function(sample) {
  colnames(sample) <- c("chrom","chromStart","chromEnd","fst")
  index <- sample[,4] > quantile(sample[,4], 0.95, na.rm = T)
  print(sample[index,])
}

L1_S1_95th <- na.omit(fst_95th(L1_S1))
L1_S2_95th <- na.omit(fst_95th(L1_S2))
L2_S1_95th <- na.omit(fst_95th(L2_S1))
L2_S2_95th <- na.omit(fst_95th(L2_S2))

write.table(L1_S1_95th, "../data/F100_L1S1_top95_fst.txt")
write.table(L1_S2_95th, "../data/F100_L1S2_top95_fst.txt")
write.table(L2_S1_95th, "../data/F100_L2S1_top95_fst.txt")
write.table(L2_S2_95th, "../data/F100_L2S2_top95_fst.txt")

belowAverage <- function(sample) {
  colnames(sample) <- c("chrom","chromStart","chromEnd","SNPs","coverage","absPi","pi")
  index <- sample[,7] < quantile(sample[,7], 0.1)
  print(sample[index,])
}


C1.pi <- F100_pi[,c(1,2,3,4,5,6,7)]
C2.pi <- F100_pi[,c(1,2,3,8,9,10,11)]
L1.pi <- F100_pi[,c(1,2,3,12,13,14,15)]
L2.pi <- F100_pi[,c(1,2,3,16,17,18,19)]
S1.pi <- F100_pi[,c(1,2,3,20,21,22,23)]
S2.pi <- F100_pi[,c(1,2,3,24,25,26,27)]

low_C1.pi <- belowAverage(C1.pi)
low_C2.pi <- belowAverage(C2.pi)
low_L1.pi <- belowAverage(L1.pi)
low_L2.pi <- belowAverage(L2.pi)
low_S1.pi <- belowAverage(S1.pi)
low_S2.pi <- belowAverage(S2.pi)

middleChr <- function(ddat, chrom) {
  a <- dim(subset(ddat, chrom == "X"))[1]
  b <- dim(subset(ddat, chrom == "2L"))[1]
  c <- dim(subset(ddat, chrom == "2R"))[1]
  d <- dim(subset(ddat, chrom == "3L"))[1]
  e <- dim(subset(ddat, chrom == "3R"))[1]
  #f <- dim(subset(dat, chr == "4"))[1]
  
  result <- as.vector(rep(NA, 6))
  result[1] <- a/2
  result [2] <- a + (b/2)
  result[3] <- a + b + (c/2)
  result[4] <- a + b + c + (d/2)
  result[5] <- a + b + c + d + (e/2)
  #result[6] <- a + b + c + d + e + (f/2)
  return(result)
  
}

manhattan <- function(ddat) {
  colnames(ddat) <- c("chrom","start","end","fst")
  
  ddat$fst <- pmax(ddat$fst, 0)
  
  ddat <- ddat[ddat$chrom != "4",]
  
  ddat$position <- 1:nrow(ddat)
  
  chrlabel<-middleChr(ddat)
  
  plot <- ggplot(ddat, aes(x=position,y=fst, colour = chrom)) +
    geom_point(size=0.65, show.legend = F, alpha = 0.2) +
    theme(panel.background = element_blank()) +
    scale_colour_manual(values=c("black", "darkgrey", "black", "darkgrey", "black")) +
    stat_smooth(aes(group=as.factor(chrom)), colour = "red", size = 0.5, span= 1, method = "gam", method.args = list(family = binomial)) +
    ylab("Genetic Differentiation (Fst)") +
    scale_x_discrete(limits=c(chrlabel),
                     labels = c("2L", "2R", '3L', '3R', "X")) +
    #geom_smooth(colour = chrom) +
    xlab("Chromosome") +
    ylim(0,1) +
    theme(text = element_text(size=9),
          axis.text.x= element_text(size=7),
          axis.text.y= element_text(size=7))
  
  plot
  
}

manhattan(L1_S1)
manhattan(L1_S2)
manhattan(L2_S1)
manhattan(L2_S2)




