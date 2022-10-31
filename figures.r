#Example for one
data<-read.table("/path/to/AVE.fst",header=TRUE)

# OR #

#data<-AVE

# Break everything in to chromosomes again so that two different chromosomes are not included in the same window

data2L<-data[data$ID.Chromosome=="chr2L",]
data2R<-data[data$ID.Chromosome=="chr2R",]
data3L<-data[data$ID.Chromosome=="chr3L",]
data3R<-data[data$ID.Chromosome=="chr3R",]
data4<-data[data$ID.Chromosome=="chr4",]
dataX<-data[data$ID.Chromosome=="chrX",]

# Sliding window function

sliding_window<-function(data,window){
temp.vec=vector("list",0)
results.vec=vector("list",0)
x=0
y=5000
i=1
while (i <= length(data$Means)){
 if (data$ID.Position[i] >= x & data$ID.Position[i] < y){
    temp.vec<-append(temp.vec,data$Means[i])
    i<-i+1
  } else {
    if (data$ID.Position[i] <y) {
      temp.vec<-append(temp.vec,NA)
      i<-i+1
      } else {
    temp.vec<-as.numeric(temp.vec)
    results.vec <- append(results.vec,mean(temp.vec, na.omit=TRUE))
    temp.vec=vector("list",0)
    x <- x + window
    y <- y + window
      }
  }
}
print(results.vec)
}

# Running sliding window function on all chromosomes

fst2L<-sliding_window(data2L,5000)
fst2R<-sliding_window(data2R,5000)
fst3L<-sliding_window(data3L,5000)
fst3R<-sliding_window(data3R,5000)
fst4<-sliding_window(data4,5000)
fstX<-sliding_window(dataX,5000)

# Making matrices which are easier to work with

fst2L <- as.matrix(fst2L)
fst2R <- as.matrix(fst2R)
fst3L <- as.matrix(fst3L)
fst3R <- as.matrix(fst3R)
fst4 <- as.matrix(fst4)
fstX <- as.matrix(fstX)

fst2L <- fst2L[which(fst2L!="NaN"),]
fst2R <- fst2R[which(fst2R!="NaN"),]
fst3L <- fst3L[which(fst3L!="NaN"),]
fst3R <- fst3R[which(fst3R!="NaN"),]
fst4 <- fst4[which(fst4!="NaN"),]
fstX <- fstX[which(fstX!="NaN"),]


# Merging all in to one file

fst2L<-cbind(rep("2L",length(fst2L)),fst2L)
fst2R<-cbind(rep("2R",length(fst2R)),fst2R)
fst3L<-cbind(rep("3L",length(fst3L)),fst3L)
fst3R<-cbind(rep("3R",length(fst3R)),fst3R)
fst4<-cbind(rep("4",length(fst4)),fst4)
fstX<-cbind(rep("X",length(fstX)),fstX)
fst_data<-rbind(fstX,fst2L,fst2R,fst3L,fst3R,fst4)



# Changing column names

crap<-c("chr","fst")
colnames(fst_data)<-crap
fst_data<-as.data.frame(fst_data)

# This code adds numeric indicators for each row to make sure the order in the plots are correct

g <- nrow(fst_data[which(fst_data$chr=='2L'),])
h <- nrow(fst_data[which(fst_data$chr=='2R'),])
i <- nrow(fst_data[which(fst_data$chr=='3L'),])
j <- nrow(fst_data[which(fst_data$chr=='3R'),])
k <- nrow(fst_data[which(fst_data$chr=='4'),])
l <- nrow(fst_data[which(fst_data$chr=='X'),])

fst_data$number <-  c((1:l),
                   (l+1):(l+g), 
                   (l+g+1):(l+g+h), 
                   (l+g+h+1):(l+g+h+i),
                   (l+g+h+i+1):(l+g+h+i+j),
                   (l+g+h+i+j+1):(l+g+h+i+j+k))

fst_data <- as.data.frame(lapply(fst_data, unlist))
fst_data$number<-as.numeric(fst_data$number)
fst_data$chr<-as.factor(fst_data$chr)

# Plotting

plot<-ggplot(fst_data, aes(x=as.numeric(number), y=as.numeric(fst), color=as.factor(chr))) +
  geom_point(size=0.65, show.legend = F, alpha = 0.2) +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values=c("black", "darkgrey", 'black', 'darkgrey', 'black', 'darkgrey')) +
  geom_smooth(aes(group=as.factor(chr)), colour = "red", size = 0.5) +
  ylab("Genomic differentiation (Fst)") +
  xlab("Chromosome") +
  scale_x_discrete(labels = c("X","2L","2R","3L","3R","4")) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))
        
        
# Need to change file name !!!!!!!!!

png("/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVS_window5000.png",type="cairo")
plot
dev.off()
