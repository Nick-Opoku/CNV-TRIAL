library(HaarSeg)
library(GenomicRanges)

I <- read.table("/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/Removalofzeros.bed", sep="\t",as.is=TRUE)

seq.tab=haarSeg(I$V5) # chromPos = matrix(c(1,length(I)),  nrow = 4, ncol = 4))

logratio = I$V5
a <- plot(logratio)
print(seq.tab$Segmented)
b <- lines(seq.tab$Segmented, col="red", lwd=3)
print(seq.tab$SegmentsTable)

I <- read.table("/home/nick/Documents/TRIAL/CODES-FOR-GRAPHS/TRYS/my-pipeline/Viterbi-Code/my-vit-algorithm-using-for-thesis/z-norm.bed",
                sep="\t",as.is=TRUE)

seq.tab=haarSeg(I$V2) # chromPos = matrix(c(1,length(I)),  nrow = 4, ncol = 4))


Signal = I$V2
a <- plot(Signal, xlab='Genomic region')
print(seq.tab$Segmented)
b <- lines(seq.tab$Segmented, col="red", lwd=3)
b <- print(seq.tab$SegmentsTable)





