#!/home/jsong/R-3.3.1/bin/R

args <- commandArgs(T)
infileL <- args[1]
infileR <- args[2]

path <- strsplit(infileL, "/")
prefix <- strsplit(rev(path[[1]])[1], ".", fixed=T)[[1]][1]

dataL <- read.table(infileL, sep="\t", header=T)
dataR <- read.table(infileR, sep="\t", header=T)
data <- merge(dataL, dataR, by="ReadLen")
offset <- dataL[2:dim(dataL)[1], 2:dim(dataL)[2]]
rownames(offset) <- dataL[2:dim(dataL)[1],1]
colnames(offset) <- seq(-20,-1)
offsetM <- as.matrix(offset)

library(pheatmap)
pdf(paste('Plot/', prefix, 'metagene_HP.pdf',sep=""), onefile=F, paper="a4")
pheatmap(offsetM, cluster_rows=F, cluster_cols=F, cellheight=10, cellwidth=10, main = paste("Metagene plot of", prefix))
dev.off()

for (i in 1:length(data$ReadLen)){
    size <- data$ReadLen[i]
    Xval <- seq(-20,50)
    Yval <- data[i,-1]
    x <- seq(1,length(Xval),3)
    y <- seq(2,length(Xval),3)
    z <- seq(3,length(Xval),3)
    pdf(paste('Plot/', prefix, '_metagene.', size,'.pdf', sep=''), onefile=F, width=9,height=5)
    plot(Xval,Yval, type="h", col="white",xlim=c(-20,50),xlab="Distance to start codon (nt)", ylab="Number of RPFs", main=paste("Metagene plot of RPFs in",size, "nt"), cex.axis=1.5, cex.lab=1.5)
    lines(Xval[x],Yval[x],col="cyan",type="h",lwd=2)
    lines(Xval[y],Yval[y],col="orange",type="h", lwd=2)
    lines(Xval[z],Yval[z],col="purple",type="h",lwd=2)
    dev.off()
}



