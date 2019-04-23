args <- commandArgs(T)
infile <- args[1]
path <- strsplit(infile, "/")
prefix <- strsplit(rev(path[[1]])[1], ".", fixed=T)[[1]][1]
data <- read.table(infile, sep="\t", header=T)
freqFile <- paste (prefix,"_multitaper.Freq.txt",sep="");
cat("#size","logVmax","freq", file=freqFile, sep="\t",ecol="\n")
library(multitaper)
for(i in 1:dim(data)[1]) {
	size <- data[i,1]
	ts <- ts(as.numeric(data[i, 2:dim(data)[2]]))
	spec <- spec.mtm(ts, plot=F,nFFT=500, Ftest=T,nw=4,k=7)
	Pval <- (1-pf(spec$mtm$Ftest, 2, 12, lower.tail=T))
	logV <- -1*log10(Pval)
	logV[is.na(logV)] <- 0
	logVmax <- round(max(logV),2)
	freqMax <- round(spec$freq[which.max(logV)],2)
	freqMax[is.na(freqMax)] <- 0
	cat(size,logVmax,freqMax, file=freqFile,append=T,sep="\t",ecol="\n")
	if (logVmax>=0 && freqMax>=0){
		outpdf <- paste('Plot/',prefix,"_multitaper",".",size,".pdf", sep="")
		pdf(outpdf, onefile=F, paper="a4")
		plot(NA, NA, xlab="Frequency (Hz)", ylab="-log10(P-value)", cex.lab=1.5, cex.axis=2, main=paste("Frequency of ", size, "-nt reads"), xlim=c(0,0.5), ylim=c(0,12))
		for (j in 1:dim(data)[1]){
			spec2 <- spec.mtm(ts(as.numeric(data[j, 2:dim(data)[2]])), plot=F,nFFT=500, Ftest=T,nw=4,k=7)
			Pval2 <- (1-pf(spec2$mtm$Ftest, 2, 12, lower.tail=T))
			logV2 <- -1*log10(Pval2)
			lines(spec2$freq, logV2, lwd=2,type="l",col=rgb(0.7,0.7,0.7))
		}
		lines(spec$freq, logV, lwd=3, type="l", col="red")
		abline(v=0.3333, col="blue", lwd=4, lty=3)
		abline(h=2, col="orange",lwd=2)
		legend ("topleft",c(paste(size, "bp", sep=" "), "other sizes"), lwd=3, cex=2, bty="n",col=c("red", rgb(0.7,0.7,0.7)), lty=c(1,1)) 
		dev.off()
	}
}

