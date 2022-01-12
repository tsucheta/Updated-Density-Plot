library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(png)
library(BiocGenerics)
library(ggplot2)
library(rgl)
library("grid")
library("gridExtra")
library(fields)
library(EBImage)
source("filled.contour3.R")
source("getFeat2b.R")

# Redone by Sucheta Oct 2021
# Reading FIR data created separately for oomycetes paper


FIRdata<-read.csv(file="buscoFIR",sep="\t")
FIRdata

NumBins=40

if ((max(FIRdata$five, na.rm=TRUE)>max(FIRdata$three, na.rm=TRUE)) == TRUE) {
  FIR2Bin<-FIRdata$five
} else {
  FIR2Bin<-FIRdata$three
}


FIR2Bin=FIR2Bin[which(FIR2Bin!=0)]

FIR2Bin<-na.omit(FIR2Bin)

BinSteps<-round(length(FIR2Bin)/(NumBins-1), digits=0)

FIR2BinOrd<-sort(FIR2Bin)

TempBinLimits<-FIR2BinOrd[seq(FIR2BinOrd[2*BinSteps], length(FIR2BinOrd),BinSteps)]

TempBinLimits[length(TempBinLimits)+1]<-max(FIR2Bin, na.rm=TRUE)

x<-seq(length(TempBinLimits))
x

fit<-nls(log(TempBinLimits) ~ a*x + b, start=c(a=0, b=0), 
         algorithm='port', weights=((x-0.5*NumBins)^2))

pred=predict(fit, x)

BinLimits=c(1, round(exp(pred),0), max(FIR2Bin))



xbin=cut(FIRdata$five, breaks= c(BinLimits))


ybin=cut(FIRdata$three, breaks= c(BinLimits))
FIRdata<-cbind(FIRdata, xbin, ybin, 
               genevalue=rep(1, length(FIRdata$five)))
# Here added FUN=sum. In this version it was complaining about the function
GenValMatrix<-with(FIRdata, tapply(genevalue, list(xbin, ybin), FUN=sum))

x<-1:ncol(GenValMatrix)
y<-1:nrow(GenValMatrix)

zlim = range(as.numeric(unlist(GenValMatrix)), finite=TRUE)

mypalette <- colorRampPalette(c( "white", 
                                 "darkblue", 
                                 "forestgreen", 
                                 "goldenrod1", 
                                 "orangered", 
                                 "red3", 
                                 "darkred"), 
                              space="rgb")

mycol=mypalette(2*max(GenValMatrix, na.rm=TRUE))


mycol=mypalette(2*max(GenValMatrix, na.rm=TRUE))

mylabels <- paste(BinLimits[1:length(BinLimits)-1], 
                  BinLimits[2:length(BinLimits)], 
                  sep = " - ", collapse = NULL)

filled.contour3(x, y, z=GenValMatrix, 
                plot.title = title(main ="Busco genes of all Phytophthora genome", 
                                   xlab = "Five prime intergenic regions", 
                                   ylab = "Three prime intergenic regions", 
                                   cex.main=0.8, cex.lab=0.5), 
                key.title = title(main ="Number of genes", cex.main=0.5, line=1), 
                col=mycol, 
                levels = pretty(zlim, 2*max(GenValMatrix, na.rm=TRUE)), 
                plot.axes={axis(1,at=x, labels=mylabels, las=2, 
                                cex.axis=0.5);
                  axis(2,at=y, labels=mylabels, cex.axis=0.5)})




image_name <- paste(as.character(format(Sys.time(), "%Y%m%d%H%M%S")), "_graph", sep="")

png(filename = paste("test", ".png", sep=""))
par(mar=c(0,0,0,0))

filled.contour3(x, y, z=GenValMatrix, 
                col=mycol, 
                levels = pretty(zlim, 2*max(GenValMatrix, na.rm=TRUE)), 
                frame.plot = FALSE, 
                axes = FALSE)

dev.off()


# Read the Image again for overlaying with RXLRs

img <- readPNG(paste("test",".png", sep=""))
g <- rasterGrob(img, interpolate=TRUE)
rxlrdata<-as.data.frame(read.csv(file="rxlrFIR", sep = "\t", header=TRUE))
rxlrdata

NumBins<-as.numeric(c(length(BinLimits)))

gg <-(ggplot(data=rxlrdata, aes(x=five , y=three, geom="blank")) + 
        annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
        coord_fixed(ratio=1) + 
        labs(fill = "RXLR Genes") +
        geom_point(shape=21, colour="red", size=0.5, alpha=0.9, na.rm=FALSE , position=position_jitter(width=0.01, height=0.01)) +
        scale_fill_manual(values = c('rxlr'="red")) +
        scale_y_log10(breaks=BinLimits[2:length(BinLimits)], limits=c(BinLimits[2], BinLimits[NumBins-1])) + 
        scale_x_log10(breaks=BinLimits[2:length(BinLimits)], limits=c(BinLimits[2], BinLimits[NumBins-1])) + 
        ggtitle("RXLR Gene Intergenic Flanking Distance \n vs Global Intergenic Distance ") +
        xlab("Five prime intergenic distance (nt)") + 
        ylab("Three prime intergenic distance (nt)") +
        theme( plot.title = element_text(size=14, face="bold.italic",hjust = 0.5),
               axis.text.y  = element_text(size =8, vjust=0.5),
               axis.text.x  = element_text(size=8, vjust=0.5, angle=90),
               axis.title.x = element_text(face="bold",size=12),
               axis.title.y = element_text(face="bold",size=12))
)

# Run
gg 
