library("matrixStats")

inputFile_WT1 <- "/Volumes/SEAGATE/Genewiz_no1/bow2/WT_1st_rD.bedgraph"
inputFile_WT2 <- "/Volumes/SEAGATE/Genewiz_no2/BS1602266/bow2/WT_2nd_rD.bedgraph"
inputFile_WT3 <- "/Volumes/SEAGATE/Genewiz_no2/BS1602266/bow2/WT_3rd_rD.bedgraph"
inputFile_WT4 <- "/Volumes/SEAGATE/Genewiz_no3/alighned/WT4_rD.bedgraph"

inputFile_lem2D_1 <- "/Volumes/SEAGATE/Genewiz_no3/alighned/lem2_1_rD.bedgraph"
inputFile_lem2D_2 <- "/Volumes/SEAGATE/Genewiz_no3/alighned/lem2_1_rD.bedgraph"

inputFile_14512 <- "/Volumes/SEAGATE/Genewiz_no1/bow2/bqt4D_14512_rD.bedgraph"
inputFile_11935 <- "/Volumes/SEAGATE/Genewiz_no1/bow2/bqt4D_11935_rD.bedgraph"
inputFile_14513 <- "/Volumes/SEAGATE/Genewiz_no1/bow2/bqt4D_14513_rD.bedgraph"
inputFile_14514 <- "/Volumes/SEAGATE/Genewiz_no3/alighned/bqt4D_14514_rD.bedgraph"

#inputFile_man1swi71 <- "/Users/ebrahimih/Desktop/BS1602266/bow2/man1_swi7.bedgraph"
#inputFile_man1swi72 <- "/Users/ebrahimih/Desktop/bow2/man1_swi7_11950.bedgraph"
#inputFile_man1swi73 <- "/Users/ebrahimih/Desktop/bow2/man1_swi7_14519.bedgraph"

#inputFile_man1rap11 <- "/Users/ebrahimih/Desktop/BS1602266/bow2/rap1_man1.bedgraph"
#inputFile_man1rap1bqt41 <- "/Users/ebrahimih/Desktop/BS1602266/bow2/12160_1.bedgraph"
#inputFile_man1rap1bqt42 <- "/Users/ebrahimih/Desktop/BS1602266/bow2/12161_1.bedgraph"
#inputFile_man1swi7bqt41 <- "/Users/ebrahimih/Desktop/bow2/man1_swi7_bqt4D_14517.bedgraph"
#inputFile_man1swi7bqt42 <- "/Users/ebrahimih/Desktop/bow2/man1_swi7_bqt4D_14516.bedgraph"
inputFile_dcr1 <- "/Volumes/SEAGATE/Genewiz_no2/BS1602266/bow2/dcr1D_rD.bedgraph"
inputFile_dcr1bqt4 <- "/Volumes/SEAGATE/Genewiz_no2/BS1602266/bow2/dcr1Dbqt4D_rD.bedgraph"
inputFile_dcr1lem2 <- "/Volumes/SEAGATE/Genewiz_no3/alighned/dcr1Dlem2_rD.bedgraph"

require(ggplot2)

slidingwindowplot <- function(windowsize, inputseq)
{
	starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
	n <- length(starts)
	chunkbps <- numeric(n)
	chunkstats<- numeric(n)
	for (i in 1:n) {
		chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
		chunkmean <- mean(chunk)
		chunkstdv<-sd(chunk)
		chunkbps[i] <- chunkmean
		chunkstats[i]<-chunkstdv
	}
	return (list(starts,chunkbps,chunkstats))
}


#filelist <- c("inputFile_WT1","inputFile_WT2","inputFile_WT3","inputFile_14512","inputFile_11935","inputFile_14513")
filelist <- c("inputFile_WT1","inputFile_WT2","inputFile_WT3","inputFile_WT4","inputFile_dcr1")
#filelist <- c("inputFile_WT1","inputFile_WT2","inputFile_WT3","inputFile_man1rap11")
#filelist <- c("inputFile_WT1","inputFile_WT2","inputFile_WT3","inputFile_man1swi71","inputFile_man1swi72","inputFile_man1swi73")
#filelist <- c("inputFile_WT1","inputFile_WT2","inputFile_WT3","inputFile_man1rap1bqt41","inputFile_man1rap1bqt42")
#filelist <- c("inputFile_WT1","inputFile_WT2","inputFile_WT3","inputFile_dcr1","inputFile_dcr1bqt4")
binSize<-250
maxy<-3
#column.types <- c("character", "numeric", "numeric")
#chrom <- "II"
#chromstart <- 0
#chromend <- 30000
#chrom <- "I"
#chromstart <- 5520000
#chromend <- 5579133
#chrom <- "I"
#chromstart <- 0
#chromend <- 50000
#chrom <- "II"
#chromstart <- 4494964
#chromend <- 4554128
## centromere I
#chrom <- "I"
#chromstart <- 3743687
#chromend <- 3799421
## centromere III
chrom <- "III"
chromstart <- 1055790
chromend <- 1152207
## MAT locus
#chrom <- "II"
#chromstart <- 2127509
#chromend <- 2150471
## ade6
#chrom <- "III"
#chromstart <- 1313650
#chromend <- 1320000
##mei4
#chrom <- "II"
#chromstart <- 1471625
#chromend <- 1477785

for (file in filelist) {
	file.data <- read.csv(file = get(file), head=FALSE, sep="\t",numerals = "no.loss")
	colnames(file.data )<-c("chromosome","position","value")
	signaltrack <- subset(file.data, chromosome == chrom & position > chromstart & position < chromend )
    coverage_signaltrack <- signaltrack[,3]
    assign(paste(file,'.vector', sep = ""), coverage_signaltrack)
    results <- slidingwindowplot(binSize,coverage_signaltrack)
    assign(paste(file,'.window', sep = ""), results)
}

df_WT <- data.frame(inputFile_WT1.window[[1]],apply(cbind(inputFile_WT1.window[[2]],inputFile_WT2.window[[2]],inputFile_WT3.window[[2]],inputFile_WT4.window[[2]]),1,mean),apply(cbind(inputFile_WT1.window[[2]],inputFile_WT2.window[[2]],inputFile_WT3.window[[2]],inputFile_WT4.window[[2]]),1,max),apply(cbind(inputFile_WT1.window[[2]],inputFile_WT2.window[[2]],inputFile_WT3.window[[2]],inputFile_WT4.window[[2]]),1,min))
#df_WT4 <- data.frame(inputFile_WT4.window[[1]],apply(cbind(inputFile_WT4.window[[2]]),1,mean),apply(cbind(inputFile_WT4.window[[2]]),1,max),apply(cbind(inputFile_WT4.window[[2]]),1,min))
#df_bqt4 <- data.frame(inputFile_14512.window[[1]],apply(cbind(inputFile_14512.window[[2]],inputFile_11935.window[[2]],inputFile_14513.window[[2]],inputFile_14514.window[[2]]),1,mean),apply(cbind(inputFile_14512.window[[2]],inputFile_11935.window[[2]],inputFile_14513.window[[2]],inputFile_14514.window[[2]]),1,max),apply(cbind(inputFile_14512.window[[2]],inputFile_11935.window[[2]],inputFile_14513.window[[2]],inputFile_14514.window[[2]]),1,min))
#df_lem2D<- data.frame(inputFile_lem2D_1.window[[1]],apply(cbind(inputFile_lem2D_1.window[[2]],inputFile_lem2D_2.window[[2]]),1,mean),apply(cbind(inputFile_lem2D_1.window[[2]],inputFile_lem2D_2.window[[2]]),1,max),apply(cbind(inputFile_lem2D_1.window[[2]],inputFile_lem2D_2.window[[2]]),1,min))
#df_man1swi7 <- data.frame(inputFile_man1swi71.window[[1]],apply(cbind(inputFile_man1swi71.window[[2]],inputFile_man1swi72.window[[2]],inputFile_man1swi73.window[[2]]),1,mean),apply(cbind(inputFile_man1swi71.window[[2]],inputFile_man1swi72.window[[2]],inputFile_man1swi73.window[[2]]),1,max),apply(cbind(inputFile_man1swi71.window[[2]],inputFile_man1swi72.window[[2]],inputFile_man1swi73.window[[2]]),1,min))
#df_man1rap11 <- data.frame(inputFile_man1rap11.window[[1]],apply(cbind(inputFile_man1rap11.window[[2]],inputFile_man1rap11.window[[2]],inputFile_man1rap11.window[[2]]),1,mean),apply(cbind(inputFile_man1rap11.window[[2]],inputFile_man1rap11.window[[2]],inputFile_man1rap11.window[[2]]),1,max),apply(cbind(inputFile_man1rap11.window[[2]],inputFile_man1rap11.window[[2]],inputFile_man1rap11.window[[2]]),1,min))
#df_man1rap1bqt4 <- data.frame(inputFile_man1rap1bqt41.window[[1]],apply(cbind(inputFile_man1rap1bqt41.window[[2]],inputFile_man1rap1bqt42.window[[2]]),1,mean),apply(cbind(inputFile_man1rap1bqt41.window[[2]],inputFile_man1rap1bqt42.window[[2]]),1,max),apply(cbind(inputFile_man1rap1bqt41.window[[2]],inputFile_man1rap1bqt42.window[[2]]),1,min))
#df_man1swi7bqt4<-data.frame(inputFile_man1swi7bqt41.window[[1]],apply(cbind(inputFile_man1swi7bqt41.window[[2]],inputFile_man1swi7bqt42.window[[2]]),1,mean),apply(cbind(inputFile_man1swi7bqt41.window[[2]],inputFile_man1swi7bqt42.window[[2]]),1,max),apply(cbind(inputFile_man1swi7bqt41.window[[2]],inputFile_man1swi7bqt42.window[[2]]),1,min))
df_dcr1 <- data.frame(inputFile_dcr1.window[[1]],apply(cbind(inputFile_dcr1.window[[2]],inputFile_dcr1.window[[2]]),1,mean),apply(cbind(inputFile_dcr1.window[[2]],inputFile_dcr1.window[[2]]),1,max),apply(cbind(inputFile_dcr1.window[[2]],inputFile_dcr1.window[[2]]),1,min))
#df_dcr1bqt4 <- data.frame(inputFile_dcr1bqt4.window[[1]],apply(cbind(inputFile_dcr1bqt4.window[[2]],inputFile_dcr1bqt4.window[[2]]),1,mean),apply(cbind(inputFile_dcr1bqt4.window[[2]],inputFile_dcr1bqt4.window[[2]]),1,max),apply(cbind(inputFile_dcr1bqt4.window[[2]],inputFile_dcr1bqt4.window[[2]]),1,min))
#df_dcr1lem2 <- data.frame(inputFile_dcr1lem2.window[[1]],apply(cbind(inputFile_dcr1lem2.window[[2]],inputFile_dcr1lem2.window[[2]]),1,mean),apply(cbind(inputFile_dcr1lem2.window[[2]],inputFile_dcr1lem2.window[[2]]),1,max),apply(cbind(inputFile_dcr1lem2.window[[2]],inputFile_dcr1lem2.window[[2]]),1,min))

#colname<-c("x","mean","sd")
colnames(df_WT)<-c("x","mean","max","min")
#colnames(df_bqt4)<-c("x","mean","max","min") #color="red"
#colnames(df_lem2D)<-c("x","mean","max","min") #color="red"
#colnames(df_man1swi7)<-c("x","mean","max","min") # color "darkorchid"
#colnames(df_man1rap11)<-c("x","mean","max","min") # color "goldenrod1"
#colnames(df_man1rap1bqt4)<-c("x","mean","max","min") # color "darkorange3"
#colnames(df_man1swi7bqt4)<-c("x","mean","max","min") # color "forestgreen"
colnames(df_dcr1)<-c("x","mean","max","min") #color="deeppink4"
#colnames(df_dcr1bqt4)<-c("x","mean","max","min") #dimgray
#colnames(df_dcr1lem2)<-c("x","mean","max","min") #dimgray

mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line
eb <- aes(ymax = mean + 0, ymin = mean - 0)
svg("mapping_all.svg",height=80, width=800)
dev.new(width=3, height=3)
ggplot(data = df_WT, aes(x = x, y = mean)) +geom_line(colour="#0066CC",size=1)+geom_ribbon(data = df_WT, aes(ymin=min, ymax=max) , fill="#0066CC", alpha=0.4)+geom_line(data = df_dcr1, aes(x=x,y=mean),color="red")+geom_ribbon(data = df_dcr1, aes(ymin=min, ymax=max) , fill="red", alpha=0.4)+mytheme+xlab("Reference Start Position")+ylab("coverage")+ scale_x_continuous(expand = c(0,0))+scale_y_continuous(limits = c(0, maxy))+ggtitle("Coverage Across Reference")+coord_cartesian(expand = c(0))
dev.off()