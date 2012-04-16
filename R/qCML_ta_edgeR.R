Args <- commandArgs(TRUE)
options(error=expression(NULL))
library(edgeR) 
targets <- read.delim(file = "ta_targets.txt", stringsAsFactors = FALSE) 
d <- readDGE(targets, columns=c(2,1))
d$counts <- round(d$counts)
#d <- d[rowSums(d$counts) > 2*nrow(targets), ]
#d$counts <- replace(d$counts, d$counts ==0,1)
d$counts <- d$counts +1

if (nrow(targets) >2){
    p3 <- prcomp(d$counts)
    png('ta_PCA_non_normalized.png', width=1280, height=1024)
    biplot(p3, cex=1:2, col=c(40,1), main='PCA Plot non normalized')
    #plotMDS.dge(d, main='MDS Plot non normalized on ta-sites', xlim=c(-3,3), ylim=c(-3,3))
    dev.off()
} else {
    png('ta_PCA_non_normalized.png', width=1280, height=1024)
    plot(log(d$counts))
    title("three samples needed for PCA. x-y plot generated")
    dev.off()
}
                

if (Args[1] != "Totalcounts") d <- calcNormFactors(d, method=c(Args[1]))
d <- estimateCommonDisp(d) 

if (nrow(targets) >2){
    p3 <- prcomp(d$pseudo.alt)
    png('ta_PCA_normalized.png', width=1280, height=1024)
    biplot(p3, cex=1:2, col=c(40,1), main='PCA Plot normalized')
    #plotMDS.dge(d, main='MDS Plot normalized on ta-sites', xlim=c(-3,3), ylim=c(-3,3))
    dev.off()                    
} else {
    png('ta_PCA_normalized.png', width=1280, height=1024)
    plot(log(d$pseudo.alt))
    title("three samples needed for PCA. x-y plot generated")
    dev.off()
}
                

if (Args[3] == "tagwise"){
    if (nrow(targets) >2){
	d <- estimateTagwiseDisp(d, prior.n=as.numeric(Args[4]), trend=FALSE)
	de.tagwise <-exactTest(d, common.disp = FALSE)
    } else {de.tagwise <-exactTest(d, common.disp = TRUE)}
}
    
if (Args[3] == "common") de.tagwise <-exactTest(d, common.disp = TRUE)

toptags_tgw.out <-topTags(de.tagwise, n=nrow(de.tagwise)+1, adjust.method=Args[2]) 
write.table(d$samples, 'ta_samples.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
detags200 <- rownames(topTags(de.tagwise, n = 200)$table)
png('ta_plot.png', width=1280, height=1024)
plotSmear(d, de.tags = detags200, main = 'MA plot on ta sites')
abline(h = c(-2, 2), col = 'dodgerblue')
dev.off()
d$counts <- d$counts-1

samplenames <- read.delim("samplenames.txt", header=F)
samplenumber <- 1
while (samplenumber <= nrow(samplenames)) {
     colnames(d$counts)[samplenumber] <-  as.matrix(samplenames)[samplenumber]
     colnames(d$pseudo.alt)[samplenumber] <-  as.matrix(samplenames)[samplenumber]
     samplenumber <- samplenumber + 1
}

allcounts.merged <- merge(d$counts, round(d$pseudo.alt), by=0)
row.names(allcounts.merged) <- allcounts.merged[, 1]
allcounts.merged <- allcounts.merged[, 2:(2 * nrow(targets)+1)]
alldata.merged <- merge(allcounts.merged, toptags_tgw.out$table, by=0)
genefunctions <- read.delim("ta.ptt", header=F) 

name <- as.numeric(as.matrix(row.names(as.data.frame(d$conc))))
start <- abs(as.numeric(as.matrix(row.names(as.data.frame(d$conc)))))
stop <- abs(as.numeric(as.matrix(row.names(as.data.frame(d$conc)))))

forward <-c()
for (i in name) {
    if(i >= 0) forward <- c(forward, 1)
    if(i < 0) forward <- c(forward, 0)
}

merge1 <- merge(name,start, by=0)
rownames(merge1) <-merge1[, 1]
merge1 <- merge1[c(2,3)]
merge1 <- merge(merge1, stop, by=0)
rownames(merge1) <-merge1[, 1]
merge1 <- merge1[c(2:4)]
merge1 <- merge(merge1, forward, by=0)
rownames(merge1) <-merge1[, 1]
merge1 <- merge1[c(2:5)]
rownames(merge1) <-merge1[, 1]
controlsequence <- round((as.data.frame(d$conc)[c(2)] *  mean(colSums(d$pseudo.alt))))
targetsequence <- round((as.data.frame(d$conc)[c(3)] *  mean(colSums(d$pseudo.alt))))
merge1[ , "type"] <- "promoter"
controlsites <- merge(merge1,controlsequence, by=0)
targetsites <- merge(merge1,targetsequence, by=0)
controlsites <- controlsites[c(2:7)]
controlsites <- controlsites[controlsites[,6]>Args[5], ]
targetsites <- targetsites[c(2:7)]
targetsites <- targetsites[targetsites[,6]>Args[5], ]
controlsites[ , "0"] <- "0"
controlsites[ , "1"] <- "0"
controlsites[ , "2"] <- "0"
controlsites[ , "3"] <- "0"
targetsites[ , "0"] <- "0"
targetsites[ , "1"] <- "0"
targetsites[ , "2"] <- "0"
targetsites[ , "3"] <- "0"
colnames(controlsites)[1] <- "#name"
colnames(controlsites)[2] <- "start"
colnames(controlsites)[3] <- "stop"
colnames(controlsites)[4] <- "forward"
colnames(controlsites)[5] <- "type"
colnames(controlsites)[6] <- "sequence"
colnames(controlsites)[7] <- "0"
colnames(controlsites)[8] <- "0"
colnames(controlsites)[9] <- "0"
colnames(controlsites)[10] <- "0"
colnames(targetsites)[1] <- "#name"
colnames(targetsites)[2] <- "start"
colnames(targetsites)[3] <- "stop"
colnames(targetsites)[4] <- "forward"
colnames(targetsites)[5] <- "type"
colnames(targetsites)[6] <- "sequence"
colnames(targetsites)[7] <- "0"
colnames(targetsites)[8] <- "0"
colnames(targetsites)[9] <- "0"
colnames(targetsites)[10] <- "0"
write.table(controlsites, 'control_minomics.tsv', sep = "\t", quote=FALSE, row.names = FALSE)
write.table(targetsites, 'target_minomics.tsv', sep = "\t", quote=FALSE, row.names = FALSE)
alloutput.merged <- merge (alldata.merged, genefunctions, by=1, all.x=TRUE)
alloutput.merged2 <- alloutput.merged[, c(1,1,2:ncol(alloutput.merged))]
alloutput.merged2$Row.names.1 <- abs(as.numeric(alloutput.merged2$Row.names.1))
colnames(alloutput.merged2)[2] <- "Position"
write.table(alloutput.merged2, 'ta_alloutputmerged.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
    
    #density plot and possible cutoff calculation
    library(zoo)

    nr.x.minima <- 100
    densityadjust <- 0.8

    while(nr.x.minima>4) {
        d<-density(toptags_tgw.out$table$logFC, kernel = c("gaussian"), bw="nrd", n=2048, adjust=densityadjust)
        y <- d$y
        x <- d$x
        xz <- as.zoo(y)
        rxzmin <- rollapply(xz, 3, function(x) which.min(x)==2) #local minima
        rxzmax <- rollapply(xz, 3, function(x) which.max(x)==2) #local maxima
        x.minima<- x[index(rxzmin)[coredata(rxzmin)]]
        y.minima<- y[index(rxzmin)[coredata(rxzmin)]]
        nr.x.minima=nrow(as.matrix(x.minima))
        densityadjust <- densityadjust+0.1
    }

    y.range <- max(y)-min(y)
    y.shift <- y.range/20
    png('ta_densityplot.png', width=1280, height=1024)
    plot(d, main="Density plot of Log Fold change. Red dot depics putative fold change cutoff(s)", xlab="Log2 fold change")
    points(x.minima,y.minima, col="red", pch=19:21 )
    text(x.minima,y.minima+y.shift, labels=round(x.minima,2))
    dev.off()

q()
