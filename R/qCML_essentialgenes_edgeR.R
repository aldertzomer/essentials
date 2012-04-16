    Args <- commandArgs(TRUE)
    options(error=expression(NULL))
    #load all count data with ta sites included

    library(edgeR)
    targets <- read.delim(file = "essentialgenes_targets.txt", stringsAsFactors = FALSE)
    d <- readDGE(targets, columns=c(2,1))
    
    d <- d[rowSums(d$counts) > 10, ]
    d$counts <- round(d$counts)
    counts.old<-d$counts
    d$counts <- d$counts + 1
    
    if (nrow(targets) >2){
        p3 <- prcomp(d$counts)
        png('essentialgenes_PCA_non_normalized.png', width=1280, height=1024)
        biplot(p3, cex=1:2, col=c(40,1), main='PCA Plot on measured and expected counts')
        #plotMDS.dge(d, main='MDS Plot non normalized on measured and expected counts', xlim=c(-3,3), ylim=c(-3,3))
        dev.off()
    } else {
        png('essentialgenes_PCA_non_normalized.png', width=1280, height=1024)
        plot(log(d$counts))
        title("three samples needed for PCA. x-y plot generated")
        dev.off()
    }


    #median normalization to 2^10

    #channel.medians<-apply(log2(d$counts),2,median)
    #channel.norm.factor<-channel.medians-10
    #normalized.log.x=sweep(log2(d$counts),2,channel.norm.factor)
    #d$counts <- 2^normalized.log.x
    #d$counts <- round(d$counts)
    
    #normalization changes ta sites to expected frequencies

    if (Args[1] != "Totalcounts") d <- calcNormFactors(d, refColumn=nrow(targets)-1, method=c(Args[1]))
       
    d <- estimateCommonDisp(d)
    
    if (nrow(targets) >2){
	p3 <- prcomp(d$pseudo.alt)
	png('essentialgenes_PCA_normalized.png', width=1280, height=1024)
	biplot(p3, cex=1:2, col=c(40,1), main='PCA Plot on normalized measured and expected counts')
	#plotMDS.dge(d, main='PCA Plot normalized on measured and expected counts', xlim=c(-3,3), ylim=c(-3,3))
	dev.off()
    } else {
        png('essentialgenes_PCA_normalized.png', width=1280, height=1024)
        plot(log(d$pseudo.alt))
        title("three samples needed for PCA. x-y plot generated")
        dev.off()
    }
	                

    #test for differences between measured and expected frequencies

    if (Args[3] == "tagwise"){ 
	if (nrow(targets) >2){    
	    d <- estimateTagwiseDisp(d, prior.n=as.numeric(Args[4]), trend=TRUE)
	    de.tagwise <-exactTest(d, common.disp = FALSE)
	} else {de.tagwise <-exactTest(d, common.disp = TRUE)}
    }
    
    if (Args[3] == "common") de.tagwise <-exactTest(d, common.disp = TRUE)

    toptags_tgw.out <-topTags(de.tagwise, n=nrow(de.tagwise)+1, adjust.method=Args[2])
    
    #generate output

    essentialgenes_minomics <-toptags_tgw.out$table[toptags_tgw.out$table[,1]>-30, ]
    colnames(essentialgenes_minomics)[2] <- "Essential genes Log2 fold change"
    colnames(essentialgenes_minomics)[4] <- "Essential genes FDR corrected P"
    essentialgenes_minomics <-essentialgenes_minomics[, c(2,4)]
    write.table(essentialgenes_minomics, 'essentialgenes_minomics.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

    write.table(d$samples, 'essentialgenes_samples.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
    write.table(d$counts, 'essentialgenes_data.tsv', sep = "\t")
    detags200 <- rownames(topTags(de.tagwise, n = 200)$table)
    png('essentialgenes_plot.png', width=1280, height=1024)
    plotSmear(d, de.tags = detags200, main = 'Fold change vs signal plot on measured and expected counts')
    abline(h = c(-2, 2), col = 'dodgerblue')
    dev.off()
    d$counts<-counts.old
    
    samplenames <- read.delim("samplenames.txt", header=F)
    samplenumber <- 1
    while (samplenumber <= nrow(samplenames)) {
	colnames(d$counts)[samplenumber] <-  as.matrix(samplenames)[samplenumber]
	colnames(d$pseudo.alt)[samplenumber] <-  as.matrix(samplenames)[samplenumber]
	samplenumber <- samplenumber + 1
    }
    colnames(d$counts)[(nrow(samplenames)+1)] <-  "unique_flanking_sequences"
    colnames(d$pseudo.alt)[(nrow(samplenames)+1)] <-  "expected_reads"
    
    allcounts.merged <- merge(d$counts, round(d$pseudo.alt), by=0)
    row.names(allcounts.merged) <- allcounts.merged[, 1]
    allcounts.merged <- allcounts.merged[, 2:(2 * nrow(targets)+1)]
    alldata.merged <- merge(allcounts.merged, toptags_tgw.out$table, by=0)

    genefunctions <- read.delim("genome.ptt", header=T)
    alloutput.merged <- merge (alldata.merged, genefunctions, by=1, all=TRUE)
    write.table(alloutput.merged, 'essentialgenes_alloutputmerged.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

    #density plot and possible cutoff calculation
    library(zoo)
    
    nr.x.minima <- 100
    densityadjust <- 0.1
    while(nr.x.minima>1) { 
	d.density<-density(toptags_tgw.out$table$logFC, kernel = c("gaussian"), bw="nrd", n=2048, adjust=densityadjust)
	y <- d.density$y
	x <- d.density$x
	xz <- as.zoo(y)
	rxzmin <- rollapply(xz, 3, function(x) which.min(x)==2) #local minima
	rxzmax <- rollapply(xz, 3, function(x) which.max(x)==2) #local maxima
	x.minima<- x[index(rxzmin)[coredata(rxzmin)]]
	y.minima<- y[index(rxzmin)[coredata(rxzmin)]]
	nr.x.minima=nrow(as.matrix(x.minima[x.minima< -0.5]))
	densityadjust <- densityadjust+0.1
    }

    
    y.range <- max(y)-min(y)
    y.shift <- y.range/20
    png('essentialgenes_densityplot.png', width=1280, height=1024)
    plot(d.density, main="Density plot of Log Fold change. Red dot depics putative fold change cutoff(s)", xlab="Log2 fold change")
    points(x.minima,y.minima, col="red", pch=19:21 )
    text(x.minima,y.minima+y.shift, labels=round(x.minima,2))
    dev.off()
q()
