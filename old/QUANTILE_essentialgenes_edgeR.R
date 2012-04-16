
    #load all count data with ta sites included

    library(edgeR)
    targets <- read.delim(file = "essentialgenes_targets.txt", stringsAsFactors = FALSE)
    d <- readDGE(targets, columns=c(2,1))
    d <- d[rowSums(d$counts) > 5*nrow(targets), ]
    png('essentialgenes_MDS_non_normalized.png')
    plotMDS.dge(d, main='MDS Plot non normalized on measured and expected counts', xlim=c(-3,3), ylim=c(-3,3))
    dev.off()

    #normalization changes ta sites to expected frequencies

    d <- calcNormFactors(d, refColumn=nrow(targets)-1, method=c("quantile"))
    d <- estimateCommonDisp(d)
    png('essentialgenes_MDS_normalized.png')
    plotMDS.dge(d, main='MDS Plot normalized on measured and expected counts', xlim=c(-3,3), ylim=c(-3,3))
    dev.off()


    #test for differences between measured and expected frequencies
    d <- estimateTagwiseDisp(d, prior.n=10)
    de.tagwise <-exactTest(d, common.disp = FALSE)
    toptags_tgw.out <-topTags(de.tagwise, n=nrow(de.tagwise)+1)
    
    #generate output
    write.table(d$samples, 'essentialgenes_samples.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

    essentialgenes_minomics <-toptags_tgw.out$table[toptags_tgw.out$table[,1]>-17, ]
    essentialgenes_minomics <-essentialgenes_minomics[, c(2,4)]
    write.table(essentialgenes_minomics, 'essentialgenes_minomics.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

    write.table(d$counts, 'essentialgenes_data.tsv', sep = "\t")
    detags200 <- rownames(topTags(de.tagwise, n = 200)$table)
    png('essentialgenes_plot.png')
    plotSmear(d, de.tags = detags200, main = 'Fold change vs signal plot on measured and expected counts')
    abline(h = c(-2, 2), col = 'dodgerblue')
    dev.off()
    
    allcounts.merged <- merge(d$counts, d$pseudo.alt, by=0)
    row.names(allcounts.merged) <- allcounts.merged[, 1]
    allcounts.merged <- allcounts.merged[, 2:(2 * nrow(targets)+1)]
    alldata.merged <- merge(allcounts.merged, toptags_tgw.out$table, by=0)

    genefunctions <- read.delim("genome.ptt", header=T)
    alloutput.merged <- merge (alldata.merged, genefunctions, by=1, all=TRUE)
    write.table(alloutput.merged, 'essentialgenes_alloutputmerged.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)


    #spatial normalization using loess. quite tricky
    #load locations of genes and couple this to fold change data
    genelocations <- read.delim("genome.ptt", header=T)
    row.names(genelocations)<-(genelocations[,(1)])
    fcandlocations <- merge (toptags_tgw.out$table, genelocations, by=0)
    row.names(fcandlocations) <- fcandlocations[, 1]
    attach(fcandlocations)
    sort1.fcandlocations <- fcandlocations[order(Start) , ]

    #actual loess normalization
    x <-sort1.fcandlocations[,(7)]
    y <-sort1.fcandlocations[,(3)]
    detach(fcandlocations)
    y.loess <- loess(y ~ x, span=1, data.frame(x=x, y=y))
    y.predict <- predict(y.loess, data.frame(x=x))
    y.ratio <- 2^y.predict
    y.ratio <- as.data.frame(y.ratio)
    row.names(y.ratio) <- row.names(sort1.fcandlocations)
    raw.data <- read.delim("essentialgenes_data.tsv")
    dataratio.merged<- merge(raw.data,y.ratio, by=0)

    #this loop normalizes the raw data counts
    i <- 2
    while(i<ncol(dataratio.merged)-1) {
         dataratio.merged[c(i)] <- dataratio.merged[c(i)] / dataratio.merged$y.ratio
         i <- i+1
    }

    # do the statistical test again, but on loess normalized data, generate output etc
    d <- dataratio.merged[, 3:ncol(dataratio.merged)-1]
    rownames(d) <- dataratio.merged[, 1]
    group <- c(rep("measured", ncol(dataratio.merged)-3), rep("expected", 1))
    d <- DGEList(counts = d, group = group)
    d <- d[rowSums(d$counts) > 5*nrow(d$samples), ]
    png('loess_essentialgenes_MDS_non_normalized.png')
    plotMDS.dge(d, main='MDS Plot non normalized on measured and expected counts', xlim=c(-3,3), ylim=c(-3,3))
    dev.off()
    d <- calcNormFactors(d, refColumn=nrow(d$samples), method=c("quantile"))
    d <- estimateCommonDisp(d)
    png('loess_essentialgenes_MDS_normalized.png')
    plotMDS.dge(d, main='MDS Plot normalized on measured and expected counts', xlim=c(-3,3), ylim=c(-3,3))
    dev.off()
    d <- estimateTagwiseDisp(d, prior.n=10)
    de.tagwise <-exactTest(d, common.disp = FALSE)
    toptags_tgw.out <-topTags(de.tagwise, n=nrow(de.tagwise)+1)
   
    loess_essentialgenes_minomics <-toptags_tgw.out$table[toptags_tgw.out$table[,1]>-17, ]
    loess_essentialgenes_minomics <-loess_essentialgenes_minomics[, c(2,4)]
    write.table(loess_essentialgenes_minomics, 'loess_essentialgenes_minomics.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

    write.table(d$samples, 'loess_essentialgenes_samples.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
    detags200 <- rownames(topTags(de.tagwise, n = 200)$table)
    png('loess_essentialgenes_plot.png')
    plotSmear(d, de.tags = detags200, main = 'Fold change vs signal plot on measured and expected counts')
    abline(h = c(-2, 2), col = 'dodgerblue')
    dev.off()
    allcounts.merged <- merge(d$counts, d$pseudo.alt, by=0)
    row.names(allcounts.merged) <- allcounts.merged[, 1]
    allcounts.merged <- allcounts.merged[, 2:(2 * nrow(targets)+1)]
    alldata.merged <- merge(allcounts.merged, toptags_tgw.out$table, by=0)
    genefunctions <- read.delim("genome.ptt", header=T)
    alloutput.merged <- merge (alldata.merged, genefunctions, by=1, all=TRUE)
    write.table(alloutput.merged, 'loess_essentialgenes_alloutputmerged.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
    q()
