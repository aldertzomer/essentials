library(edgeR) 
targets <- read.delim(file = "gene_targets.txt", stringsAsFactors = FALSE) 
d <- readDGE(targets, columns=c(2,1))
d <- d[rowSums(d$counts) > 5*nrow(targets), ]
png('gene_MDS_non_normalized.png')
plotMDS.dge(d, main='MDS Plot non normalized on genes', xlim=c(-3,3), ylim=c(-3,3))
dev.off()
d <- calcNormFactors(d, method=c("quantile")) 
d <- estimateCommonDisp(d) 
png('gene_MDS_normalized.png')
plotMDS.dge(d, main='MDS Plot normalized on genes', xlim=c(-3,3), ylim=c(-3,3))
dev.off()
d <- estimateTagwiseDisp(d, prior.n=10) 
de.tagwise <-exactTest(d, common.disp = FALSE) 
toptags_tgw.out <-topTags(de.tagwise, n=nrow(de.tagwise)+1)

gene_minomics <-toptags_tgw.out$table[toptags_tgw.out$table[,1]>-17, ]
gene_minomics <-gene_minomics[, c(2,4)]
write.table(gene_minomics, 'gene_minomics.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
 
write.table(d$samples, 'gene_samples.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
detags200 <- rownames(topTags(de.tagwise, n = 200)$table)
png('gene_plot.png')
plotSmear(d, de.tags = detags200, main = 'Fold change vs signal plot')
abline(h = c(-2, 2), col = 'dodgerblue')
dev.off()
allcounts.merged <- merge(d$counts, d$pseudo.alt, by=0)
row.names(allcounts.merged) <- allcounts.merged[, 1]
allcounts.merged <- allcounts.merged[, 2:(2 * nrow(targets)+1)] 
alldata.merged <- merge(allcounts.merged, toptags_tgw.out$table, by=0)
genefunctions <- read.delim("genome.ptt", header=T) 
alloutput.merged <- merge (alldata.merged, genefunctions, by=1, all=TRUE) 
write.table(alloutput.merged, 'gene_alloutputmerged.tsv', sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE) 
q() 
