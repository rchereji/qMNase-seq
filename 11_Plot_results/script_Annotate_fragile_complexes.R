library(rtracklayer)
library(GenomeInfoDb)

## S2_exp1
fragile_complexes <- import("S2_exp1/Dmel_fragile_complexes.bed")

require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
require(ChIPseeker)
require(org.Dm.eg.db)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

peak <- fragile_complexes
seqlevelsStyle(peak) = "UCSC"
peakAnno <- annotatePeak(peak,
                         tssRegion=c(-3000,3000),
                         TxDb = txdb,
                         annoDb="org.Dm.eg.db")

pdf("S2_exp1/Pie_chart_fragile_complexes.pdf",width=7,height=3,paper='special') 
plotAnnoPie(peakAnno)
dev.off()

## S2_exp2
fragile_complexes <- import("S2_exp2/Dmel_fragile_complexes.bed")

require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
require(ChIPseeker)
require(org.Dm.eg.db)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

peak <- fragile_complexes
seqlevelsStyle(peak) = "UCSC"
peakAnno <- annotatePeak(peak,
                         tssRegion=c(-3000,3000),
                         TxDb = txdb,
                         annoDb="org.Dm.eg.db")

pdf("S2_exp2/Pie_chart_fragile_complexes.pdf",width=7,height=3,paper='special') 
plotAnnoPie(peakAnno)
dev.off()

