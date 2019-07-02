library(rtracklayer)

###################
# Load nuc. calls #
###################
# S2 cells
nuc_S1 <- import("../../data/S2_exp1/Typical_nucleosome_positions_in_S2_exp1.bed")
nuc_S2 <- import("../../data/S2_exp2/Typical_nucleosome_positions_in_S2_exp2.bed")

# Kc167 cells
nuc_K1 <- import("../../data/Kc167_exp1/Typical_nucleosome_positions_in_Kc167_exp1.bed")

# Number of nuc. calls
no_S1 <- length(nuc_S1)   # no_S1 = 681925 nucleosomes
no_S2 <- length(nuc_S2)   # no_S2 = 677654 nucleosomes
no_K1 <- length(nuc_K1)   # no_K1 = 687281 nucleosomes

# Compute the total length of DNA covered by these nucleosomes
sum(width(reduce(nuc_S1)))
# [1] 97513123
sum(width(reduce(nuc_S2)))
# [1] 96740885
sum(width(reduce(nuc_K1)))
# [1] 97983761

###########################################
# Check the reproducibility of nuc. calls #
###########################################
# Check the common nucleosome calls (that are shifted by at most 10 bp)
OL_S1_S2 <- findOverlaps(nuc_S1, nuc_S2, minoverlap=127)
length(OL_S1_S2)
# [1] 593112
length(OL_S1_S2) / no_S1
# [1] 0.8697613
length(OL_S1_S2) / no_S2
# [1] 0.8752431
# About 87% of the nucleosome calls are found in both replicates 
# in S2 cells, shifted by at most 10 bp


# Check the common nucleosome calls (that are shifted by at most 20 bp)
OL_S1_S2 <- findOverlaps(nuc_S1, nuc_S2, minoverlap=107)
length(OL_S1_S2)
# [1] 630937
length(OL_S1_S2) / no_S1
# [1] 0.9252293
length(OL_S1_S2) / no_S2
# [1] 0.9310607
# About 93% of the nucleosome calls are found in both replicates 
# in S2 cells, shifted by at most 20 bp


# Check the common nucleosome calls (that are shifted by at most 30 bp)
OL_S1_S2 <- findOverlaps(nuc_S1, nuc_S2, minoverlap=87)
length(OL_S1_S2)
# [1] 646594
length(OL_S1_S2) / no_S1
# [1] 0.9481893
length(OL_S1_S2) / no_S2
# [1] 0.9541654
# About 95% of the nucleosome calls are found in both replicates 
# in S2 cells, shifted by at most 20 bp



#######################################################
# Check the common nuc. calls from S2 and Kc167 cells #
#######################################################
# Check the overlaps between nuc. calls in S2 and Kc167 cells 
# (shifted by at most 20 bp)
OL_S1_K1 <- findOverlaps(nuc_S1, nuc_K1, minoverlap=107)
length(OL_S1_K1)
# [1] 548554
length(OL_S1_K1) / no_S1
# [1] 0.8044198
length(OL_S1_K1) / no_K1
# [1] 0.798151
# About 80% of the nucleosome calls are found in both cell types 
# (shifted by at most 20 bp)


# Check the overlaps between nuc. calls in S2 and Kc167 cells 
# (shifted by at most 20 bp)
OL_S2_K1 <- findOverlaps(nuc_S2, nuc_K1, minoverlap=107)
length(OL_S2_K1)
# [1] 536504
length(OL_S2_K1) / no_S2
# [1] 0.7917079
length(OL_S2_K1) / no_K1
# [1] 0.7806181
# About 80% of the nucleosome calls are found in both cell types 
# (shifted by at most 20 bp)


# Check the overlaps between nuc. calls in S2 and Kc167 cells 
# (shifted by at most 30 bp)
OL_S1_K1 <- findOverlaps(nuc_S1, nuc_K1, minoverlap=87)
length(OL_S1_K1)
# [1] 598964
length(OL_S1_K1) / no_S1
# [1] 0.8783429
length(OL_S1_K1) / no_K1
# [1] 0.8783429
# About 87% of the nucleosome calls are found in both cell types 
# (shifted by at most 30 bp)



################################################################
# Check the nucleosomes that were called only in one cell type #
################################################################
library(GenomeInfoDb)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(ChIPseeker)
library(org.Dm.eg.db)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene


# Annotate nucleosomes that were found in both sets (shifted by at most 20 bp)
OL_S1_K1 <- findOverlaps(nuc_S1, nuc_K1, minoverlap=107)
common_nucs <- nuc_S1[queryHits(OL_S1_K1)]
gr <- common_nucs
seqlevelsStyle(gr) = "UCSC"
anno <- annotatePeak(gr,
                     tssRegion=c(-3000,3000),
                     TxDb = txdb,
                     annoDb="org.Dm.eg.db")

pdf("Pie_chart_common_nucleosomes.pdf", width=7, height=3, paper='special') 
plotAnnoPie(anno)
dev.off()

# Annotate nucleosomes that are unique or remodeled (shifted by > 20 bp)
OL_S1_K1 <- findOverlaps(nuc_S1, nuc_K1, minoverlap=107)
unique_nuc_S1 <- nuc_S1[-queryHits(OL_S1_K1)]
unique_nuc_K1 <- nuc_S1[-subjectHits(OL_S1_K1)]
unique_nucs <- c(unique_nuc_S1, unique_nuc_K1)
gr <- unique_nucs
seqlevelsStyle(gr) = "UCSC"
anno <- annotatePeak(gr,
                     tssRegion=c(-3000,3000),
                     TxDb = txdb,
                     annoDb="org.Dm.eg.db")

pdf("Pie_chart_unique_nucleosomes.pdf", width=7, height=3, paper='special') 
plotAnnoPie(anno)
dev.off()
