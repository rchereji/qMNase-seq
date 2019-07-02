library(GenomicRanges)

###################
# Load data files #
###################
# Function to read GRanges from CSV files
convertDFtoGR <- function(df) {
  gr <- makeGRangesFromDataFrame(df,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=TRUE,
                                 seqnames.field="ChrLabel",
                                 start.field="StartBP",
                                 end.field="EndBP")
  return(gr)
}


# S2 cells
nuc_S1 <- convertDFtoGR(read.csv(file="../../data/S2_exp1/Fit_results.100-200.NF40.S2_exp1.csv", header=TRUE, sep=","))
nuc_S2 <- convertDFtoGR(read.csv(file="../../data/S2_exp2/Fit_results.100-200.NF35.S2_exp2.csv", header=TRUE, sep=","))

# Kc167 cells
nuc_K1 <- convertDFtoGR(read.csv(file="../../data/Kc167_exp1/Fit_results.100-200.NF35.Kc167_exp1.csv", header=TRUE, sep=","))

# Number of nuc. calls
no_S1 <- length(nuc_S1)
no_S1
# no_S1 = 664095 nucleosomes

no_S2 <- length(nuc_S2)
no_S2
# no_S2 = 660838 nucleosomes

no_K1 <- length(nuc_K1)   
no_K1
# no_K1 = 670976 nucleosomes

#########################################
# Check the reproducibility of the fits #
#########################################
# Check the common nucleosome calls (that are shifted by at most 10 bp)
OL_S1_S2 <- findOverlaps(nuc_S1, nuc_S2, minoverlap=127)
length(OL_S1_S2)
# [1] 578314
length(OL_S1_S2) / no_S1
# [1] 0.8708302
length(OL_S1_S2) / no_S2
# [1] 0.8751222
# About 87% of the nucleosome calls are found in both replicates 
# in S2 cells, shifted by at most 10 bp


# Check the common nucleosome calls (that are shifted by at most 20 bp)
OL_S1_S2 <- findOverlaps(nuc_S1, nuc_S2, minoverlap=107)
length(OL_S1_S2)
# [1] 615261
length(OL_S1_S2) / no_S1
# [1] 0.9264653
length(OL_S1_S2) / no_S2
# [1] 0.9310315
# About 93% of the nucleosome calls are found in both replicates 
# in S2 cells, shifted by at most 20 bp


# Check the fits
fits_S1 <- nuc_S1[queryHits(OL_S1_S2)]
fits_S2 <- nuc_S2[subjectHits(OL_S1_S2)]

cor(fits_S1$O_vector, fits_S2$O_vector, method="pearson") 
# [1] 0.7869686

cor(fits_S1$k1_vector, fits_S2$k1_vector, method="pearson") 
# [1] 0.6642765

cor(fits_S1$k2_vector, fits_S2$k2_vector, method="pearson") 
# [1] 0.8410524



#######################################################
# Nucleosomes found in both cell types (S2 and Kc167) #
#######################################################
# Check the common nucleosome calls (that are shifted by at most 20 bp)
OL_S1_K1 <- findOverlaps(nuc_S1, nuc_K1, minoverlap=107)
fits_S1_ <- nuc_S1[queryHits(OL_S1_K1)]
fits_K1_ <- nuc_K1[subjectHits(OL_S1_K1)]
cor(fits_S1_$O_vector, fits_K1_$O_vector, method="pearson") 
# [1] 0.3617206

OL_S2_K1 <- findOverlaps(nuc_S2, nuc_K1, minoverlap=107)
fits_S2_ <- nuc_S2[queryHits(OL_S2_K1)]
fits_K1_ <- nuc_K1[subjectHits(OL_S2_K1)]
cor(fits_S2_$O_vector, fits_K1_$O_vector, method="pearson") 
# [1] 0.3022945

###################
# Plot histograms #
###################
#########################
# Plot histograms of k1 #
#########################
x <- nuc_S1$k1_vector
hist.data <- hist(x, plot=FALSE, breaks = 1001)
hist.data$counts <- log10(hist.data$counts+1)
plot(hist.data, 
     main="S2 cells, replicate 1",
     xlab="Release rate constant (k1)",
     ylab='Number of nucleosomes (log10 scale)')


x <- nuc_S2$k1_vector
hist.data <- hist(x, plot=FALSE, breaks = 1001)
hist.data$counts <- log10(hist.data$counts+1)
plot(hist.data, 
     main="S2 cells, replicate 2",
     xlab="Release rate constant (k1)",
     ylab='Number of nucleosomes (log10 scale)')


x <- nuc_K1$k1_vector
hist.data <- hist(x, plot=FALSE, breaks = 1001)
hist.data$counts <- log10(hist.data$counts+1)
plot(hist.data, 
     main="Kc167 cells, replicate 1",
     xlab="Release rate constant (k1)",
     ylab='Number of nucleosomes (log10 scale)')



##########################################
# Check the fragile complexes (k_1 > 10) #
##########################################
library(GenomeInfoDb)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(ChIPseeker)
library(org.Dm.eg.db)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene


# Annotate fragile complexes
gr <- nuc_S1[nuc_S1$k1_vector > 10]
seqlevelsStyle(gr) = "UCSC"
anno <- annotatePeak(gr,
                     tssRegion=c(-3000,3000),
                     TxDb = txdb,
                     annoDb="org.Dm.eg.db")

pdf("Pie_chart_fragile_complexes_S2_exp1.pdf", width=7, height=3, paper='special') 
plotAnnoPie(anno)
dev.off()

gr <- nuc_S2[nuc_S2$k1_vector > 10]
seqlevelsStyle(gr) = "UCSC"
anno <- annotatePeak(gr,
                     tssRegion=c(-3000,3000),
                     TxDb = txdb,
                     annoDb="org.Dm.eg.db")

pdf("Pie_chart_fragile_complexes_S2_exp2.pdf", width=7, height=3, paper='special') 
plotAnnoPie(anno)
dev.off()

gr <- nuc_K1[nuc_K1$k1_vector > 10]
seqlevelsStyle(gr) = "UCSC"
anno <- annotatePeak(gr,
                     tssRegion=c(-3000,3000),
                     TxDb = txdb,
                     annoDb="org.Dm.eg.db")

pdf("Pie_chart_fragile_complexes_Kc167_exp1.pdf", width=7, height=3, paper='special') 
plotAnnoPie(anno)
dev.off()
