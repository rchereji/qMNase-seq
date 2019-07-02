library(rtracklayer)
library(GenomeInfoDb)

# Function to align GRanges
Align_GRanges <- function(Profile, ReferenceGRanges)
{
  # Obtain Views for all GRanges that we wish to align
  myViews <- Views(Profile, ReferenceGRanges)
  
  # Convert the RleViewsList (myViews) into a matrix
  AlignedProfilesList <- lapply(myViews, function(gr) t(viewApply(gr, as.vector)))
  AlignedProfiles <- do.call("rbind", AlignedProfilesList)
  
  return(AlignedProfiles)
}

###########
# S2_exp1 #
###########
# Import fragile complexes
fragile_complexes <- import("S2_exp1/Dmel_fragile_complexes.bed")

# Get 4kb windows centered on the fragile complexes
ref_GRs <- resize(fragile_complexes, width=4001, fix="center")

# Create GRanges that contain the entire chromosomes
chrLen <- c(23513712,25286936,28110227,32079331,1348131,23542271,3667352)
names(chrLen) <- c("2L","2R","3L","3R","4","X","Y")
whole_chrom <- GRanges(seqnames = names(chrLen),
                       ranges = IRanges(start = 1, width = chrLen),
                       strand = rep("+",length(chrLen)),
                       seqlengths = chrLen)

# Eliminate the windows that extend outside of the chromosome limits
ref_GRs <- subsetByOverlaps(ref_GRs, whole_chrom, type="within")
seqlengths(ref_GRs) <- seqlengths(whole_chrom)

# Import the ChIP Atlas annotations of TF binding sites
gr_TFs <- GRanges(import("ChIPAtlas_S2cells_TFs.bed"))
seqlevelsStyle(gr_TFs) <- "NCBI"
gr_TFs <- subsetByOverlaps(gr_TFs, whole_chrom, type="within")
seqlevels(gr_TFs) <- names(chrLen)
seqlengths(gr_TFs) <- seqlengths(whole_chrom)

# Compute the number of annotated TF binding sites at each position
cov_TFs <- coverage(gr_TFs)

# TFs
matrix_TFs <- Align_GRanges(cov_TFs, ref_GRs)
total_TFs <- colSums(matrix_TFs)

x <- seq(-2000, 2000)

library(ggplot2)
df <- data.frame(Position = x/1000, Number_of_annotated_TFBS = Avg_TFs)
p <- ggplot(df, aes(x = Position, y = Number_of_annotated_TFBS)) + 
  geom_line(colour="#56B4E9") +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0), 
                     breaks = seq(-2, 2, 0.5)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_classic() +
  xlab("Position relative to the fragile complex (kb)") + 
  ylab("Number of annotated TF binding sites")
p

ggsave(filename = "S2_exp1/ChIP_atlas_TFBS_S2_cells.pdf", 
       plot = p, width = 4, height = 3, units = "in")



###########
# S2_exp1 #
###########
# Import fragile complexes
fragile_complexes <- import("S2_exp2/Dmel_fragile_complexes.bed")

# Get 4kb windows centered on the fragile complexes
ref_GRs <- resize(fragile_complexes, width=4001, fix="center")

# Create GRanges that contain the entire chromosomes
chrLen <- c(23513712,25286936,28110227,32079331,1348131,23542271,3667352)
names(chrLen) <- c("2L","2R","3L","3R","4","X","Y")
whole_chrom <- GRanges(seqnames = names(chrLen),
                       ranges = IRanges(start = 1, width = chrLen),
                       strand = rep("+",length(chrLen)),
                       seqlengths = chrLen)

# Eliminate the windows that extend outside of the chromosome limits
ref_GRs <- subsetByOverlaps(ref_GRs, whole_chrom, type="within")
seqlengths(ref_GRs) <- seqlengths(whole_chrom)

# Import the ChIP Atlas annotations of TF binding sites
gr_TFs <- GRanges(import("ChIPAtlas_S2cells_TFs.bed"))
seqlevelsStyle(gr_TFs) <- "NCBI"
gr_TFs <- subsetByOverlaps(gr_TFs, whole_chrom, type="within")
seqlevels(gr_TFs) <- names(chrLen)
seqlengths(gr_TFs) <- seqlengths(whole_chrom)

# Compute the number of annotated TF binding sites at each position
cov_TFs <- coverage(gr_TFs)

# TFs
matrix_TFs <- Align_GRanges(cov_TFs, ref_GRs)
total_TFs <- colSums(matrix_TFs)

x <- seq(-2000, 2000)

library(ggplot2)
df <- data.frame(Position = x/1000, Number_of_annotated_TFBS = Avg_TFs)
p <- ggplot(df, aes(x = Position, y = Number_of_annotated_TFBS)) + 
  geom_line(colour="#56B4E9") +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0), 
                     breaks = seq(-2, 2, 0.5)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_classic() +
  xlab("Position relative to the fragile complex (kb)") + 
  ylab("Number of annotated TF binding sites")
p

ggsave(filename = "S2_exp2/ChIP_atlas_TFBS_S2_cells.pdf", 
       plot = p, width = 4, height = 3, units = "in")