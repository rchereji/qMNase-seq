#!/bin/bash

module load ucsc

wigToBigWig data/Spikein_normalized_nucleosome_counts.S2_1min.wig data/dm6.chrom.sizes data/Spikein_normalized_nucleosome_counts.S2_1min.bw
wigToBigWig data/Spikein_normalized_nucleosome_counts.S2_2min.wig data/dm6.chrom.sizes data/Spikein_normalized_nucleosome_counts.S2_2min.bw
wigToBigWig data/Spikein_normalized_nucleosome_counts.S2_5min.wig data/dm6.chrom.sizes data/Spikein_normalized_nucleosome_counts.S2_5min.bw
wigToBigWig data/Spikein_normalized_nucleosome_counts.S2_15min.wig data/dm6.chrom.sizes data/Spikein_normalized_nucleosome_counts.S2_15min.bw
wigToBigWig data/Spikein_normalized_nucleosome_counts.S2_40min.wig data/dm6.chrom.sizes data/Spikein_normalized_nucleosome_counts.S2_40min.bw
wigToBigWig data/Spikein_normalized_nucleosome_counts.S2_60min.wig data/dm6.chrom.sizes data/Spikein_normalized_nucleosome_counts.S2_60min.bw
