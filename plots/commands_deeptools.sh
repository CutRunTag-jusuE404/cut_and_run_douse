computeMatrix scale-regions -S ../../../TASOR_CutAndRun/bigwigs/mutT-H3K9me3-rep1-norm.bw ../../../TASOR_CutAndRun/bigwigs/WTT-H3K9me3-rep1-norm.bw ../../../TASOR_CutAndRun/bigwigs/mutT-H3K9me3-rep2-norm.bw ../../../TASOR_CutAndRun/bigwigs/WTT-H3K9me3-rep2-norm.bw -R ../output/lostPeaksInTKOvsD3.bed -b 1000 -a 1000 --skipZeros --samplesLabel mutT1 WTT1  mutT2WTT2 -o WT_mut_Chip_over_TKOvsD3.gz --outFileSortedRegions WT_mut_Chip_over_TKOvsD3.bed

plotHeatmap --matrixFile WT_mut_Chip_over_TKOvsD3.gz --outFileName out
