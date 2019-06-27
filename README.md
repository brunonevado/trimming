
Scripts for trimming alignments and trees before running analysis of selection, as used in Nevado et al. (under review) "Adaptive evolution is common in rapid evolutionary radiations".  

Input should be sequence (fasta) and tree (newick) files for each gene of interest, and it assumes the same set of species are sampled across genes (but not all species have to be present in all genes). Use at your own risk and carefully check your results, and if useful please cite the paper.


Scripts included:

      (1) trimTrees.R

This R script will read a list of fasta files and associated tree files, and remove tips that are too long*. After that, it will cut the tree on internal edges that are too long* and return the subtree that has most tips. It returns new fasta and tree files for each gene, with the subset of species surviving the trimming. 

* Too long means above a defined quantile of the distribution observed across all genes, with values standardised by total tree length and, for tips, calculated for each species. Species names should be the same across trees, or  identical in the first 20 characters.

usage: Rscript --vanilla trimTrees.R treefile fastafile thr suffix minspecies plotN

  treefile - list of tree files in newick format, one per line
  fastafile - list of sequence files in fasta format, one per line, same order as treefile
  thr - quantile to use to trim tips and edges (use e.g. 99 for 99%)
  suffix - will add this suffix to name of output sequence and tree files
  minspecies - will not output result data if fewer than this number of species survive trimming
  plotN - will make plotN plots for genes before and after trimming instead of running the trimming pipeline. For exploratory use only. Set to 0 to run the actual trimming pipeline.

      (2) trimAlign.R
  
This R script will trim beginning and ends of sequences that show high divergence. For each gene, it will calculate the average genetic distance between each sequence and all others along windows of defined size, and cut the beginning/end windows for each sequence that show average divergence above a defined threshold. 

usage: Rscript --vanilla trimAlign.R fastafile thr winlen ta

  fastafile - list of sequence files in fasta format
  thr - threshold used to exclude window due to having too high divergence. Use e.g. 20 to exclude windows that have 20x higher divergence than the average calculated along the entire gene.
  winlen - length (in base pairs) of windows to use
  suffix - will add this suffix to name of output files


      (3) maskshortseqs.pl
      
This perl script will take a list of fasta files, and check each sequence in each file for short stretches of continuous nucleotides flanked by gaps on both sides, and mask these regions (replace with Ns).


usage: perl maskshortseqs.pl -filelist=fastafile -suffix=.ms.fas -minlen=N -mindist=N

  fastafile - list of sequence files in fasta format, one per line
  suffix - will add this suffix to name of output files
  minLen - maximum length of stretches to remove – stretches equal or longer than this will not be masked
  mindist - number of gaps on each side – masking will occur if more than this number of continuous gaps exist on each side

