This is for test. Please run the piple using commands of:

../RiboNT --genome Yeast_chr.fasta --gtf Yeast.gtf --bam YeastRPF.bam --outdir TEST

The results will be output to a fold of "TEST". You can further compare the outputs
in "TEST" and those in "Results".

Note:

It may report some error of "ERROR: Database file ribont.start.bed contains chromosome
chrIII, but the query file does not. Please re-reun with the -g option for a genome file.
See documentation for details". This is because the dataset is too small to include
the alignments to all the chromosomes. SIMPLY IGNORE THIS ERROR and continue.
