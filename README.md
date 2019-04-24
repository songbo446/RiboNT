# RiboNT
RiboNT, a noise-tolerant predictor of ORF from RPFs

####################################################################################

RiboNT: a noise-tolerant ORF predictor

Bo Song (songbo446@yeah.net)

Shenzhen University

####################################################################################

0.	DEPENDENCIES

    0.1 R packages:

        multitaper
        pheatmap

    0.2 Python packages:

        math
        pickle
        scipy
        multiprocessing
        functools
    
    0.3 Other tools:

        bedtools
        samtools
    
1.	INTRODUCTION

    RiboNT is coded in Python3 and is developed for Linux users. Make sure Python3 has been installed.

2.	INSTALLATION

    Git clone or download zip files and unzip:
    
        cd RiboNT
        chmod +x RiboNT

    Now you can excute directly by typing "YOURPATH/RiboNT" in command line Alternatively, you can also run the program by typing "python3 YOURPATH/RiboNT".

3.	USAGE

    RiboNT requires the inputs of the sequence of reference genome (in fasta format), the genome annotation file (in gtf format) and the alignment of RPFs (in sorted bam)

    There are eight options

    --start  RiboNT use AUG as the start codon by default. As many non-canonical start codons have been identified in many studies. 
    The users can select other triples as start codons. Codons should be seperated by comma if more than one is used. For example, 
    AUG,UUG,CUG,GUG

    --pcov      The minimum coverage of RPF at P-sites for a canidate ORF. This value should be set between 0-1. By default, RiboNT 
    skip if there is no RPF-supported P-site in a candiate ORF.

    --nCores  The number of multiprocesssors. The identification of ORFs could be speed up by using more processors simultaneously. 
    Five tasks will be processed parellely by default. Please note that there is a trade-off between the speed and RAM cost. It's 
    hihgly recommanded that users should test the required RAM for a single task by setting nCores at 1 and run the program to the 
    step of "Predicting ORFs", then decide how many processors should be used according to the capability of their machines.

    --outdir  The output directory of results.

    --prefix  The prefix of output files.
    
    --Rscript   The path to the excutable file of Rscript. RiboNT uses "Rscript" under /usr/bin/ by default. If "Rscript" is installed in your local directory, or you want to use your own path of Rscript, the path should be provided.
    
    --bedtools  The path to the excutable file of bedtools. RiboNT uses "bedtools" under /usr/bin/ by default. If "bedtools" is installed in your local directory, the path should be provided.
    
    --samtools  The path to the excutable file of samtools. RiboNT uses "samtools" under /usr/bin/ by default. If "samtools" is installed in your local directory, the path should be provided.

4.	RESULTS

    RiboNT outputs 8 files and a directory "Plot" in which metagene and multitaper plots are included.

    4.1 prefix_OffSet.txt
    
      It lists the offsets for RPFs in different sizes 
        Column 1: Size of RPF 
        Column 2: Offset 
        Column 3: Proportion

    4.2 prefix_multitaper.Freq.txt

      The results of multitaper tests of RPFs in various sizes Column 1: Size of RPF Column 2: Minus Log10 of P-values Column 3: The peaked frequency

    4.3 DataSumForPrediction.txt

      It includes the RPFs and their amounts that were selected for ORF prediction, as well as the entropy used in this prediction.

    4.4 prefix_summary.txt

      In this file, the numbers of all the types of predicted ORFs are listed.

    4.5 prefix_orf.fa

      The DNA sequences of predicted ORFs in fasta format. The ID line is formmated as: '>'ID_of_ORF Chromosome Strand Name_of_transcript Combined_P-value P1:P2:P3:P4 The ID_of_ORF is consisted of the positional information as formmated as: Name_of_transcript:Chromsome:Strand:Coordinate_of_start:Coordinate_of_stop

    4.6 prefix_orf_aa.fa

      The amino acid sequences translated from prefix_orf.fa

    4.7 prefix_ORF.gff

      The annotation of predicted ORFs in .gff format

    4.8 prefix_orf_table.txt

      The ID, class, DNA and AA sequences of all the predicted ORFs

    4.9 Plot A directory including metagene and multitaper plots.

    
    NOTE:

    Classification of predicted ORFs:

      (i) annotated ORF, the ORFs identical with the annotated ones

      (ii) truncated ORF, the ORFs with the same start or stop codon but shorter than the annotated ones

      (iii) extended ORF, the ORFs with the same start or stop codon but longer than the annotated ones

      (iv) uORF, upstream ORF, ORFs located in 5'UTR

      (v) ouORF, overlapped uORF, ORFs located in 5'UTR and overlapped with the annotated start codon

      (vi) dORF, downstream ORF, ORFs located in 3'UTR

      (vii) odORF, overlapped dORF, ORFs located in 3'UTR and overlapped with the annotated stop codon

      (ix) ncsORF, ORFs located in non-coding RNAs, ORFs predicted from a gene without any annotated CDS will also be classified into ncsORF

      (x) teORF, ORFs located in transponsable elements

      (xi) pORF, ORFs on pseudo genes.

