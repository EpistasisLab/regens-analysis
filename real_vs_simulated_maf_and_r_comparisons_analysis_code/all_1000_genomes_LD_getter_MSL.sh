#!/bin/bash
#BSUB -J all_1000_genomes_LD_getter
#BSUB -o all_1000_genomes_LD_getter.out
#BSUB -e all_1000_genomes_LD_getter.error
#BSUB -R "rusage[mem=50000MB]"
#BSUB -M 50000MB
source activate regens
module load plink/1.90Beta6.18

plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 1 --make-bed --out all_1000_genomes_MSL_processed_chr1
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 2 --make-bed --out all_1000_genomes_MSL_processed_chr2
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 3 --make-bed --out all_1000_genomes_MSL_processed_chr3
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 4 --make-bed --out all_1000_genomes_MSL_processed_chr4
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 5 --make-bed --out all_1000_genomes_MSL_processed_chr5
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 6 --make-bed --out all_1000_genomes_MSL_processed_chr6
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 7 --make-bed --out all_1000_genomes_MSL_processed_chr7
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 8 --make-bed --out all_1000_genomes_MSL_processed_chr8
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 9 --make-bed --out all_1000_genomes_MSL_processed_chr9
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 10 --make-bed --out all_1000_genomes_MSL_processed_chr10
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 11 --make-bed --out all_1000_genomes_MSL_processed_chr11
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 12 --make-bed --out all_1000_genomes_MSL_processed_chr12
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 13 --make-bed --out all_1000_genomes_MSL_processed_chr13
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 14 --make-bed --out all_1000_genomes_MSL_processed_chr14
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 15 --make-bed --out all_1000_genomes_MSL_processed_chr15
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 16 --make-bed --out all_1000_genomes_MSL_processed_chr16
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 17 --make-bed --out all_1000_genomes_MSL_processed_chr17
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 18 --make-bed --out all_1000_genomes_MSL_processed_chr18
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 19 --make-bed --out all_1000_genomes_MSL_processed_chr19
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 20 --make-bed --out all_1000_genomes_MSL_processed_chr20
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 21 --make-bed --out all_1000_genomes_MSL_processed_chr21
plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --chr 22 --make-bed --out all_1000_genomes_MSL_processed_chr22

plink --bfile all_1000_genomes_MSL_processed_chr1 --merge-list all_1000_genomes_MSL_processed_merge.txt --make-bed --out all_1000_genomes_MSL_processed

plink --bfile all_1000_genomes_MSL_processed --chr 1 --make-bed --out all_1000_genomes_MSL_processed_chr1
plink --bfile all_1000_genomes_MSL_processed --chr 2 --make-bed --out all_1000_genomes_MSL_processed_chr2
plink --bfile all_1000_genomes_MSL_processed --chr 3 --make-bed --out all_1000_genomes_MSL_processed_chr3
plink --bfile all_1000_genomes_MSL_processed --chr 4 --make-bed --out all_1000_genomes_MSL_processed_chr4
plink --bfile all_1000_genomes_MSL_processed --chr 5 --make-bed --out all_1000_genomes_MSL_processed_chr5
plink --bfile all_1000_genomes_MSL_processed --chr 6 --make-bed --out all_1000_genomes_MSL_processed_chr6
plink --bfile all_1000_genomes_MSL_processed --chr 7 --make-bed --out all_1000_genomes_MSL_processed_chr7
plink --bfile all_1000_genomes_MSL_processed --chr 8 --make-bed --out all_1000_genomes_MSL_processed_chr8
plink --bfile all_1000_genomes_MSL_processed --chr 9 --make-bed --out all_1000_genomes_MSL_processed_chr9
plink --bfile all_1000_genomes_MSL_processed --chr 10 --make-bed --out all_1000_genomes_MSL_processed_chr10
plink --bfile all_1000_genomes_MSL_processed --chr 11 --make-bed --out all_1000_genomes_MSL_processed_chr11
plink --bfile all_1000_genomes_MSL_processed --chr 12 --make-bed --out all_1000_genomes_MSL_processed_chr12
plink --bfile all_1000_genomes_MSL_processed --chr 13 --make-bed --out all_1000_genomes_MSL_processed_chr13
plink --bfile all_1000_genomes_MSL_processed --chr 14 --make-bed --out all_1000_genomes_MSL_processed_chr14
plink --bfile all_1000_genomes_MSL_processed --chr 15 --make-bed --out all_1000_genomes_MSL_processed_chr15
plink --bfile all_1000_genomes_MSL_processed --chr 16 --make-bed --out all_1000_genomes_MSL_processed_chr16
plink --bfile all_1000_genomes_MSL_processed --chr 17 --make-bed --out all_1000_genomes_MSL_processed_chr17
plink --bfile all_1000_genomes_MSL_processed --chr 18 --make-bed --out all_1000_genomes_MSL_processed_chr18
plink --bfile all_1000_genomes_MSL_processed --chr 19 --make-bed --out all_1000_genomes_MSL_processed_chr19
plink --bfile all_1000_genomes_MSL_processed --chr 20 --make-bed --out all_1000_genomes_MSL_processed_chr20
plink --bfile all_1000_genomes_MSL_processed --chr 21 --make-bed --out all_1000_genomes_MSL_processed_chr21
plink --bfile all_1000_genomes_MSL_processed --chr 22 --make-bed --out all_1000_genomes_MSL_processed_chr22

python -m regens --in ../real_1000_genomes_input_for_analysis/all_1000_genomes_MSL_processed --out all_1000_genomes_MSL_simulated --simulate_nbreakpoints 4 --simulate_nsamples 100000 --population_code MSL --human_genome_version hg19

plink --bfile all_1000_genomes_MSL_processed_chr1 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr1
plink --bfile all_1000_genomes_MSL_processed_chr2 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr2
plink --bfile all_1000_genomes_MSL_processed_chr3 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr3
plink --bfile all_1000_genomes_MSL_processed_chr4 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr4
plink --bfile all_1000_genomes_MSL_processed_chr5 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr5
plink --bfile all_1000_genomes_MSL_processed_chr6 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr6
plink --bfile all_1000_genomes_MSL_processed_chr7 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr7
plink --bfile all_1000_genomes_MSL_processed_chr8 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr8
plink --bfile all_1000_genomes_MSL_processed_chr9 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr9
plink --bfile all_1000_genomes_MSL_processed_chr10 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr10
plink --bfile all_1000_genomes_MSL_processed_chr11 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr11
plink --bfile all_1000_genomes_MSL_processed_chr12 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr12
plink --bfile all_1000_genomes_MSL_processed_chr13 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr13
plink --bfile all_1000_genomes_MSL_processed_chr14 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr14
plink --bfile all_1000_genomes_MSL_processed_chr15 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr15
plink --bfile all_1000_genomes_MSL_processed_chr16 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr16
plink --bfile all_1000_genomes_MSL_processed_chr17 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr17
plink --bfile all_1000_genomes_MSL_processed_chr18 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr18
plink --bfile all_1000_genomes_MSL_processed_chr19 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr19
plink --bfile all_1000_genomes_MSL_processed_chr20 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr20
plink --bfile all_1000_genomes_MSL_processed_chr21 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr21
plink --bfile all_1000_genomes_MSL_processed_chr22 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_processed_chr22

plink --bfile all_1000_genomes_MSL_simulated_chr1 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr1
plink --bfile all_1000_genomes_MSL_simulated_chr2 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr2
plink --bfile all_1000_genomes_MSL_simulated_chr3 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr3
plink --bfile all_1000_genomes_MSL_simulated_chr4 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr4
plink --bfile all_1000_genomes_MSL_simulated_chr5 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr5
plink --bfile all_1000_genomes_MSL_simulated_chr6 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr6
plink --bfile all_1000_genomes_MSL_simulated_chr7 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr7
plink --bfile all_1000_genomes_MSL_simulated_chr8 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr8
plink --bfile all_1000_genomes_MSL_simulated_chr9 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr9
plink --bfile all_1000_genomes_MSL_simulated_chr10 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr10
plink --bfile all_1000_genomes_MSL_simulated_chr11 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr11
plink --bfile all_1000_genomes_MSL_simulated_chr12 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr12
plink --bfile all_1000_genomes_MSL_simulated_chr13 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr13
plink --bfile all_1000_genomes_MSL_simulated_chr14 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr14
plink --bfile all_1000_genomes_MSL_simulated_chr15 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr15
plink --bfile all_1000_genomes_MSL_simulated_chr16 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr16
plink --bfile all_1000_genomes_MSL_simulated_chr17 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr17
plink --bfile all_1000_genomes_MSL_simulated_chr18 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr18
plink --bfile all_1000_genomes_MSL_simulated_chr19 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr19
plink --bfile all_1000_genomes_MSL_simulated_chr20 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr20
plink --bfile all_1000_genomes_MSL_simulated_chr21 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr21
plink --bfile all_1000_genomes_MSL_simulated_chr22 --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out all_1000_genomes_MSL_simulated_chr22

python REALGenomeSIM_LD_getter.py --ref all_1000_genomes_MSL_processed --sim all_1000_genomes_MSL_simulated

rm all_1000_genomes_MSL_processed.bim
rm all_1000_genomes_MSL_processed.fam
rm all_1000_genomes_MSL_processed.bed
rm all_1000_genomes_MSL_processed.log
rm all_1000_genomes_MSL_processed_chr1.bed
rm all_1000_genomes_MSL_simulated_chr1.bed
rm all_1000_genomes_MSL_processed_chr1.bim
rm all_1000_genomes_MSL_simulated_chr1.bim
rm all_1000_genomes_MSL_processed_chr1.fam
rm all_1000_genomes_MSL_simulated_chr1.fam
rm all_1000_genomes_MSL_processed_chr1.log
rm all_1000_genomes_MSL_simulated_chr1.log
rm all_1000_genomes_MSL_processed_chr1.ld
rm all_1000_genomes_MSL_simulated_chr1.ld
rm all_1000_genomes_MSL_processed_chr2.bed
rm all_1000_genomes_MSL_simulated_chr2.bed
rm all_1000_genomes_MSL_processed_chr2.bim
rm all_1000_genomes_MSL_simulated_chr2.bim
rm all_1000_genomes_MSL_processed_chr2.fam
rm all_1000_genomes_MSL_simulated_chr2.fam
rm all_1000_genomes_MSL_processed_chr2.log
rm all_1000_genomes_MSL_simulated_chr2.log
rm all_1000_genomes_MSL_processed_chr2.ld
rm all_1000_genomes_MSL_simulated_chr2.ld
rm all_1000_genomes_MSL_processed_chr3.bed
rm all_1000_genomes_MSL_simulated_chr3.bed
rm all_1000_genomes_MSL_processed_chr3.bim
rm all_1000_genomes_MSL_simulated_chr3.bim
rm all_1000_genomes_MSL_processed_chr3.fam
rm all_1000_genomes_MSL_simulated_chr3.fam
rm all_1000_genomes_MSL_processed_chr3.log
rm all_1000_genomes_MSL_simulated_chr3.log
rm all_1000_genomes_MSL_processed_chr3.ld
rm all_1000_genomes_MSL_simulated_chr3.ld
rm all_1000_genomes_MSL_processed_chr4.bed
rm all_1000_genomes_MSL_simulated_chr4.bed
rm all_1000_genomes_MSL_processed_chr4.bim
rm all_1000_genomes_MSL_simulated_chr4.bim
rm all_1000_genomes_MSL_processed_chr4.fam
rm all_1000_genomes_MSL_simulated_chr4.fam
rm all_1000_genomes_MSL_processed_chr4.log
rm all_1000_genomes_MSL_simulated_chr4.log
rm all_1000_genomes_MSL_processed_chr4.ld
rm all_1000_genomes_MSL_simulated_chr4.ld
rm all_1000_genomes_MSL_processed_chr5.bed
rm all_1000_genomes_MSL_simulated_chr5.bed
rm all_1000_genomes_MSL_processed_chr5.bim
rm all_1000_genomes_MSL_simulated_chr5.bim
rm all_1000_genomes_MSL_processed_chr5.fam
rm all_1000_genomes_MSL_simulated_chr5.fam
rm all_1000_genomes_MSL_processed_chr5.log
rm all_1000_genomes_MSL_simulated_chr5.log
rm all_1000_genomes_MSL_processed_chr5.ld
rm all_1000_genomes_MSL_simulated_chr5.ld
rm all_1000_genomes_MSL_processed_chr6.bed
rm all_1000_genomes_MSL_simulated_chr6.bed
rm all_1000_genomes_MSL_processed_chr6.bim
rm all_1000_genomes_MSL_simulated_chr6.bim
rm all_1000_genomes_MSL_processed_chr6.fam
rm all_1000_genomes_MSL_simulated_chr6.fam
rm all_1000_genomes_MSL_processed_chr6.log
rm all_1000_genomes_MSL_simulated_chr6.log
rm all_1000_genomes_MSL_processed_chr6.ld
rm all_1000_genomes_MSL_simulated_chr6.ld
rm all_1000_genomes_MSL_processed_chr7.bed
rm all_1000_genomes_MSL_simulated_chr7.bed
rm all_1000_genomes_MSL_processed_chr7.bim
rm all_1000_genomes_MSL_simulated_chr7.bim
rm all_1000_genomes_MSL_processed_chr7.fam
rm all_1000_genomes_MSL_simulated_chr7.fam
rm all_1000_genomes_MSL_processed_chr7.log
rm all_1000_genomes_MSL_simulated_chr7.log
rm all_1000_genomes_MSL_processed_chr7.ld
rm all_1000_genomes_MSL_simulated_chr7.ld
rm all_1000_genomes_MSL_processed_chr8.bed
rm all_1000_genomes_MSL_simulated_chr8.bed
rm all_1000_genomes_MSL_processed_chr8.bim
rm all_1000_genomes_MSL_simulated_chr8.bim
rm all_1000_genomes_MSL_processed_chr8.fam
rm all_1000_genomes_MSL_simulated_chr8.fam
rm all_1000_genomes_MSL_processed_chr8.log
rm all_1000_genomes_MSL_simulated_chr8.log
rm all_1000_genomes_MSL_processed_chr8.ld
rm all_1000_genomes_MSL_simulated_chr8.ld
rm all_1000_genomes_MSL_processed_chr9.bed
rm all_1000_genomes_MSL_simulated_chr9.bed
rm all_1000_genomes_MSL_processed_chr9.bim
rm all_1000_genomes_MSL_simulated_chr9.bim
rm all_1000_genomes_MSL_processed_chr9.fam
rm all_1000_genomes_MSL_simulated_chr9.fam
rm all_1000_genomes_MSL_processed_chr9.log
rm all_1000_genomes_MSL_simulated_chr9.log
rm all_1000_genomes_MSL_processed_chr9.ld
rm all_1000_genomes_MSL_simulated_chr9.ld
rm all_1000_genomes_MSL_processed_chr10.bed
rm all_1000_genomes_MSL_simulated_chr10.bed
rm all_1000_genomes_MSL_processed_chr10.bim
rm all_1000_genomes_MSL_simulated_chr10.bim
rm all_1000_genomes_MSL_processed_chr10.fam
rm all_1000_genomes_MSL_simulated_chr10.fam
rm all_1000_genomes_MSL_processed_chr10.log
rm all_1000_genomes_MSL_simulated_chr10.log
rm all_1000_genomes_MSL_processed_chr10.ld
rm all_1000_genomes_MSL_simulated_chr10.ld
rm all_1000_genomes_MSL_processed_chr11.bed
rm all_1000_genomes_MSL_simulated_chr11.bed
rm all_1000_genomes_MSL_processed_chr11.bim
rm all_1000_genomes_MSL_simulated_chr11.bim
rm all_1000_genomes_MSL_processed_chr11.fam
rm all_1000_genomes_MSL_simulated_chr11.fam
rm all_1000_genomes_MSL_processed_chr11.log
rm all_1000_genomes_MSL_simulated_chr11.log
rm all_1000_genomes_MSL_processed_chr11.ld
rm all_1000_genomes_MSL_simulated_chr11.ld
rm all_1000_genomes_MSL_processed_chr12.bed
rm all_1000_genomes_MSL_simulated_chr12.bed
rm all_1000_genomes_MSL_processed_chr12.bim
rm all_1000_genomes_MSL_simulated_chr12.bim
rm all_1000_genomes_MSL_processed_chr12.fam
rm all_1000_genomes_MSL_simulated_chr12.fam
rm all_1000_genomes_MSL_processed_chr12.log
rm all_1000_genomes_MSL_simulated_chr12.log
rm all_1000_genomes_MSL_processed_chr12.ld
rm all_1000_genomes_MSL_simulated_chr12.ld
rm all_1000_genomes_MSL_processed_chr13.bed
rm all_1000_genomes_MSL_simulated_chr13.bed
rm all_1000_genomes_MSL_processed_chr13.bim
rm all_1000_genomes_MSL_simulated_chr13.bim
rm all_1000_genomes_MSL_processed_chr13.fam
rm all_1000_genomes_MSL_simulated_chr13.fam
rm all_1000_genomes_MSL_processed_chr13.log
rm all_1000_genomes_MSL_simulated_chr13.log
rm all_1000_genomes_MSL_processed_chr13.ld
rm all_1000_genomes_MSL_simulated_chr13.ld
rm all_1000_genomes_MSL_processed_chr14.bed
rm all_1000_genomes_MSL_simulated_chr14.bed
rm all_1000_genomes_MSL_processed_chr14.bim
rm all_1000_genomes_MSL_simulated_chr14.bim
rm all_1000_genomes_MSL_processed_chr14.fam
rm all_1000_genomes_MSL_simulated_chr14.fam
rm all_1000_genomes_MSL_processed_chr14.log
rm all_1000_genomes_MSL_simulated_chr14.log
rm all_1000_genomes_MSL_processed_chr14.ld
rm all_1000_genomes_MSL_simulated_chr14.ld
rm all_1000_genomes_MSL_processed_chr15.bed
rm all_1000_genomes_MSL_simulated_chr15.bed
rm all_1000_genomes_MSL_processed_chr15.bim
rm all_1000_genomes_MSL_simulated_chr15.bim
rm all_1000_genomes_MSL_processed_chr15.fam
rm all_1000_genomes_MSL_simulated_chr15.fam
rm all_1000_genomes_MSL_processed_chr15.log
rm all_1000_genomes_MSL_simulated_chr15.log
rm all_1000_genomes_MSL_processed_chr15.ld
rm all_1000_genomes_MSL_simulated_chr15.ld
rm all_1000_genomes_MSL_processed_chr16.bed
rm all_1000_genomes_MSL_simulated_chr16.bed
rm all_1000_genomes_MSL_processed_chr16.bim
rm all_1000_genomes_MSL_simulated_chr16.bim
rm all_1000_genomes_MSL_processed_chr16.fam
rm all_1000_genomes_MSL_simulated_chr16.fam
rm all_1000_genomes_MSL_processed_chr16.log
rm all_1000_genomes_MSL_simulated_chr16.log
rm all_1000_genomes_MSL_processed_chr16.ld
rm all_1000_genomes_MSL_simulated_chr16.ld
rm all_1000_genomes_MSL_processed_chr17.bed
rm all_1000_genomes_MSL_simulated_chr17.bed
rm all_1000_genomes_MSL_processed_chr17.bim
rm all_1000_genomes_MSL_simulated_chr17.bim
rm all_1000_genomes_MSL_processed_chr17.fam
rm all_1000_genomes_MSL_simulated_chr17.fam
rm all_1000_genomes_MSL_processed_chr17.log
rm all_1000_genomes_MSL_simulated_chr17.log
rm all_1000_genomes_MSL_processed_chr17.ld
rm all_1000_genomes_MSL_simulated_chr17.ld
rm all_1000_genomes_MSL_processed_chr18.bed
rm all_1000_genomes_MSL_simulated_chr18.bed
rm all_1000_genomes_MSL_processed_chr18.bim
rm all_1000_genomes_MSL_simulated_chr18.bim
rm all_1000_genomes_MSL_processed_chr18.fam
rm all_1000_genomes_MSL_simulated_chr18.fam
rm all_1000_genomes_MSL_processed_chr18.log
rm all_1000_genomes_MSL_simulated_chr18.log
rm all_1000_genomes_MSL_processed_chr18.ld
rm all_1000_genomes_MSL_simulated_chr18.ld
rm all_1000_genomes_MSL_processed_chr19.bed
rm all_1000_genomes_MSL_simulated_chr19.bed
rm all_1000_genomes_MSL_processed_chr19.bim
rm all_1000_genomes_MSL_simulated_chr19.bim
rm all_1000_genomes_MSL_processed_chr19.fam
rm all_1000_genomes_MSL_simulated_chr19.fam
rm all_1000_genomes_MSL_processed_chr19.log
rm all_1000_genomes_MSL_simulated_chr19.log
rm all_1000_genomes_MSL_processed_chr19.ld
rm all_1000_genomes_MSL_simulated_chr19.ld
rm all_1000_genomes_MSL_processed_chr20.bed
rm all_1000_genomes_MSL_simulated_chr20.bed
rm all_1000_genomes_MSL_processed_chr20.bim
rm all_1000_genomes_MSL_simulated_chr20.bim
rm all_1000_genomes_MSL_processed_chr20.fam
rm all_1000_genomes_MSL_simulated_chr20.fam
rm all_1000_genomes_MSL_processed_chr20.log
rm all_1000_genomes_MSL_simulated_chr20.log
rm all_1000_genomes_MSL_processed_chr20.ld
rm all_1000_genomes_MSL_simulated_chr20.ld
rm all_1000_genomes_MSL_processed_chr21.bed
rm all_1000_genomes_MSL_simulated_chr21.bed
rm all_1000_genomes_MSL_processed_chr21.bim
rm all_1000_genomes_MSL_simulated_chr21.bim
rm all_1000_genomes_MSL_processed_chr21.fam
rm all_1000_genomes_MSL_simulated_chr21.fam
rm all_1000_genomes_MSL_processed_chr21.log
rm all_1000_genomes_MSL_simulated_chr21.log
rm all_1000_genomes_MSL_processed_chr21.ld
rm all_1000_genomes_MSL_simulated_chr21.ld
rm all_1000_genomes_MSL_processed_chr22.bed
rm all_1000_genomes_MSL_simulated_chr22.bed
rm all_1000_genomes_MSL_processed_chr22.bim
rm all_1000_genomes_MSL_simulated_chr22.bim
rm all_1000_genomes_MSL_processed_chr22.fam
rm all_1000_genomes_MSL_simulated_chr22.fam
rm all_1000_genomes_MSL_processed_chr22.log
rm all_1000_genomes_MSL_simulated_chr22.log
rm all_1000_genomes_MSL_processed_chr22.ld
rm all_1000_genomes_MSL_simulated_chr22.ld
