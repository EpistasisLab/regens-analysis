# About REALGenomeSIM

**REALGenomeSIM** simulates GWAS data with **R**ealistic **E**pistatic effects, **A**dditive effects, and **L**D patterns with **Genome** **S**ubsampling and **I**ntegrated **M**odeling. It divides each simulated genome into chunks such that the probability that any genomic position demarates two chunks equals the probability that this position would demarcate a real (given) recombination event. Each chunk is then filled with the homologous genomic subset of a real whole genome that is sampled uniformly (with replacement) from an input plink bed file. Since REALGenomeSIM's in-silico recombinations occur in the same proportions as equivalent real recombination events, each simulated individual is likely similar to at least one member of the input population whose genome is not included in the data. 

<p align="center">
<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=750/>
</p>

The Triadsim algorithm has used this method to simulate LD patterns that are almost indistinguishable from those of the input Dataset. REALGenomeSIM simulates equally realistic data, and it has the following improvements:

<details>
  <summary>REALGenomeSIM is approximately 50 times faster than Triadsim.</summary>
<p align="center">
<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
</p>
</details>

<details>
  <summary>REALGenomeSIM can use any whole genomes as an input dataset. This allows REALGenomeSIM to simulate new individuals from any publically available GWAS dataset.
    The following figures compare the publically available 1000 genomes project subpopulations to datasets that REALGenomeSIM simulated from each subpopulation.</summary>
  
   <details>
     <summary>These figures compare r values and mafs between every real and simulated 1000 genomes subpopulation for SNP pairs less than  200 kilobases apart.</summary>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_ACB.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_ASW.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_BEB.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_CDX.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_CEU.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_CHB.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_CHS.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_CLM.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_ESN.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_FIN.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_GBR.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_GIH.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_GWD.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_IBS.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_ITU.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_JPT.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_ASW.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_KHV.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_LWK.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_MSL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_MXL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_PEL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_PJL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_PUR.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_STU.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_TSI.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_val_maf_comparison_YRI.png" width=1000/>
   </p>
   </details>
   
   <details>
     <summary>These figures compare profiles of SNP pairs'absolute r values vs distance between every real and simulated 1000 genomes subpopulation for SNP pairs less than  200 kilobases apart.</summary>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_ACB.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_ASW.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_BEB.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_CDX.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_CEU.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_CHB.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_CHS.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_CLM.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_ESN.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_FIN.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_GBR.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_GIH.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_GWD.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_IBS.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_ITU.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_JPT.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_ASW.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_KHV.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_LWK.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_MSL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_MXL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_PEL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_PJL.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_PUR.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_STU.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_TSI.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/real_vs_sim_r_vs_distance_profile_comparison_YRI.png" width=1000/>
   </p>
   </details>
  
   <details>
     <summary>These figures compare TSNE plots of the first 10 principal components for real and simulated 1000 genomes subpopulations.</summary>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/TSNE1_vs_TSNE2_for_1000_genome_African_subpopulations.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/TSNE1_vs_TSNE2_for_1000_genome_European_subpopulations.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/TSNE1_vs_TSNE2_for_1000_genome_American_subpopulations.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/TSNE1_vs_TSNE2_for_1000_genome_East_Asian_subpopulations.png" width=1000/>
   </p>
   <p align="center">
   <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/TSNE1_vs_TSNE2_for_1000_genome_South_Asian_subpopulations.png" width=1000/>
   </p>
   </details>
</details>

<details>
  <summary>REALGenomeSIM uses the output recombination rate maps from () and converts them to probabilities of drawing each simulated breakpoint at a specific genomic location. Depending on your intent, the following optiona are available:</summary>

   *  If you want to simulate arbitrary datasets of realistic human genomes, then this package includes all samples from the 1000 genomes project that I have filtered down to 500000 SNPs. In summary, I kept every biallelic SNP such that every subpopulation contains at least two instances of the minor allele. Exact thinning methods are here. 

   *  If you want to simulate realistic human genomes that are similar to a specific dataset, then you should choose the 1000 genomes project subpopulation that is closest to the population of your sample in the () argument. 

   *  If you want to simulate non-human genomes, then you will need to find external recombination rate values for all genomic regions that you intend to simulate. This information will need to be formatted as follows:
</details>

# Installing REALGenomeSIM 

REALGenomeSIM and its dependencies can be installed with pip as follows: 

```
pip install bed-reader==0.0.3a0
pip install numpy
pip install pandas
pip install matplotlib
pip install scipy
pip install scikit-learn
line to install REALGenomeSIM
```

# How to Simulate Whole Genomes without Simulating Phenotypes 

The following command uses all_1000_genomes_ACB_processed.bed, all_1000_genomes_ACB_processed.bim, and all_1000_genomes_ACB_processed.fam to simulate 10000 individuals

```
python REALGenomeSIM.py --in all_1000_genomes_ACB_processed --out all_1000_genomes_ACB_simulated --simulate_nbreakpoints 4 --simulate_nsamples 10000 --population_code ACB --human_genome_version hg19
```

Each argument does the following:
* --in: the filename prefix for the input (.bed, .bim, .fam) set of Plink files
* --out: the filename prefix for the output (.bed, .bim, .fam) set of Plink files
* --simulate_nbreakpoints: specifies the number of breakpoints to use for each chromosome. Setting it equal to n means that every chromosome is divided into n+1 chunks. Note that I used 4 breakpoints, but using fewer would be reasonable and also require less RAM and computation time. With that said, I recommend choosing at least 2 or 3 breakpoints to ensure diversity because some of the chromosomes have a single recombination interval with a roughly 20% chance of being sampled as a breakpoint.   
* --simulate_nsamples: the number of samples to simulate.
* --population_code: the 1000 genomes project population code that is closest to the population of the input fileset. In this case, they are the same exact population. 
* --human_genome_version: either hg19 or hg38, depending on the reference genome to which your input dataset was mapped. All provided 1000 genomes datasets were mapped to hg19. 

# How to Simulate Precise Genotype Phenotype correlations

REALGenomeSIM can simulate correlations between a binary or continuous phenotype and any linear combination of products of f(SNPs), where f is one of five possible functions that transform the original SNP values in biologically plausible ways:

* (regular, minor) = <img src="https://render.githubusercontent.com/render/math?math=I:\{0, 1, 2\} \rightarrow \{0, 1, 2\}">. All SNP values stay the same
* (recessive, minor) = <img src="https://render.githubusercontent.com/render/math?math=R:\{0, 1, 2\} \rightarrow \{0, 0, 2\}">. All 1 values become 0, so only homozygous minor allele genotypes have an effect
* (dominant, minor) = <img src="https://render.githubusercontent.com/render/math?math=D:\{0, 1, 2\} \rightarrow \{0, 2, 2\}">. All 1 values become 2, so both homozygous minor allele genotypes and heterozygous genotypes have equivalent effects
* (heterozygous_only, minor) = <img src="https://render.githubusercontent.com/render/math?math=He:\{0, 1, 2\} \rightarrow \{0, 2, 0\}">. All original 2 values become 0, and then all 1 values become 2. Only heterozygous allele genotypes have an effect
* (homozygous_only, minor) = <img src="https://render.githubusercontent.com/render/math?math=Ho:\{0, 1, 2\} \rightarrow \{2, 0, 2\}">. All original 0 values become 2, and then all 1 values become 0. Only homozygous allele genotypes have an effect

REALGenomeSIM can also specify that the major allele is associated with an effect instead of the minor allele, which is equivant to setting <img src="https://render.githubusercontent.com/render/math?math=f_{swap}:\{0, 1, 2\} \rightarrow \{2, 1, 0\}">. Using this option prior to using either of the five functions listed above results in the following transformations:

* (regular, major) = <img src="https://render.githubusercontent.com/render/math?math=I \circ f_{swap}:\{0, 1, 2\} \rightarrow I:\{2, 1, 0\} \rightarrow \{2, 1, 0\}">.
* (recessive, major) = <img src="https://render.githubusercontent.com/render/math?math=R \circ f_{swap}:\{0, 1, 2\} \rightarrow R:\{2, 1, 0\} \rightarrow \{2, 0, 0\}">.
* (dominant, major) = <img src="https://render.githubusercontent.com/render/math?math=D \circ f_{swap}:\{0, 1, 2\} \rightarrow D:\{2, 1, 0\} \rightarrow \{2, 2, 0\}">. 
* (heterozygous_only, major) = <img src="https://render.githubusercontent.com/render/math?math=He \circ f_{swap}:\{0, 1, 2\} \rightarrow He:\{2, 1, 0\} \rightarrow \{0, 2, 0\}">. 
* (homozygous_only, major) = <img src="https://render.githubusercontent.com/render/math?math=Ho \circ f_{swap}:\{0, 1, 2\} \rightarrow Ho:\{2, 1, 0\} \rightarrow \{2, 0, 2\}">.

### Example 1: a simple additive model

let <img src="https://render.githubusercontent.com/render/math?math=y"> be an individual's phenotype, <img src="https://render.githubusercontent.com/render/math?math=S_i^m"> be the <img src="https://render.githubusercontent.com/render/math?math=i^{th}"> to influence the value of <img src="https://render.githubusercontent.com/render/math?math=y"> such that the minor allele equals 1, and <img src="https://render.githubusercontent.com/render/math?math=B"> be the bias term. The goal is to simulate the following relationship between genotype and phenotype:
  
<img src="https://render.githubusercontent.com/render/math?math=y = 0.1S_1^m %2B 0.1S_2^m %2B 0.1S_3^m %2B 0.1S_4^m %2B 0.1S_5^m %2B 0.1S_6^m %2B 0.1S_7^m %2B 0.1S_8^m %2B 0.1S_9^m %2B 0.1S_10^m %2B B %2B \epsilon">

<img src="https://render.githubusercontent.com/render/math?math=\epsilon ~ N(\mu = 0, \sigma_{\epsilon} = 0.5E[y])">

<img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75">

The following files, formated as follows, must must exist in your working directory (you can name them as you please, which must be specified later.)

REALGenomeSIM_main_causal_SNP_IDs.txt contains specified SNP IDs from the input bim file all_1000_genomes_ACB_processed.bim seperated by newline characters

```
rs113633859
rs6757623
rs5836360
rs35542336
rs34342515
rs1867634
rs5004086
rs10883077
rs2852253
rs5801463
```

For each SNP ID, REALGenomeSIM_main_betas.txt contains a corresponding beta coefficient.

```
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
```

A full command for REALGenomeSIM to simulate genomic data with correlated phenotypes would be formatted as follows (noting that REALGenomeSIM_main_SNP_phenotype_map.txt is not necessary at this point). 

```
python REALGenomeSIM_main.py --in all_1000_genomes_ACB_processed --out all_1000_genomes_ACB_simulated --simulate_nbreakpoints 4 --simulate_nsamples 10000 --phenotype continuous --mean_phenotype 5.75 --population_code ACB --human_genome_version hg19 --noise 0.5 --causal_SNP_IDs_path REALGenomeSIM_main_causal_SNP_IDs.txt  --betas_path REALGenomeSIM_main_betas.txt
```

Each argument that wasn't previously explained does the following:
* --phenotype: can be either continuous (for linear regression) or binary (for logistic regression)
* --mean_phenotype: the desired value of <img src="https://render.githubusercontent.com/render/math?math=E[y]">. REALGenomeSIM modified <img src="https://render.githubusercontent.com/render/math?math=B"> so that <img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75">. If the phenotype is binary, then <img src="https://render.githubusercontent.com/render/math?math=E[y]"> is the proportion of cases and is thusly constrained to reside in the interval <img src="https://render.githubusercontent.com/render/math?math=0 < E[y] < 1\]">. I would recommend setting <img src="https://render.githubusercontent.com/render/math?math=0.05 < E[y] < 0.95">, particularly since most real GWAS datasets have at least 5% of the less frequent binary status. 
* --noise: a percentage of <img src="https://render.githubusercontent.com/render/math?math=E[y]"> that <img src="https://render.githubusercontent.com/render/math?math=\sigma_{\epsilon}"> is equal to. In the first model equation, <img src="https://render.githubusercontent.com/render/math?math=\sigma_{\epsilon} = 0.5*5.75 = 2.875">.
* --causal_SNP_IDs_path: the name of the file formatted as REALGenomeSIM_main_causal_SNP_IDs.txt or the full path if its not in the working directory
* --betas_path: the name of the file formatted as REALGenomeSIM_main_betas.txt or the full path if its not in the working directory
* --major_minor_assignments_path: the name of the file formatted as REALGenomeSIM_main_major_minor_assignments.txt or the full path if its not in the working directory

### Example 2: inclusion of nonlinear single-SNP effects

In addition to the notation from the first example, let <img src="https://render.githubusercontent.com/render/math?math=S_i^M = f_{swap}(S_i^m)"> be the <img src="https://render.githubusercontent.com/render/math?math=i^{th}"> to influence the value of <img src="https://render.githubusercontent.com/render/math?math=y"> such that the major allele equals 1. Also recall the definitions for the four nontrivial mapping functions (R, D, He, Ho) defined prior to the first example. The second example will model phenotypes as follows:

<img src="https://render.githubusercontent.com/render/math?math=y = 0.1S_1^m %2B 0.1R(S_2^m) %2B 0.1D(S_3^m) %2B 0.1He(S_4^m) %2B 0.1Ho(S_5^m) %2B 0.1S_6^M %2B 0.1R(S_7^M) %2B 0.1D(S_8^M) %2B 0.1He(S_9^M) %2B 0.1Ho(S_10^M) %2B B %2B \epsilon">

For each SNP ID, REALGenomeSIM_main_major_minor_assignments.txt specifies whether the minor allele or the major allele equals 1, and correspondingly, whether homozygous minor or homozygous major equals 2. In REALGenomeSIM_main_major_minor_assignments.txt, values of 0 mean the minor allele equals 1, and values of 1 mean that the major allele equals 1. Recall that the first five SNPs should be minor and the last 5 SNPs should be major:

```
0
0
0
0
0
1
1
1
1
1
```

There is a fourth optional document REALGenomeSIM_main_SNP_phenotype_map.txt that specifies whether the SNP's effect is regular, recessive, dominant, heterozygous_only, or homozygous_only. You can ignore this document if you want all of the SNPs to have regular effects (as in the first model equation), but an equivalent REALGenomeSIM_main_SNP_phenotype_map.txt would be formatted as follows:

```
regular
regular
regular
regular
regular
regular
regular
regular
regular
regular
```
