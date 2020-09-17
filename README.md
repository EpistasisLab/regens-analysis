# About REALGenomeSIM

**REALGenomeSIM** simulates GWAS data with **R**ealistic **E**pistatic effects, **A**dditive effects, and **L**D patterns with **Genome** **S**ubsampling and **I**ntegrated **M**odeling. Each genome that REALGenomeSIM simulates starts out as an empty template of SNP positions without SNP values, which is divided into empty segments that are demarcated by breakpoints. The probability of drawing any genomic position for a given breakpoint is equal to the probability that this position would demarcate a given real recombination event. Once an empty simulated genome is segmented by breakpoints, the row indices of whole genome bed file rows from a real dataset are duplicated so that 1) there is one real individual for each empty segment and 2) every real individual is selected an equal number of times (minus 1 for each remainder sample if the number of segments is not divisible by the number of individuals). Then, for each empty segment, a whole genome is randomly selected without replacement from the set of genomes that correspond to the duplicated indices, and the empty simulated segment is filled with the the homologous segment from the sampled real genome. These steps are repeated for every empty simulated segment so that all of the empty simulated genomes are filled with real SNP values. This quasirandom selection minimizes maf variation between the simulated and real datasets and also maintains normal population level genetic variability by randomizing segment selection. Even though the randomly selected segments are independent from one-another, the simulated dataset will contain the input dataset's LD pattern because each breakpoint location is selected with the same probability that they would demarcate a given real recombination event (i.e. a real biological concatenation of two independent genomic segment).

The Triadsim algorithm has used this method to simulate LD patterns that are almost indistinguishable from those of the input Dataset. REALGenomeSIM simulates equally realistic data, and it was measured to be 88 times faster and require 8 times lower peak RAM than Triadsim. REALGenomeSIM is designed to easily simulate GWAS data from any of the 26 populations in the [1000 genomes project](https://www.cog-genomics.org/plink/2.0/resources), and a filtered subset of these subpopulations' genotype data is provided in the github in corresponding plink filesets. In summary, I kept every biallelic SNP such that every subpopulation contains at least two instances of the minor allele. Exact thinning methods are here. REALGenomeSIM converts output recombination rate maps from [pyrho](https://github.com/popgenmethods/pyrho) (which correspond to the twenty-six 1000 Genome populations on a one to one basis) into probabilities of drawing each simulated breakpoint at a specific genomic location. It is also possible to simulate GWAS data from a custom plink (bed, bim, bam) fileset or a custom recombination rate map (or both files can be custom). Note that recombination rate maps between populations within a superpopulation (i.e. british and italian) have pearson correlation coefficients of roughly 0.9 [(see figure 2B of the pyrho paper)](https://advances.sciencemag.org/content/advances/5/10/eaaw9206.full.pdf), so if a genotype dataset has no recombination rate map for the exact population, then map for a closely relatrf population should suffice. 

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

REALGenomeSIM_causal_SNP_IDs.txt contains specified SNP IDs from the input bim file all_1000_genomes_ACB_processed.bim seperated by newline characters

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

For each row containing one or more SNP IDs (examples 1 and 2 only contain one SNP ID per row), REALGenomeSIM_betas.txt contains a corresponding beta coefficient.

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

A full command for REALGenomeSIM to simulate genomic data with correlated phenotypes would be formatted as follows (noting that REALGenomeSIM_SNP_phenotype_map.txt is not necessary at this point). 

```
python REALGenomeSIM.py --in all_1000_genomes_ACB_processed --out all_1000_genomes_ACB_simulated --simulate_nbreakpoints 4 --simulate_nsamples 10000 --phenotype continuous --mean_phenotype 5.75 --population_code ACB --human_genome_version hg19 --noise 0.5 --causal_SNP_IDs_path REALGenomeSIM_causal_SNP_IDs.txt --betas_path REALGenomeSIM_betas.txt
```

Each argument that wasn't previously explained does the following:
* --phenotype: can be either continuous (for linear regression) or binary (for logistic regression)
* --mean_phenotype: the desired value of <img src="https://render.githubusercontent.com/render/math?math=E[y]">. REALGenomeSIM modified <img src="https://render.githubusercontent.com/render/math?math=B"> so that <img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75">. If the phenotype is binary, then <img src="https://render.githubusercontent.com/render/math?math=E[y]"> is the proportion of cases and is thusly constrained to reside in the interval <img src="https://render.githubusercontent.com/render/math?math=0 < E[y] < 1\]">. I would recommend setting <img src="https://render.githubusercontent.com/render/math?math=0.05 < E[y] < 0.95">, particularly since most real GWAS datasets have at least 5% of the less frequent binary status. 
* --noise: a percentage of <img src="https://render.githubusercontent.com/render/math?math=E[y]"> that <img src="https://render.githubusercontent.com/render/math?math=\sigma_{\epsilon}"> is equal to. In the first model equation, <img src="https://render.githubusercontent.com/render/math?math=\sigma_{\epsilon} = 0.5*5.75 = 2.875">.
* --causal_SNP_IDs_path: the name of the file formatted as REALGenomeSIM_causal_SNP_IDs.txt or the full path if its not in the working directory
* --betas_path: the name of the file formatted as REALGenomeSIM_betas.txt or the full path if its not in the working directory
* --major_minor_assignments_path: the name of the file formatted as REALGenomeSIM_major_minor_assignments.txt or the full path if its not in the working directory

REALGenomeSIM also generates two files that characterize the simulated relationship between genotypes and phenptypes:

* The distribution of phenotypes (notice that it's normally distributed, which is expected under an additive model by the central limit theorem)

<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/example1_all_1000_genomes_ACB_simulated_phenotype_profile.png">

* A file containing the <img src="https://render.githubusercontent.com/render/math?math=R^2"> value of the phenotype/genotype correlation and the inferred beta coefficients (which will most likely be close to but not equal to the input beta coefficients).
```
measured R^2 of model fit: 0.21840881212783192
measured beta value of feature1: 0.10375469993274405
measured beta value of feature2: 0.09158917815769968
measured beta value of feature3: 0.09621670127541823
measured beta value of feature4: 0.10099616024893379
measured beta value of feature5: 0.10595962889622762
measured beta value of feature6: 0.09990693835458762
measured beta value of feature7: 0.09343695082947306
measured beta value of feature8: 0.1028013611752096
measured beta value of feature9: 0.09132853270999804
measured beta value of feature10: 0.10450671675638493
measured beta value of intercept: 4.92101167192938
```

### Example 2: inclusion of nonlinear single-SNP effects

In addition to the notation from the first example, let <img src="https://render.githubusercontent.com/render/math?math=S_i^M = f_{swap}(S_i^m)"> be the <img src="https://render.githubusercontent.com/render/math?math=i^{th}"> to influence the value of <img src="https://render.githubusercontent.com/render/math?math=y"> such that the major allele equals 1. Also recall the definitions for the four nontrivial mapping functions (R, D, He, Ho) defined prior to the first example. The second example will model phenotypes as follows:

<img src="https://render.githubusercontent.com/render/math?math=y = 0.1S_1^m %2B 0.1R(S_2^m) %2B 0.1D(S_3^m) %2B 0.1He(S_4^m) %2B 0.1Ho(S_5^m) %2B 0.1S_6^M %2B 0.1R(S_7^M) %2B 0.1D(S_8^M) %2B 0.1He(S_9^M) %2B 0.1Ho(S_10^M) %2B B %2B \epsilon">

<img src="https://render.githubusercontent.com/render/math?math=\epsilon ~ N(\mu = 0, \sigma_{\epsilon} = 0.5E[y])">

<img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75">

Specifying these components is optional and requires two additional optional documents, examples of which are named REALGenomeSIM_major_minor_assignments.txt and REALGenomeSIM_SNP_phenotype_map.txt.

For each SNP ID, REALGenomeSIM_major_minor_assignments.txt specifies whether the minor allele or the major allele equals 1, and correspondingly, whether homozygous minor or homozygous major equals 2. In REALGenomeSIM_major_minor_assignments.txt, values of 0 mean the minor allele equals 1, and values of 1 mean that the major allele equals 1. In accordance with the second example model equation, REALGenomeSIM_major_minor_assignments.txt specifies that the first five SNPs have their minor allele equal 1, and the last five SNPs have their major allele equal 1:

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

For each SNP ID, REALGenomeSIM_SNP_phenotype_map.txt specifies whether the SNP's effect is regular, recessive, dominant, heterozygous_only, or homozygous_only. You can ignore this document if you want all of the SNPs to have regular effects (as in the first model equation), but an equivalent REALGenomeSIM_SNP_phenotype_map.txt would be formatted as follows:

```
regular
dominant
recessive
heterozygous_only
homozygous_only
regular
dominant
recessive
heterozygous_only
homozygous_only
```

python REALGenomeSIM.py --in all_1000_genomes_ACB_processed --out all_1000_genomes_ACB_simulated --simulate_nbreakpoints 4 --simulate_nsamples 10000 --phenotype continuous --mean_phenotype 5.75 --population_code ACB --human_genome_version hg19 --noise 0.5 --causal_SNP_IDs_path REALGenomeSIM_causal_SNP_IDs.txt --major_minor_assignments_path REALGenomeSIM_major_minor_assignments.txt --betas_path REALGenomeSIM_betas.txt --SNP_phenotype_map_path REALGenomeSIM_SNP_phenotype_map.txt 

* --major_minor_assignments_path: the name of the file formatted as REALGenomeSIM_major_minor_assignments.txt or the full path if its not in the working directory
* --SNP_phenotype_map_path: the name of the file formatted as REALGenomeSIM_SNP_phenotype_map.txt or the full path if its not in the working directory

Here is the distribution of phenotypes.

<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/example2_all_1000_genomes_ACB_simulated_phenotype_profile.png">

Here is the file containing the <img src="https://render.githubusercontent.com/render/math?math=R^2"> value of the phenotype/genotype correlation and the inferred beta coefficients.
```
measured R^2 of model fit: 0.26062585493232215
measured beta value of feature1: 0.09532940357383092
measured beta value of feature2: 0.09838413538443105
measured beta value of feature3: 0.10356136855704572
measured beta value of feature4: 0.100622751774058
measured beta value of feature5: 0.09881877516455176
measured beta value of feature6: 0.10015672892458262
measured beta value of feature7: 0.09781134124143136
measured beta value of feature8: 0.10654323045769121
measured beta value of feature9: 0.10531031279802351
measured beta value of feature10: 0.10091331861759821
measured beta value of intercept: 4.800871008690798
```

### Example 3: inclusion of epistatic effects

REALGenomeSIM models epistasis between an arbitrary number of SNPs as the product of SNP values in an individual. Example three uses the following model:

<img src="https://render.githubusercontent.com/render/math?math=y = 0.1S_1^m %2B 0.1R(S_2^m) %2B 0.1D(S_3^m) %2B 0.1He(S_4^m) %2B 0.1Ho(S_5^m) %2B 0.1S_6^M %2B 0.1R(S_7^M) %2B 0.1D(S_8^M) %2B 0.1He(S_9^M) %2B 0.1Ho(S_10^M) %2B 0.15S_11^mS_12^m %2B 0.2S_13^mS_14^mD(S_15^M) %2B 0.3S_16^m(S_17^M)^2S_2^MD(S_4^m) %2B B %2B \epsilon">

<img src="https://render.githubusercontent.com/render/math?math=\epsilon ~ N(\mu = 0, \sigma_{\epsilon} = 0.5E[y])">

<img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75">

REALGenomeSIM_causal_SNP_IDs2.txt contains sets of specified SNP IDs from the input bim file all_1000_genomes_ACB_processed.bim, where the sets are seperated by newline characters. Each set (i.e. row) only contained one SNP in examples 1 and two, but rows corresponding to SNP products contain multiple tab seperated SNPs. A couple of other things are worth noticing:
* The <img src="https://render.githubusercontent.com/render/math?math=(S_17^M)^2"> factor in the last additive term squares the effect of that SNP. This is represented in REALGenomeSIM_causal_SNP_IDs2.txt by including that SNP twice. 
* <img src="https://render.githubusercontent.com/render/math?math=S_2 = rs6757623"> and <img src="https://render.githubusercontent.com/render/math?math=S_4 = rs35542336"> both appear twice. Also (if you look ahead to REALGenomeSIM_major_minor_assignments.txt), <img src="https://render.githubusercontent.com/render/math?math=rs35542336 = S_2^m"> in the monoallelic effect and <img src="https://render.githubusercontent.com/render/math?math=rs35542336 = S_2^M"> in the epistatic effect, demonstrating how the minor allele might cause a simulated monoallelic effect, while the major allele may participate in a different simulated epistatic effect.

```
rs113633859
rs6757623
rs5836360
rs35542336
rs34342515
rs7836990
rs1111818
rs7895376
rs35834549
rs5801463
rs11852537	rs10883077
rs1867634	rs545673871	rs10416530
rs2066224	rs62240045	rs62240045	rs35542336	rs6757623
```

For each row containing one or more SNP IDs, REALGenomeSIM_betas2.txt contains a corresponding beta coefficient. (Giving each SNP that participates in a multiplicative interaction its own beta coefficient would be pointless).

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
0.15
0.2
0.3
```

As before, both REALGenomeSIM_major_minor_assignments2.txt and REALGenomeSIM_SNP_phenotype_map.txt are both optional. If they are not specified, then all SNPs will have the minor allele equal 1 so that the homozygous minor genotype equals 2, and all SNP_phenotype maps will be regular. 

For each SNP ID, REALGenomeSIM_major_minor_assignments2.txt specifies whether the minor allele or the major allele equals 1, and correspondingly, whether homozygous minor or homozygous major equals 2. CAUTION: In general, if a SNP appears in two different effects, then it may safely have different major/minor assignments in different effects. However, if a SNP appears twice in the same effect, then make sure it has the same major/minor assignment within that effect, or that effect may equal 0 depending on the map functions that are used on the SNP. 

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
0	0
0	0	1
0	1	1	1	0
```
For each SNP ID, REALGenomeSIM_SNP_phenotype_map.txt specifies whether the SNP's effect is regular, recessive, dominant, heterozygous_only, or homozygous_only. In context to an epistatic interaction, first the map functions are applied to their respective SNPs, and then the resulting mapped SNP values are multiplied together as specified by REALGenomeSIM_causal_SNP_IDs2.txt

```
regular
dominant
recessive
heterozygous_only
homozygous_only
regular
dominant
recessive
heterozygous_only
homozygous_only
regular	regular
regular	regular	dominant
regular	regular	regular	regular	dominant
```

Here is the distribution of phenotypes (notice that epistatic interactions cause skewing because their effect sizes are asymetric. SNP value multiplication causes most values to equal 0, but a small percentage of them are very large.)

<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/example3_all_1000_genomes_ACB_simulated_phenotype_profile.png">

Here is the file containing the <img src="https://render.githubusercontent.com/render/math?math=R^2"> value of the phenotype/genotype correlation and the inferred beta coefficients.
```
measured R^2 of model fit: 0.7462774511213254
measured beta value of feature1: 0.09748408274096143
measured beta value of feature2: 0.1037345030486333
measured beta value of feature3: 0.10747729400786989
measured beta value of feature4: 0.10356087444337082
measured beta value of feature5: 0.08897293145796997
measured beta value of feature6: 0.114037362198526
measured beta value of feature7: 0.08469133818043673
measured beta value of feature8: 0.09587043261083873
measured beta value of feature9: 0.10945717146197069
measured beta value of feature10: 0.10210806247400407
measured beta value of feature11: 0.15534388074372998
measured beta value of feature12: 0.20040913226285006
measured beta value of feature13: 0.30116564210405655
measured beta value of intercept: 3.6428125632607453
```
