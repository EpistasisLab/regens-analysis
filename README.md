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
    <summary>These figures compare TSNE plots of the first 10 principal components for real and simulated 1000 genomes subpopulations.</summary>
  <p align="center">
  <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
  </p>
  <p align="center">
  <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
  </p>
  <p align="center">
  <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
  </p>
  <p align="center">
  <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
  </p>
  <p align="center">
  <img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
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
pip install git+https://github.com/fastlmm/PySnpTools.git@9d8eed0a
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
python triadsim_main.py --bfile all_1000_genomes_ACB_processed --out all_1000_genomes_ACB_simulated --simulate-nbreakpoints 2 --simulate-nsamples 10000 --population_code ACB --human_genome_version hg19
```

Each argument does the following:
* --bfile: the filename prefix for the input (.bed, .bim, .fam) set of Plink files
* --out: the filename prefix for the output (.bed, .bim, .fam) set of Plink files
* --simulate-nbreakpoints: specifies the number of breakpoints to use for each chromosome. Setting it equal to n means that every chromosome is divided into n+1 chunks. 
* --simulate-nsamples: the number of samples to simulate.
* --population_code: the 1000 genomes project population code that is closest to the population of the input fileset. In this case, they are the same exact population. 
* --human_genome_version: either hg19 or hg38, depending on the reference genome to which your input dataset was mapped. All provided 1000 genomes datasets were mapped to hg19. 

# How to Simulate Precise Genotype Phenotype correlations
