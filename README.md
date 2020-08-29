# About REALGenomeSIM

**REALGenomeSIM** simulates GWAS data with **R**ealistic **E**pistatic effects, **A**dditive effects, and **L**D patterns with **Genome** **S**ubsampling and **I**ntegrated **M**odeling. It divides each simulated genome into chunks such that the probability that any genomic position demarates two chunks equals the probability that this position would demarcate a real (given) recombination event. Each chunk is then filled with the homologous genomic subset of a real whole genome that is sampled uniformly (with replacement) from an input plink bed file. Since REALGenomeSIM's in-silico recombinations occur in the same proportions as equivalent real recombination events, each simulated individual is likely similar to at least one member of the input population whose genome is not included in the data. 

<p align="center">
<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=750/>
</p>

The Triadsim algorithm has used this method to simulate LD patterns that are almost indistinguishable from those of the input Dataset. REALGenomeSIM simulates equally realistic data, and it has the following improvements:

1. REALGenomeSIM is approximately 50 times faster than Triadsim.

<p align="center">
<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
</p>

2. REALGenomeSIM can use any whole genomes as an input dataset. This allows REALGenomeSIM to simulate new individuals from any publically available GWAS dataset, including non-human organisms that do not reproduce sexually. The following figures compare the publically available 1000 genomes project subpopulations to datasets that Triadsim simulated from each subpopulation. 

<p align="center">
<img src="https://github.com/EpistasisLab/REALGenomeSIM/blob/master/images/REALGenomeSIM.png" width=500/>
</p>

3. REALGenomeSIM uses the output recombination rate maps from () and converts them to probabilities of drawing each simulated breakpoint at a specific genomic location. Depending on your intent, the following optiona are available:

..*  If you want to simulate arbitrary datasets of realistic human genomes, then this package includes all samples from the 1000 genomes project that I have filtered down to 500000 SNPs. In summary, I kept every biallelic SNP such that every subpopulation contains at least two instances of the minor allele. Exact thinning methods are here. 

..*  If you want to simulate realistic human genomes that are similar to a specific dataset, then you should choose the 1000 genomes project subpopulation that is closest to the population of your sample in the () argument. 

..*  If you want to simulate non-human genomes, then you will need to find external recombination rate values for all genomic regions that you intend to simulate. This information will need to be formatted as follows:

# How to Simulate Whole Genomes without Simulating Phenotypes 

# How to Simulate Precise Genotype Phenotype correlations
