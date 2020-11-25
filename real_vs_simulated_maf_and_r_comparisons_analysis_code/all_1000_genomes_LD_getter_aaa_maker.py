subpops = ['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'JPT', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']

for pop in subpops:
    file = open("all_1000_genomes_LD_getter_" + pop + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#BSUB -J all_1000_genomes_LD_getter\n")
    file.write("#BSUB -o all_1000_genomes_LD_getter.out\n")
    file.write("#BSUB -e all_1000_genomes_LD_getter.error\n")
    file.write('#BSUB -R "rusage[mem=50000MB]"\n')
    file.write("#BSUB -M 50000MB\n")
    file.write("source activate regens\n")
    file.write("module load plink/1.90Beta6.18\n\n")

    for chr in range(1,23):
        file.write("plink --bfile ../real_1000_genomes_input_for_analysis/all_1000_genomes_" + pop + "_processed --chr " + str(chr))
        file.write(" --make-bed --out all_1000_genomes_" + pop + "_processed_chr" + str(chr) + "\n")
    
    merge_file_name = "all_1000_genomes_" + pop + "_processed_merge.txt"
    merge_file = open(merge_file_name, "w")
    for chr in range(2,23):
        merge_file.write("all_1000_genomes_" + pop + "_processed_chr" + str(chr) + "\n")
    merge_file.close()
    file.write("\nplink --bfile all_1000_genomes_" + pop + "_processed_chr1 --merge-list ")
    file.write(merge_file_name + " --make-bed --out all_1000_genomes_" + pop + "_processed\n\n")

    for chr in range(1,23):
        file.write("plink --bfile all_1000_genomes_" + pop + "_processed --chr " + str(chr))
        file.write(" --make-bed --out all_1000_genomes_" + pop + "_processed_chr" + str(chr) + "\n")

    file.write("\npython -m regens --in ../real_1000_genomes_input_for_analysis/all_1000_genomes_" + pop + "_processed ")
    file.write("--out all_1000_genomes_" + pop + "_simulated --simulate_nbreakpoints 4 ")
    file.write("--simulate_nsamples 100000 --population_code " + pop + " --human_genome_version hg19\n\n")

    for chr in range(1,23):
        file.write("plink --bfile all_1000_genomes_" + pop + "_processed_chr" + str(chr))
        file.write(" --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0")
        file.write(" --out all_1000_genomes_" + pop + "_processed_chr" + str(chr) + "\n")

    file.write("\n")
    for chr in range(1,23):
        file.write("plink --bfile all_1000_genomes_" + pop + "_simulated_chr" + str(chr))
        file.write(" --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0")
        file.write(" --out all_1000_genomes_" + pop + "_simulated_chr" + str(chr) + "\n")

    file.write("\n")
    file.write("python REALGenomeSIM_LD_getter.py --ref all_1000_genomes_" + pop)
    file.write("_processed --sim all_1000_genomes_" + pop + "_simulated")

    file.write("\n\n")

    file.write("rm all_1000_genomes_" + pop + "_processed.bim\n")
    file.write("rm all_1000_genomes_" + pop + "_processed.fam\n")
    file.write("rm all_1000_genomes_" + pop + "_processed.bed\n")
    file.write("rm all_1000_genomes_" + pop + "_processed.log\n")
    for chr in range(1,23):
        for end in [".bed\n", ".bim\n", ".fam\n", ".log\n", ".ld\n"]:
            file.write("rm all_1000_genomes_" + pop + "_processed_chr" + str(chr) + end)
            file.write("rm all_1000_genomes_" + pop + "_simulated_chr" + str(chr) + end)


    file.close()