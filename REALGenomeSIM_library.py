import numpy as np
import pandas as pd
import pdb
from time import time
from bed_reader import open_bed
from bed_reader import to_bed
from scipy.optimize import root
from copy import deepcopy as COPY
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
pd.options.mode.chained_assignment = None
#from triadsim_testers import test_drawn_breakpoints
#from triadsim_testers import test_breakpoint_SNP_mapping

def reduce_recomb_rate_info(rcmb_rate_info, bim_SNP_positions):
    
    rcmb_rate_intervals = rcmb_rate_info["Position(bp)"].to_numpy()
    SNP_pos_rcmb_interval_map = SNP_positions_to_rcmb_intervals(rcmb_rate_intervals, 
                                                                COPY(bim_SNP_positions), 
                                                                np.zeros(len(bim_SNP_positions), dtype = np.int64),
                                                                warnings = True)
    all_rcmb_intervals = np.arange(len(rcmb_rate_intervals))
    occupied_rcmb_intervals = np.isin(all_rcmb_intervals, SNP_pos_rcmb_interval_map)
    occupied_rcmb_intervals[np.min(np.where(occupied_rcmb_intervals == True)) - 1] = True
    reduced_rcmb_rate_info = rcmb_rate_info[occupied_rcmb_intervals]
    new_positions = reduced_rcmb_rate_info["Position(bp)"].to_numpy()
    new_lengths = new_positions[1:] - new_positions[0:-1]
    new_cumulative_cM = reduced_rcmb_rate_info["Map(cM)"].to_numpy()
    new_cM = new_cumulative_cM[1:] - new_cumulative_cM[0:-1]
    new_rcmb_rate_intervals = np.append(new_cM*1000000/new_lengths, [0])
    reduced_rcmb_rate_info.loc[:, "Rate(cM/Mb)"] = new_rcmb_rate_intervals
    return(reduced_rcmb_rate_info)

def SNP_positions_to_rcmb_intervals(rcmb_interval_boundaries, SNP_positions, SNP_rcmb_interval_indices, warnings):
    prev_start_pos = 0 
    too_large_indices = SNP_positions >= np.max(rcmb_interval_boundaries)
    too_small_indices = SNP_positions <= np.min(rcmb_interval_boundaries)
    SNP_positions[too_large_indices] = np.max(rcmb_interval_boundaries) - 1
    SNP_positions[too_small_indices] = np.min(rcmb_interval_boundaries) + 1
    num_outer_SNPs = np.sum(too_large_indices) + np.sum(too_small_indices) 
    if num_outer_SNPs > 0 and warnings == True:
        print("\nCAUTION: " + str(num_outer_SNPs) + "SNP(s) in the bim file are outside of the positional range considered by the recombination rate info file.\n")
        print("CAUTION: bim file SNPs earlier than the first recombination interval are assumed to be inside of the first recombination interval.\n")
        print("CAUTION: bim file SNPs farther than the last recombination interval are assumed to be inside of the last recombination interval.\n")
        print("CAUTION: This is unlikely to be an issue if the number of such SNPs is small. Otherwise, remove those SNPs from the input plink files, or find a more comprehensive recombination rate file\n")
    for i in range(len(SNP_positions)):
        for j in range(len(rcmb_interval_boundaries) - prev_start_pos):
            if rcmb_interval_boundaries[prev_start_pos + j] >= SNP_positions[i]:
                SNP_rcmb_interval_indices[i] = prev_start_pos + j
                prev_start_pos += j
                break
    return SNP_rcmb_interval_indices

def choice_with_periodic_replacement(length, width, probabilities):
    samples = np.zeros((width, length)).astype(int)
    cumulative_probabilities = np.repeat([np.cumsum(np.append([0], probabilities))], length, axis = 0)
    max_cumulative_vals = cumulative_probabilities[:, -1]
    for i in range(width):
        adjusted_uniform_sample = max_cumulative_vals*np.random.rand(length)
        samples[i] = np.array([np.searchsorted(probs, sample) for probs, sample in zip(cumulative_probabilities, adjusted_uniform_sample)])
        if i < (width - 1):
            for j, pos in enumerate(samples[i]):
                cumulative_probabilities[j][pos:] -= cumulative_probabilities[j][pos] - cumulative_probabilities[j][pos - 1]
            max_cumulative_vals = cumulative_probabilities[:, -1]
    return(samples.T - 1)
     
def centimorgans_to_probabilities(recomb_rate_info):

    # imports recombination rates and modifies them to centimorgan units.
    # The first recombination rate (in cM/mB) is measured is from the first position to the second position.
    # The recombination rate in the row of the last position is 0 cM/m because no position comes after it. 
    # This is the format by which there must be n positions and n - 1 intervals. 

    positions = recomb_rate_info["Position(bp)"].to_numpy()
    non_standard_rcmb_rates = (recomb_rate_info["Rate(cM/Mb)"].to_numpy())[0:-1]
    lengths = positions[1:] - positions[0:-1]
    rcmb_rates = non_standard_rcmb_rates*lengths/1000000

    # converts centimorgan units into into p(Ri = True|R = 1) values
    expected_Ri = (1 - np.exp(-rcmb_rates/50))/2
    expected_R = np.sum(expected_Ri)
    p_Ri_true_given_R_one = expected_Ri/expected_R
    return(p_Ri_true_given_R_one)

def draw_breakpoints(rcmb_rate_info, bim_SNP_pos, num_breakpoints, simulation_sample_size):
    SNP_count = len(bim_SNP_pos)
    probabilities = centimorgans_to_probabilities(rcmb_rate_info)
    rcmb_rate_intervals = rcmb_rate_info["Position(bp)"].to_numpy()
    breakpoints = choice_with_periodic_replacement(simulation_sample_size,  
                                                   num_breakpoints, 
                                                   probabilities)

    # NOTE: the values of rcmb_rate_intervals are genomic positions of the intervals' left bounds
    # NOTE: the first value in "probabilities" is the probability of drawing the interval in between indices 0 and 1.
    # NOTE: SNP_pos_rcmb_interval_map maps each SNP index to the left bound of the interval containing the SNP position.
    # NOTE: It is also impossible to select the last index in the cumulative distribution, so if the last index is 18021, the max is 18020. 
    # NOTE: the domain of original breakpoints is [0, last_index - 1], so I subtract 1 from the result below. 

    SNP_pos_rcmb_interval_map = SNP_positions_to_rcmb_intervals(rcmb_rate_intervals, 
                                                                COPY(bim_SNP_pos), 
                                                                np.zeros(SNP_count, dtype = np.int64),
                                                                warnings = False) - 1
   
    rcmb_interval_SNP_pos_map = {}
    for rcmb_interval in np.unique(SNP_pos_rcmb_interval_map):
        rcmb_interval_SNP_pos_map[rcmb_interval] = np.where(SNP_pos_rcmb_interval_map == rcmb_interval)[0]
        
    for jj in range(len(breakpoints)):
        for k in range(num_breakpoints):
            interval_index = breakpoints[jj][k]
            SNP_indices = rcmb_interval_SNP_pos_map[interval_index]
            if len(SNP_indices) == 1:
                breakpoints[jj][k] = SNP_indices[0]
            else:
                breakpoints[jj][k] = SNP_indices[int(len(SNP_indices)*np.random.rand() - 0.5)]

    return(breakpoints)

def get_samples_fast(simulation_sample_size, bed_col_bounds, plink_file_name_prefix, num_breakpoints):
    bed_file_path = plink_file_name_prefix + ".bed"
    bed_reader = open_bed(bed_file_path, count_A1 = True, num_threads = 1)
    num_repeats = int(simulation_sample_size*(num_breakpoints + 1)/int(bed_reader.iid_count) + 1)
    bed_row_indices = np.repeat(range(bed_reader.iid_count), num_repeats)
    np.random.shuffle(bed_row_indices)
    bed_row_indices = bed_row_indices[:simulation_sample_size*(num_breakpoints + 1)]
    bed_file_samples = bed_reader.read((bed_row_indices, slice(bed_col_bounds[0], bed_col_bounds[1])), dtype = 'int8')
    bed_file_samples_dimensions = (int(len(bed_row_indices)/(num_breakpoints + 1)),  
                                   num_breakpoints + 1, 
                                   len(bed_file_samples[0]))
    return(bed_file_samples.reshape(bed_file_samples_dimensions))

def breakpoints_to_simulated_individuals(breakpoints, SNP_count, sampled_individuals):
    breakpoints = np.append(breakpoints, SNP_count*np.ones((len(breakpoints), 1), dtype = np.int64), axis = 1)
    breakpoints_masks = np.zeros((len(breakpoints), len(breakpoints[0]), SNP_count), dtype = bool)
    for j, breakpoint_vec in enumerate(breakpoints):
        startpoint = 0
        for k, breakpoint in enumerate(np.sort(breakpoint_vec)):
            breakpoints_masks[j][k][startpoint:breakpoint] = True
            startpoint = breakpoint
    simulated_individuals = sampled_individuals[breakpoints_masks].reshape(-1, SNP_count)
    return(simulated_individuals)

def write_bed_file(simulated_individuals, bim_SNP_names, output_name, bim_SNP_complete_pos, bim_SNP_nucleotides, population_ID):
    simulated_IDs = np.array([population_ID + "_" + str(i) for i in range(1, len(simulated_individuals) + 1)])
    metadata = {"fid": simulated_IDs,
                "iid": simulated_IDs,
                "sex": np.array([2]*len(simulated_IDs)),
                "pheno": np.array([-9]*len(simulated_IDs)),
                "chromosome": bim_SNP_complete_pos.T[0],
                "sid": bim_SNP_names,
                "cm_position": bim_SNP_complete_pos.T[1],
                "bp_position": bim_SNP_complete_pos.T[2],
                "allele_1": bim_SNP_nucleotides.T[0],
                "allele_2": bim_SNP_nucleotides.T[1]}
    to_bed(output_name, simulated_individuals, properties = metadata, count_A1 = True)
        
def get_feature_SNPs(feature_SNP_IDs, cumulative_SNP_counts, output_file_names, sample_size, bim_SNP_names):
    feature_SNPs = np.zeros((sample_size, len(feature_SNP_IDs)))
    for i, SNP_ID in enumerate(feature_SNP_IDs):
        index = np.where(SNP_ID == bim_SNP_names)[0][0]
        output_file_index = np.min(np.where(index <= cumulative_SNP_counts)) - 1
        output_file_reader = open_bed(output_file_names[output_file_index], count_A1 = True, num_threads = 1)
        adjusted_SNP_index = index - cumulative_SNP_counts[output_file_index]
        feature_SNPs[:, [i]] += output_file_reader.read(index = adjusted_SNP_index, dtype = 'float32')
    return(feature_SNPs)

def simulate_phenotypes(output_file_names, causal_SNP_IDs_path, cumulative_SNP_counts, 
                        major_minor_assignments_path, betas_path, mean_phenotype, sample_size,
                        bim_SNP_names, phenotype, SNP_phenotype_map_path, noise = 0):

    # imports required model components (SNPs and beta values).
    github_link = "blahblahblah.com"
    causal_SNP_IDs = open(causal_SNP_IDs_path, "r").readlines()
    try:
        betas = np.array(open(betas_path, "r").readlines()).astype(np.float64)
    except:
        print("\nerror: The beta coefficients file at " + betas_path + " is incorrectly formatted. Visit " + github_link + " for examples of correct formatting.\n")
        exit()
    if len(betas) != len(causal_SNP_IDs):
        print("\nerror: The causal_SNP_IDs and betas files must have the same number of rows. Visit " + github_link + " for examples of correct formatting.\n")
        exit()

    # imports optional model components (major/minor assignments and SNP_phenotype_maps).
    if major_minor_assignments_path != "standard":
        major_minor_assignments = open(major_minor_assignments_path, "r").readlines()
        if len(major_minor_assignments) != len(causal_SNP_IDs):
            print("\nerror: The causal_SNP_IDs and major_minor_assignments files must have the same number of rows. Visit " + github_link + " for examples of correct formatting.\n")
            exit()
    if SNP_phenotype_map_path != "standard":
        SNP_phenotype_map = open(SNP_phenotype_map_path, "r").readlines()
        if len(SNP_phenotype_map) != len(causal_SNP_IDs):
            print("\nerror: The causal_SNP_IDs and SNP_phenotype_map files must have the same number of rows. Visit " + github_link + " for examples of correct formatting.\n")
            exit()

    # simulates phenotypes based on model specifications
    feature_size = len(betas)
    features = np.zeros((sample_size, feature_size))
    for p in range(feature_size):
        feature_SNP_IDs = causal_SNP_IDs[p].strip().split('\t')
        try:
            feature_SNPs = get_feature_SNPs(feature_SNP_IDs, cumulative_SNP_counts, output_file_names, sample_size, bim_SNP_names)
        except:
            print("\nerror: The causal SNP IDs on row " + str(p + 1) + " are either incorrectly formatted or they do not exist in the input bim file:\n")
            print("Visit " + github_link + " for examples of correct formatting.\n")
        if major_minor_assignments_path != "standard":
            if np.all(np.isin(major_minor_assignments[p].strip().split('\t'), ['0', '1'])):            
                feature_major_minor_assignments = np.array(major_minor_assignments[p].strip().split('\t')).astype(np.int64)
                feature_major_minor_assignments_alt = COPY(feature_major_minor_assignments)
                feature_major_minor_assignments_alt[feature_major_minor_assignments_alt == 0] = -1
                feature_SNPs_with_assignments = (feature_SNPs - 2*feature_major_minor_assignments)*(-1*feature_major_minor_assignments_alt)
            else:
                print("\nerror: The major minor assignments on row " + str(p + 1) + " are either incorrectly formatted or they do not exist in the input bim file:\n")
                print("Visit " + github_link + " for examples of correct formatting.\n")
        if major_minor_assignments_path == "standard":
            feature_SNPs_with_assignments = feature_SNPs

        if SNP_phenotype_map_path != "standard":
            feature_SNP_phenotype_map = SNP_phenotype_map[p].strip().split('\t')
            if feature_SNPs_with_assignments.shape[1] > 1:
                for m in range(len(feature_SNP_phenotype_map)):
                    if feature_SNP_phenotype_map[m] == "recessive":
                        feature_SNPs_with_assignments[:, m][feature_SNPs_with_assignments[:, m] == 1] = 0
                    elif feature_SNP_phenotype_map[m] == "dominant":
                        feature_SNPs_with_assignments[:, m][feature_SNPs_with_assignments[:, m] == 1] = 2
                    elif feature_SNP_phenotype_map[m] == "heterozygous_only":
                        feature_SNPs_with_assignments[:, m][feature_SNPs_with_assignments[:, m] == 2] = 0
                        feature_SNPs_with_assignments[:, m][feature_SNPs_with_assignments[:, m] == 1] = 2
                    elif feature_SNP_phenotype_map[m] == "homozygous_only":
                        feature_SNPs_with_assignments[:, m][feature_SNPs_with_assignments[:, m] == 0] = 2
                        feature_SNPs_with_assignments[:, m][feature_SNPs_with_assignments[:, m] == 1] = 0
                    elif feature_SNP_phenotype_map[m] == "regular":
                        pass
                    else:
                        print("\nerror: all SNP_phenotype labels must be 'regular', 'recessive', 'dominant', 'heterozygous_only', or 'homozygous_only'.\n")
                        print("Visit " + github_link + " for examples of correct formatting.\n")
                        exit()
            if feature_SNPs_with_assignments.shape[1] == 1:
                for m in range(len(feature_SNP_phenotype_map)):
                    if feature_SNP_phenotype_map[m] == "recessive":
                        feature_SNPs_with_assignments[feature_SNPs_with_assignments == 1] = 0
                    elif feature_SNP_phenotype_map[m] == "dominant":
                        feature_SNPs_with_assignments[feature_SNPs_with_assignments == 1] = 2
                    elif feature_SNP_phenotype_map[m] == "heterozygous_only":
                        feature_SNPs_with_assignments[feature_SNPs_with_assignments == 2] = 0
                        feature_SNPs_with_assignments[feature_SNPs_with_assignments == 1] = 2
                    elif feature_SNP_phenotype_map[m] == "homozygous_only":
                        feature_SNPs_with_assignments[feature_SNPs_with_assignments == 0] = 2
                        feature_SNPs_with_assignments[feature_SNPs_with_assignments == 1] = 0
                    elif feature_SNP_phenotype_map[m] == "regular":
                        pass
                    else:
                        print("\nerror: all SNP_phenotype labels must be 'regular', 'recessive', 'dominant', 'heterozygous_only', or 'homozygous_only'.\n")
                        print("Visit " + github_link + " for examples of correct formatting.\n")
                        exit()

        features[:, p] = np.product(feature_SNPs_with_assignments, axis = 1)
    weighted_feature_sums = np.sum(betas*features, axis = 1, keepdims = True) 
    weighted_feature_sums += np.random.normal(loc = 0, scale = noise*np.mean(weighted_feature_sums), size = weighted_feature_sums.shape)
    
    if phenotype == "binary":
        def logistic_with_unknown_intercept(intercept, weighted_feature_sums, mean_phenotype):
            disease_probabilities = 1/(1 + np.exp(-1*(weighted_feature_sums + intercept)))
            return(np.mean(disease_probabilities) - mean_phenotype)
        intercept = root(fun = logistic_with_unknown_intercept, x0 = np.array([0]), args = (weighted_feature_sums,  mean_phenotype)).x[0]
        disease_probabilities = 1/(1 + np.exp(-1*(weighted_feature_sums + intercept)))
        simulated_phenotypes = (np.random.rand(len(disease_probabilities)) <= disease_probabilities.reshape(-1)).astype(np.int8)
        model = LogisticRegression(C = 1E100, tol = 1E-100, max_iter = 1000000).fit(features, simulated_phenotypes)

    elif phenotype == "continuous":
        def linear_with_unknown_intercept(intercept, weighted_feature_sums, mean_phenotype):
            return(np.mean(weighted_feature_sums + intercept) - mean_phenotype)
        intercept = root(fun = linear_with_unknown_intercept, x0 = np.array([0]), args = (weighted_feature_sums,  mean_phenotype)).x[0]
        simulated_phenotypes = weighted_feature_sums + intercept
        model = LinearRegression().fit(features, simulated_phenotypes)

    else:
        print("error: phenotype must be either 'binary' or 'continuous'.")
        exit()
 
    model_profile = open(output_file_names[0][:-8] + "model_profile.txt", "w")
    model_profile.write("measured R^2 of model fit: " + str(model.score(features, simulated_phenotypes)) + "\n")
    for i, b in enumerate(model.coef_[0]):
        model_profile.write("measured beta value of feature" + str(i+1) + ": " + str(b) + "\n")
    model_profile.write("measured beta value of intercept: " + str(model.intercept_[0]))
    model_profile.close()
    return(simulated_phenotypes)
  



    