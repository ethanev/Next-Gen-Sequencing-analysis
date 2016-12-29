# Author: Ethan D Evans

# Welcome! This file provide a number of useful functions to run during an interactive python session which was done
# for all data in the accompanying manuscript. Many functions developed here can be generalized to other NGS data

import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
import copy
import re
from Bio import pairwise2
import matplotlib.ticker as ticker
from sklearn.preprocessing import normalize
import matplotlib.font_manager as font_manager

def find_values_greater_than(dictionary,value):
    '''
Return a list of keys and values for keys with values above a desired number
    '''
    new_dict = dict()
    for k in dictionary:
        dict_val = dictionary[k]
        if dict_val > value:
            new_dict[k] = dict_val
    return new_dict

def cys_dict(library):
    '''
    This function is split the library into a cys containing portions and one without cys
    '''
    C_dict = dict()
    No_C_dict = dict()
    for pep in library:
        if 'C' in pep:
            C_dict[pep] = library[pep]
        else:
            No_C_dict[pep] = library[pep]
    return C_dict, No_C_dict

def givencopynum_vs_cpnum(dictionary):
    '''
    This function is meant to return a dictionary in which the keys are the numbers of seqeuece
    counts observed and the values are the number of sequences with that count

    This was used to gather the data for Figure 2 in the manuscript
    '''
    copynum_dict = dict()
    for peptide in dictionary:
        if dictionary[peptide] not in copynum_dict:
            copynum_dict[dictionary[peptide]] = 1
        else:
            copynum_dict[dictionary[peptide]] += 1
    return copynum_dict

def AA_histo_maker(dictionary):
    """
This make a histogram out of all AA frequencies at each position in the library
    """
    AA_histo = dict()
    protein_list = []
    for key in dictionary:
        protein_list.append(key)
    for sequence in protein_list:
        for i in range(1,len(sequence)-1):
            if (sequence[i],i) in AA_histo:
                AA_histo[(sequence[i],i)] += 1
            else:
                AA_histo[(sequence[i],i)] = 1
    return AA_histo

def surrounding_cys(dictionary):
    '''
    This function is used to find short motifs around the cysteines, it can be changed by altering
    the length of peptide you wish to investigate (must change the range based on the motif size)
    with 3rd and 5th lines
    '''
    near_cys = dict()
    for peptide in dictionary:
        for i in range(0,len(peptide)-6):
            if peptide[i] =='C':
                seq = peptide[(i-1):(i+6)]
                if seq in near_cys:
                    near_cys[seq] += dictionary[peptide]
                else:
                    near_cys[seq] = dictionary[peptide]
    return near_cys

def levenshtein_test(seed, dictionary,value):
    """Intended to be a simple levenshtein sequence clustering, currently only looking at sequences that are the
    same length as the seed"""
    global family
    for goodpeptide in seed:
        if len(goodpeptide)>50:
            print("already here")
        else:
            seed_cluster = dict()
            for peptide in dictionary:
                counter = 0
                if len(peptide) < len(goodpeptide):
                    pass
                if len(peptide) == len(goodpeptide):
                    for i in range(0,len(peptide)-1):
                        if peptide[i] == goodpeptide[i]:
                            pass
                        else:
                            counter += 1
                    if counter <= value:
                        seed_cluster[peptide] = dictionary[peptide]
                if len(peptide) > len(goodpeptide):
                    buffer = len(peptide) - len(goodpeptide)
                    temp_seed = goodpeptide + buffer*"_"
                    for i in range(0,len(peptide)-1):
                        if peptide[i] == temp_seed[i]:
                            pass
                        else:
                            counter += 1
                    if counter <= value:
                        seed_cluster[peptide] = dictionary[peptide]
    return seed_cluster

def levn_cluster(cluster,dictionary,value):
    """
    Intended to be a simple levenshtein sequence clustering, currently only looking at sequences that are the
    same length as the seed, this will continue and perform levenshtein analysis on each of the
    peptides added to the list within a given edit distance
    """
    seed_cluster = dict()
    while cluster:
        seed = cluster.pop()
        print seed
        seed_cluster[seed] = dictionary[seed]
        print cluster
        for peptide in dictionary:
            counter = 0
            if len(peptide[:30]) < len(seed):
                pass
            if len(peptide[:30]) == len(seed):
                for i in range(0,len(peptide[:30])-1):
                    if peptide[i] == seed[i]:
                        pass
                    else:
                        counter += 1
                if counter <= value:
                    if peptide[:30] not in cluster and peptide[:30] not in seed_cluster:
                        cluster.append(peptide[:30])
            if len(peptide[:30]) > len(seed):
                buffer = len(peptide[:30]) - len(seed)
                temp_seed = seed + buffer*"_"
                for i in range(0,len(peptide[:30])-1):
                    if peptide[i] == temp_seed[i]:
                        pass
                    else:
                        counter += 1
                if counter <= value:
                    if peptide[:30] not in cluster and peptide[:30] not in seed_cluster:
                        cluster.append(peptide[:30])
    return seed_cluster, cluster

def count_cys_peptides(library):
    '''
    Will simply count the number of peptides with cysteines in them and also will find
    peptides with more than 5 cysteines - not mentioned in the manuscript, but interesting!
    Check out some of the high number cys containing peptides!
    '''
    count_dict = dict()
    very_unique_peps = dict()
    for peptide in library:
        count = 0
        for aa in peptide:
            if aa == 'C':
                count += 1
        if count in count_dict:
            count_dict[count] += 1
        else:
            count_dict[count] = 1
        if count >= 5:
            if count not in very_unique_peps:
                very_unique_peps[count] = []
            else:
                very_unique_peps[count].append(peptide)
    return count_dict, very_unique_peps

def find_neighbor(peptides, lib):
    '''
    This will do a local alignment with a given peptide to a library of your choice, returning
    the best sequence.
    '''
    for peptide in peptides:
        neighbors = []
        best_score = 20
        i = 0
        for pep in lib:
            i += 1
            if pep != peptide:
                if i%1000 == 0:
                    print i, best_score, neighbors
                try:
                    align = pairwise2.align.localxx(peptide, pep)
                    if align[0][2] >= best_score:
                        print align
                        neighbors.append(pep)
                except:
                    pass
    return neighbors, best_score

def read_file(in_file):
    '''
    This function will read in the input peptide file, trim to remove the constant region leaving
    only the c-term constant G if there.
    '''
    library = dict()
    for line in in_file:
        temp = line.strip()
        if len(temp) > 8:
            if temp.count('_') > 1:
                continue
            rand_region = temp[:(temp.find('GSGS',5)+1)]
            if rand_region == '':
                if temp[-12:-1] == 'HHHHHHRL_VA':
                    rand_region = temp[:-17]
            if rand_region not in library:
                    library[rand_region] = 1
            else:
                    library[rand_region] += 1
        else:
            pass
    print 'the size of the library is %s' %len(library)
    return library

def write_dict():
    '''
    Simple function for writing the sequence-count dictionary
    '''
    out_file = open('Sequence.csv', 'w')
    out_file.write('Sequence, Count \n')
    for pep in pep_dict:
        out_file.write(pep + ',' + str(pep_dict[pep]) + '\n')
    out_file.close()

def get_AA_freq(pep_dict):
    '''
    This will make a frequency dictionary of the library passed to this function
    '''
    freq_dict = dict()
    row_labels = list('ACDEFGHIKLMNPQRSTVWY')
    for AA in row_labels:
        freq_dict[AA] = np.zeros((1,32))
    try:
        for pep in pep_dict:
            for i, AA in enumerate(pep):
                if AA in freq_dict:
                    freq_dict[AA][0][i] += 1
    except:
        pass
    return freq_dict

def AA_heatmap(freq_dict, sample_name):
    '''
    This will plot an dictionary input where the rows are the AAs and the columns are the
    positions and each cell is the count or in this case, normalized frequency
    '''
    #change data from a dictionary to a sorted matrix for plotting with matplotlib
    data_matrix = np.zeros((1,32))
    for AA in sorted(freq_dict):
        data_matrix = np.vstack((data_matrix, freq_dict[AA]))
    data_matrix = data_matrix[1:, 1:-1]
    data_matrix = normalize(data_matrix,axis=0, norm='l1')

    #build the column and row labels
    column_labels = [str(x) for x in range(1,31)]
    row_labels = list('ACDEFGHIKLMNPQRSTVWY')

    #make the heatmap plot
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data_matrix, cmap=plt.cm.Greens)

    # Set fonts and titles after making a FontProperties object which allow you to change font
    font_path = 'C:\Windows\Fonts\\times.ttf'
    font_prop = font_manager.FontProperties(fname=font_path, size=50)

    plt.xlabel('Residue Position', fontproperties=font_prop, size=50)
    ax.xaxis.set_label_position('top')
    plt.ylabel('Amino Acid', fontproperties=font_prop, size=50)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    title = 'Amino Acid Frequencies for %s' %sample_name
    plt.title(title, y=1.08, fontproperties=font_prop, size=60)

    #Set the tick labels to be in the middle of the box and rename them
    ax.set_xticklabels([str(x-1) for x in [i for i in range(1,32)]], fontproperties=font_prop, size=35)
    ax.set_yticklabels(row_labels, fontproperties=font_prop, size=35)
    ax.xaxis.set_major_locator(ticker.FixedLocator([-0.5+i for i in range(33)]))
    ax.yaxis.set_major_locator(ticker.FixedLocator([0.5+i for i in range(20)]))

    #make the colorbar and name it
    cb = plt.colorbar(heatmap)
    cb.set_label('Amino Acid Frequency', fontproperties=font_prop, size=50)
    cb.ax.tick_params(labelsize=35)

    #plot it all!
    plt.show()

################# END FUNCTIONS ##################
def main():
    DTT_file = open('DTT_peps_final.txt')
    NoDTT_file = open('NoDTT_peps_final.txt')
    DTT_library = read_file(DTT_file)
    #lib_C, lib_no_C = cys_dict(DTT_library)
    NoDTT_library = read_file(NoDTT_file)
    DTT_file.close()
    NoDTT_file.close()
    return DTT_library, NoDTT_library

if __name__ == "__main__":
    DTT_library, NoDTT_library = main()



################### HEATMAPS ####################


# First create objects storing the names you want for the heatmaps

# sample_name = 'complete No DTT library'
# sample_name_C = 'No DTT library - Cysteine containing portion'
# sample_name_No_C = 'No DTT library - No cysteine containing portion'
# sample_name_1 = 'No DTT library, 2+ copy number'
# sample_name_2 = 'No DTT library, 3+ copy number'
# sample_name_3 = 'No DTT library, 4+ copy number'
# sample_name_4 = 'No DTT library, 5+ copy number'
# sample_name_5 = 'No DTT library, 6+ copy number'
# sample_name_6 = 'No DTT library, 7+ copy number'
# sample_name_30 = 'No DTT library, 31+ copy number'
# sample_name_1_noc = 'No DTT library, 2+ copy number, No Cys'
# sample_name_2_noc = 'No DTT library, 3+ copy number, No Cys'
# sample_name_3_noc = 'No DTT library, 4+ copy number, No Cys'
# sample_name_4_noc = 'No DTT library, 5+ copy number, No Cys'
# sample_name_5_noc = 'No DTT library, 6+ copy number, No Cys'
# sample_name_6_noc = 'No DTT library, 7+ copy number, No Cys'


#Next, use the get_AA_freq on the dictionary of interest, below I am showing it used on DTT_library, the cysteine
#containing sequences in the DTT_library and then on dictionaries made from sequences having copy numbers greater than
#a specific value

# matrix = get_AA_freq(DTT_library)
# matrix_C = get_AA_freq(lib_C)
# great_1 = find_values_greater_than(DTT_library,1)
# great_matrix_1 = get_AA_freq(great_1)
# great_2 = find_values_greater_than(DTT_library,2)
# great_matrix_2 = get_AA_freq(great_2)
# great_3 = find_values_greater_than(DTT_library,3)
# great_matrix_3 = get_AA_freq(great_3)
# great_4 = find_values_greater_than(DTT_library,4)
# great_matrix_4 = get_AA_freq(great_4)
# great_5 = find_values_greater_than(DTT_library,5)
# great_matrix_5 = get_AA_freq(great_5)
# great_6 = find_values_greater_than(DTT_library,6)
# great_matrix_6 = get_AA_freq(great_6)
# great_30 = find_values_greater_than(DTT_library,30)
# great_matrix_30 = get_AA_freq(great_30)



#Same deal as above just different libraries

# matrix_No_C = get_AA_freq(lib_no_C)
# great_1_noc = find_values_greater_than(lib_no_C,1)
# great_matrix_1_noc = get_AA_freq(great_1_noc)
# great_2_noc = find_values_greater_than(lib_no_C,2)
# great_matrix_2_noc = get_AA_freq(great_2_noc)
# great_3_noc = find_values_greater_than(lib_no_C,3)
# great_matrix_3_noc = get_AA_freq(great_3_noc)
# great_4_noc = find_values_greater_than(lib_no_C,4)
# great_matrix_4_noc = get_AA_freq(great_4_noc)
# great_5_noc = find_values_greater_than(lib_no_C,5)
# great_matrix_5_noc = get_AA_freq(great_5_noc)
# great_6_noc = find_values_greater_than(lib_no_C,6)
# great_matrix_6_noc = get_AA_freq(great_6_noc)


# Now pass the frequency dictionaries to the AA_heatmap function and out pops the heatmap!

# AA_heatmap(matrix, sample_name)
# AA_heatmap(matrix_C, sample_name_C)
# AA_heatmap(great_matrix_1, sample_name_1)
# AA_heatmap(great_matrix_2, sample_name_2)
# AA_heatmap(great_matrix_3, sample_name_3)
# AA_heatmap(great_matrix_4, sample_name_4)
# AA_heatmap(great_matrix_5, sample_name_5)
# AA_heatmap(great_matrix_6, sample_name_6)
# AA_heatmap(great_matrix_30, sample_name_30)
#
# AA_heatmap(matrix_No_C, sample_name_No_C)
# AA_heatmap(great_matrix_1_noc, sample_name_1_noc)
# AA_heatmap(great_matrix_2_noc, sample_name_2_noc)
# AA_heatmap(great_matrix_3_noc, sample_name_3_noc)
# AA_heatmap(great_matrix_4_noc, sample_name_4_noc)
# AA_heatmap(great_matrix_5_noc, sample_name_5_noc)
# AA_heatmap(great_matrix_6_noc, sample_name_6_noc)


