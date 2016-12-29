import sys
import numpy as np
import scipy
import matplotlib as plt
import copy
from sklearn.feature_extraction import DictVectorizer
from sklearn.cluster import MiniBatchKMeans
import time
import re
from Bio import pairwise2
sys.setrecursionlimit(1000000)


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

def cleanlib(dictionary):
    pep_lib = dict()
    bad_seq_dict = dict()
    bad_seq = 0
    good_seq = 0
    for peptide in dictionary:
        if '_' in peptide[1:25] or len(peptide)<25:
            bad_seq += 1
            bad_seq_dict[peptide] = dictionary[peptide]
        else:
            good_seq += 1
            pep_lib[peptide] = dictionary[peptide]
    return pep_lib, good_seq, bad_seq, bad_seq_dict

def givencopynum_vs_cpnum(dictionary):
    '''

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

def cys_pattern(dictionary):
    '''
This is to analyze the AA composition surrounding the cysteines
    '''
    C_minus_2 = dict()
    C_minus_1 = dict()
    C_plus_1 = dict()
    C_plus_2 = dict()
    C_plus_3 = dict()
    C_minus_3 = dict()
    AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for peptide in dictionary:
        for i in range(0,len(peptide)):
            if peptide[i] == 'C':
                if peptide[i-3] in (C_minus_3):
                    C_minus_3[peptide[i-3]] +=1
                if peptide[i-3] not in (C_minus_3):
                    C_minus_3[peptide[i-3]] = 1
                if peptide[i-2] in (C_minus_2):
                    C_minus_2[peptide[i-2]] +=1
                if peptide[i-2] not in (C_minus_2):
                    C_minus_2[peptide[i-2]] = 1
                if peptide[i-1] in (C_minus_1):
                    C_minus_1[peptide[i-1]] +=1
                if peptide[i-1] not in (C_minus_1):
                    C_minus_1[peptide[i-1]] = 1
                if i != (len(peptide)-1) and i != (len(peptide)-2) and i != (len(peptide)-3):
                    if peptide[i+3] in (C_plus_3):
                        C_plus_3[peptide[i+3]] +=1
                    if peptide[i+3] not in (C_plus_3):
                        C_plus_3[peptide[i+3]] = 1
                if i != (len(peptide)-1) and i != (len(peptide)-2):
                    if peptide[i+2] in (C_plus_2):
                        C_plus_2[peptide[i+2]] +=1
                    if peptide[i+2] not in (C_plus_2):
                        C_plus_2[peptide[i+2]] = 1
                if i != (len(peptide)-1):
                    if peptide[i+1] in (C_plus_1):
                        C_plus_1[peptide[i+1]] +=1
                    if peptide[i+1] not in (C_plus_1):
                        C_plus_1[peptide[i+1]] = 1
                else:
                    pass
    return C_minus_3, C_minus_2, C_minus_1, C_plus_1, C_plus_2, C_plus_3

def new_cys_pattern(dictionary):
    near_cys_pattern = dict()
    error_pep = []
    for sequence in dictionary:
        for _ in re.finditer('C', sequence):
            try:
                pattern = sequence[_.start()-2:_.start()+4]
                if pattern in near_cys_pattern:
                    near_cys_pattern[pattern] += 1
                else:
                    near_cys_pattern[pattern] = 1
            except:
                error_pep.append(sequence)
    return near_cys_pattern, error_pep

def surrounding_cys(dictionary):
    near_cys = dict()
    for peptide in dictionary:
        for i in range(0,len(peptide)):
            if peptide[i] =='C':
                if i !=len(peptide[i-1]) and i !=len(peptide[i-2]) and i !=len(peptide[i-3]):
                    seq = peptide[(i-3):(i+4)]
                    if seq in near_cys:
                        near_cys[seq] += 1
                    else:
                        near_cys[seq] = 1
    return near_cys

def surrounding_cys_longer(dictionary):
    near_cys = dict()
    for peptide in dictionary:
        for i in range(0,len(peptide)):
            if peptide[i] =='C':
                if i !=len(peptide[i-1]) and i !=len(peptide[i-2]) and i !=len(peptide[i-3]) and i != len(peptide[i-4]):
                    seq = peptide[(i-4):(i+5)]
                    if seq in near_cys:
                        if len(seq) == 9:
                            near_cys[seq] += dictionary[peptide]
                        else:
                            pass
                    else:
                        if len(seq) == 9:
                            near_cys[seq] = dictionary[peptide]
                        else:
                            pass
    return near_cys

def one_mut_seq(dictionary):
    AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    check_seq = 'MHQKYKMTKDCFFSFLAHHKKRKLYPMSG'
    for i in range(0,len(check_seq)-1):
        for amino_acid in AA_list:
            trial_seq = check_seq[0:i]+amino_acid+check_seq[i+1:]
            if trial_seq in dictionary and trial_seq != check_seq:
                print(trial_seq, dictionary[trial_seq])

master_list = ['MHQKYKMTKDCFFSFLAHHKKRKLYPMSG']
def one_mut_seq_large(dictionary, peptide_list):
    '''
    takes a list of peptides you wish to find find other sequences which differ by one AA and builds a large list of all
    of these peptides and then their neighbors.
    '''
    AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for peptide in peptide_list:
        for i in range(0,len(peptide)):
            for amino_acid in AA_list:
                trial_seq = peptide[0:i]+amino_acid+peptide[i+1:]
                if trial_seq in dictionary and trial_seq not in peptide_list:
                    counter = 1
                    peptide_list.append(trial_seq)
                    if counter == 1:
                        one_mut_seq_large(dictionary, peptide_list)
                    return peptide_list

def whole_dictionary_tree(dictionary):
    '''
Take the starting library and breaks it into
smaller libraries with C,without C and bad sequences
takes about 49 min to run on DTT lib
    '''
    C_dict = dict()
    C_dict_trees = dict()
    No_C_dict = dict()
    No_C_dict_trees = dict()
    badseq_dict = dict()
    dict2 = dict()
    i = 0
    j = 0
    k = 0
    for peptide in dictionary:
        if peptide not in dict2:
            if len(peptide)<18 or len(peptide)>35 or "_" in peptide:
                pass
##                temp = []
##                temp.append(peptide)
##                tree = one_mut_seq_large(dictionary,temp)
##                if tree == None:
##                    badseq_dict[peptide] = dictionary[peptide]
##                    dict2[peptide] = 1
##                elif len(tree)>1:
##                    badseq_dict[i] = tree
##                    i += 1
##                    for peptide in tree:
##                        badseq_dict[peptide] = dictionary[peptide]
##                        dict2[peptide] = 1
            elif 'C' not in peptide:
                temp = []
                temp.append(peptide)
                tree = one_mut_seq_large(dictionary,temp)
                if tree == None:
                    No_C_dict[peptide] = dictionary[peptide]
                    dict2[peptide] = 1
                else:
                    No_C_dict_trees[j] = tree
                    j += 1
                    for peptide in tree:
                        No_C_dict[peptide] = dictionary[peptide]
                        dict2[peptide] = 1
            else:
                temp = []
                temp.append(peptide)
                tree = one_mut_seq_large(dictionary,temp)
                if tree == None:
                    C_dict[peptide] = dictionary[peptide]
                    dict2[peptide] = 1
                else:
                    C_dict_trees[k] = tree
                    k += 1
                    for peptide in tree:
                        C_dict[peptide] = dictionary[peptide]
                        dict2[peptide] = 1
        else:
            pass
    return C_dict, C_dict_trees, No_C_dict, No_C_dict_trees, badseq_dict, i, j, k

def one_mut_smaller_seq(dictionary, check_seq):
    AA_list = ['A', 'C', 'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for peptide in dictionary:
        if check_seq in peptide:
            print("exact match!:", peptide, trial_seq, dictionary[peptide])
        for i in range(0,len(check_seq)-1):
            for amino_acid in AA_list:
                trial_seq = check_seq[0:i]+amino_acid+check_seq[i+1:]
                if trial_seq in peptide and trial_seq != check_seq:
                    print(peptide, trial_seq, dictionary[peptide])

def small_seq(dictionary):
   # AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    small_seq = "FCPF"
    FCPF_counter = 0
    NO_FCPF = 0
    for peptide in dictionary:
        if small_seq in peptide:
            FCPF_counter += 1
        else:
            NO_FCPF +=1
    print('Seq continaining FCPF', FCPF_counter, '\n', 'seq without FCPF', NO_FCPF)

def levenshtein_test(seed, dictionary,value):
    """Intended to be a simple levenshtein sequence clustering, currently only looking at sequences that are the
    same length as the seed"""
    global family
    for goodpeptide in seed:
        #if goodpeptide in family:
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
##           # print('Seed is:', goodpeptide, '\n', 'cluster is:', seed_cluster)
##            if goodpeptide not in family:
##                family[goodpeptide] = dictionary[goodpeptide]
##                seed.remove(goodpeptide)
##                #print("family is: ", family, "seed is: ", seed)
##            else:
##               # print("GOODPEPTIDE P E P T I D E A L R E A D Y     I N    F A M I L Y")
##                seed.remove(goodpeptide)
##            for peptide in seed_cluster:
##                if peptide not in family and peptide not in seed:
##                    #print("adding ", peptide, "to seed")
##                    seed.append(peptide)
##                   # print(seed)
##                else:
##                    #print("cluster P E P T I D E A L R E A D Y     I N    F A M I L Y")
##                    pass
##                levenshtein_test(seed,dictionary,value)
       # print(family)

def rec_levn(seed,dictionary,value):
    """
    Intended to be a simple levenshtein sequence clustering, currently only looking at sequences that are the
    same length as the seed
    """
    seed_cluster = dict()
    cluster = [seed]
    while value >= 0:
        for peptide in dictionary:
            counter = 0
            if len(peptide) < len(seed):
                pass
            if len(peptide) == len(seed):
                for i in range(0,len(peptide)-1):
                    if peptide[i] == seed[i]:
                        pass
                    else:
                        counter += 1
                if counter <= value:
                    seed_cluster[peptide] = dictionary[peptide]
                    if peptide not in cluster:
                        cluster.append(peptide)
                        rec_levn(peptide,dictionary,value-1)
            if len(peptide) > len(seed):
                buffer = len(peptide) - len(seed)
                temp_seed = seed + buffer*"_"
                for i in range(0,len(peptide)-1):
                    if peptide[i] == temp_seed[i]:
                        pass
                    else:
                        counter += 1
                if counter <= value:
                    seed_cluster[peptide] = dictionary[peptide]
                    if peptide not in cluster:
                        cluster.append(peptide)
                        rec_levn(peptide,dictionary,value-1)
    return seed_cluster, cluster

def levn_cluster(cluster,dictionary,value):
    """
    Intended to be a simple levenshtein sequence clustering, currently only looking at sequences that are the
    same length as the seed
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
                        #print peptide
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
                        #print peptide
                        cluster.append(peptide[:30])
    return seed_cluster, cluster

def count_cys_peptides(library):
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

def five_and_sixmers():
    five = dict()
    six = dict()
    AAs= ['A','C','D','E','F','G','H','I','K','L','M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V','W','Y']
    for aa1 in AAs:
        for aa2 in AAs:
            for aa3 in AAs:
                for aa4 in AAs:
                    for aa5 in AAs:
                        fivemer = ''.join((aa1,aa2,aa3,aa4,aa5))
                        five[fivemer] = 0
                        for aa6 in AAs:
                            sixmer = ''.join((aa1,aa2,aa3,aa4,aa5,aa6))
                            six[sixmer] = 0
    return five, six

def vectorize_lib(library):
    vector_temp = vector_template()
    vector_lib = dict()
    matrix = np.ndarray(shape=(1,8420))
    i = 0
    for sequence in library:
        if i == 50:
            break
        seq_vector = vectorize(sequence, vector_temp)
        vec = DictVectorizer()
        vector_lib[sequence] = vec.fit_transform(seq_vector).toarray()
        v = vec.fit_transform(seq_vector).toarray()
        matrix = np.concatenate((matrix,v), axis=0)
        i += 1
    return vector_lib, matrix

def vectorize(sequence, vector_template):
    '''
    makes a 1-D vector of patterns present in the sequence
    '''
    vector = copy.deepcopy(vector_template)
    for AA in sequence:
        if AA == '_':
            continue
        vector[AA] += 1
    for i in xrange(len(sequence)-1):
        if sequence[i] == '_' or sequence[i+1] == '_':
            continue
        vector[sequence[i:i+2]] += 1
    for i in xrange(len(sequence)-2):
        if sequence[i] == '_' or sequence[i+1] == '_' or sequence[i+2] == '_':
            continue
        vector[sequence[i:i+3]] += 1
    return vector

def vector_template():
    vector_temp = dict()
    AAs = 'ACDEFGHIKLMNPQRSTVWY'
    for AA in AAs:
        vector_temp[AA] = 0
    for AA1 in AAs:
        for AA2 in AAs:
            vector_temp[AA1+AA2] = 0
    for AA1 in AAs:
        for AA2 in AAs:
            for AA3 in AAs:
                vector_temp[AA1+AA2+AA3] = 0
    return vector_temp

def vector_lookup(vector, vector_lib):
    pass

def cluster_sequences(clusters, vector_lib, vector_matrix):
    cluster_and_sequence = dict()
    for i in xrange(len(clusters)):
        seq = vector_lookup(vector_matrix[i], vector_lib)
        print seq
        cluster_and_sequence[clusters[i]].append(seq)
    return cluster_and_sequence

def reduce_seq_complexity(sequence):
    temp_seq = ''
    for aa in sequence:
        if aa == 'F' or aa == 'Y' or aa == 'W':
            temp_seq += 'A'
        elif aa == 'A' or aa == 'V' or aa == 'I' or aa =='L':
            temp_seq += 'Y'
        elif aa == 'D' or aa == 'E':
            temp_seq += '-'
        elif aa == 'K' or aa == 'R':
            temp_seq += '+'
        elif aa == 'N' or aa == 'Q':
            temp_seq += 'Q'
        else:
            temp_seq += aa
    return temp_seq

def read_file(in_file, out_file):
    counter = 0 # counter for longer than 8 AA and if there are more than 1 stop codons
    glob = 0 # count the number of sequences longer than 8 AAs
    library = dict()
    diff_lib = dict()
    for line in in_file:
        temp = line.strip()
        if len(temp) > 8:
            glob +=1
            if temp.count('_') > 1:
                counter += 1
                continue
            # rand_region = temp[:(temp.find('GSGS',22)+1)]
            rand_region = temp[:(temp.find('GSGS',5)+1)]
            if rand_region == '':
                if temp[-12:-1] == 'HHHHHHRL_VA':
                    rand_region = temp[:-17]
            if rand_region not in library:
                    library[rand_region] = 1
            else:
                    library[rand_region] += 1
            # if temp[-12:-1] == 'HHHHHHRL_VA':
            #     rand_region = temp[:-17]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-11:-1] == 'HHHHHHRL_V':
            #     rand_region = temp[:-16]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-5:-1] == 'GSGS':
            #     rand_region = temp[:-4]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-4:-1] == 'GSG':
            #     rand_region = temp[:-3]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] == 'SGSLGH':
            #     rand_region = temp[:-7]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-6:-1] == 'GSGSL':
            #     rand_region = temp[:-5]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='GSGSLG':
            #     rand_region = temp[:-6]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='HHHHHH':
            #     rand_region = temp[:-12]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='GHHHHH':
            #     rand_region = temp[:-11]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='LGHHHH':
            #     rand_region = temp[:-10]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='SLGHHH':
            #     rand_region = temp[:-9]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='GSLGHH':
            #     rand_region = temp[:-8]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='HHHRL_':
            #     rand_region = temp[:-15]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='HHRL_V':
            #     rand_region = temp[:-16]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-7:-1] =='HHHHRL':
            #     rand_region = temp[:-14]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-4:] =="HHHR":
            #     rand_region = temp[:-10]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-3:] == '_VA':
            #     rand_region = temp[:-16]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # elif temp[-2:] == 'SG':
            #     rand_region = temp[:-2]
            #     if rand_region not in library:
            #         library[rand_region] = 1
            #     else:
            #         library[rand_region] += 1
            # else:
            #     out_file.write(line)
            #     counter +=1
            #     if temp not in library:
            #         library[temp] = 1
            #     else:
            #         library[temp] += 1
            #     if temp not in diff_lib:
            #         diff_lib[temp] = 1
            #     else:
            #         diff_lib[temp] += 1
        else:
            pass
    print 'number of sequences greater than 8 and with more than 1 stop: %d' %counter
    print 'number of sequences greater than 8 AAs is: %d' %glob
    print 'the size of the library is %s' %len(library)
    return library, diff_lib

def Data_analysis(data_set, in_file, out_file):
    print '######## For', data_set, '###############\n\n'
    library, diff_lib = read_file(in_file, out_file)  #parse the in file and make master dictionaries
    in_file.close()
    out_file.close()

#     total_cp_vs_cp = givencopynum_vs_cpnum(library)
#     print('copy number statistics for total library:        ', total_cp_vs_cp, '\n')
#     greater_4 = find_values_greater_than(library,4)
#     print('values in library found more than 4 times: ','\n', greater_4, '\n')
#     print('length of total lib is: ', len(library))
#     pep_lib, good_seq, bad_seq, bad_seq_dict = cleanlib(library)
#     print('length of "trimmed" peptide library is: ', len(pep_lib))
#     print("good_seq= ", good_seq,'\n', "bad_seq= ", bad_seq, '\n')
#     trimmed_cp_vs_cp = givencopynum_vs_cpnum(pep_lib)
#     print('copy number statistics for trimmed library (after cleanlib fxn):        ', trimmed_cp_vs_cp, '\n')
#     greater_4_good = find_values_greater_than(pep_lib,4)
#     print('good peptide list peptides found more than 4 times:', '\n', greater_4_good, '\n')
#     greater_4_bad = find_values_greater_than(bad_seq_dict,4)
#     print('bad peptide list peptides found more than 4 times:', '\n', greater_4_bad, '\n')
#     cpnum = givencopynum_vs_cpnum(pep_lib)
#     print ('copy number distribution in good sequence library: ','\n',  cpnum, '\n')
#     cpnum_total = givencopynum_vs_cpnum(library)
#     print ('copy number distribution in whole library: ', '\n', cpnum_total, '\n')
#     cys_p = cys_pattern(pep_lib)
#     print("pattern around cysteines:",'\n', cys_p, '\n')
#     surr_cys = surrounding_cys_longer(pep_lib)
#     smallseq_greater_60 = find_values_greater_than(surr_cys,60)
#     print("small sequences present at greater than 60 copies:",'\n', smallseq_greater_60, '\n')
# ##
# ##    print('####### Levenshtein tests #######', '\n')
# ##    family = dict()
# ##    seed = ['MHQKYKMTKDCFFSFLAHHKKRKLYPMSG']
# ##    levn = levenshtein_test(seed,library,5)
# ##    print(levn)
#
#     C_dict, C_dict_trees, No_C_dict, No_C_dict_trees, badseq_dict, i, j, k = whole_dictionary_tree(pep_lib)
#     print('length of C lib is: ', len(C_dict))
#     print('length of C tree lib is: ', len(C_dict_trees))
#     print('length of No C lib is: ', len(No_C_dict))
#     print('length of No C trees lib is: ', len(No_C_dict_trees))
#     print('length of badseq dict is: ', len(badseq_dict))
#     NoC_cp_vs_cp = givencopynum_vs_cpnum(No_C_dict)
#     print('copy number statistics for No C library:        ', NoC_cp_vs_cp, '\n')
#     C_cp_vs_cp = givencopynum_vs_cpnum(C_dict)
#     print('copy number statistics for Cys containing library:        ', C_cp_vs_cp, '\n')
#
#     print('####### Levenshtein tests 2 #######', '\n')
#     for peptide in greater_4_good:
#         templist = [peptide]
#         levn = rec_levn(templist,library,5)
#         print (levn)

########################################################################
def main():
    in_file1 = open('NoDTT_peps_new.txt')
    out_file1 = open('NoDTTlibR6_canthandle_new_1.txt','w')
    data_set1 = 'NoDTT_new_1'
    ##in_file2 = open('EE02-DTTlibR6_S2_R1_pepsequences.txt')
    ##out_file2 = open('EE02-DTTlibR6_S2_R1_pepsequences_canthandle.txt','w')
    ##data_set2 = 'DTT_S2_R1'
    ##in_file3 = open('EE01-NoDTTlibR6_S1_R2_pepsequences.txt')
    ##out_file3 = open('EE01-NoDTTlibR6_S1_R2_pepsequences_canthandle.txt','w')
    ##data_set3 = 'NoDTT_S1_R2'
    ##in_file4 = open('EE01-NoDTTlibR6_S1_R1_pepsequences.txt')
    ##out_file4 = open('EE01-NoDTTlibR6_S1_R1_pepsequences_canthandle.txt','w')
    ##data_set4 = 'NoDTT_S1_R1'

    Data_analysis(data_set1, in_file1,out_file1)
    ##Data_analysis(data_set2, in_file2,out_file2)
    ##Data_analysis(data_set3, in_file3,out_file3)
    ##Data_analysis(data_set4, in_file4,out_file4)

# if __name__ == "__main__":
#     main()

# in_file = open('DTT_peps_final.txt')
# out_file = open('DTTlibR6_canthandle_new_2.txt','w')
in_file = open('EE02-DTTlibR6_S2_R1_pepsequences.txt')
out_file = open('DTTlibR6_canthandle_old.txt','w')

library, diff_lib = read_file(in_file, out_file)
in_file.close()
out_file.close()
###### Vectorize the library ######
#lib_vector, lib_matrix = vectorize_lib(library)
###### perform minibatch Kmean ######
# mbk = MiniBatchKMeans(init='k-means++', n_clusters=10, batch_size=10, n_init=10,
#                       max_no_improvement=None, verbose=True)
# t0 = time.time()
# mbk.fit(lib_matrix)
# centers = mbk.cluster_centers_
# predicted_cluster = mbk.fit_predict(lib_matrix)
# cluster_and_sequences = cluster_sequences(predicted_cluster, lib_vector, lib_matrix)
#
# t1 = time.time()

#######
#  BELOW IS AN ALT METHOD FOR DATA ANALYSIS AND READING IN
#######

# def Data_analysis1(data_set, in_file, out_file):
#     print('######## For', data_set, '###############',  "\n", '\n' )
#     library, diff_lib = read_file(in_file, out_file)  #parse the in file and make master dictionaries
#     in_file.close()
#     out_file.close()
#     total_cp_vs_cp = givencopynum_vs_cpnum(library)
#     print('copy number statistics for total library:        ', total_cp_vs_cp, '\n')
#     greater_4 = find_values_greater_than(library,4)
#     print('values in library found more than 4 times: ','\n', greater_4, '\n')
#     print('length of total lib is: ', len(library))
#     pep_lib, good_seq, bad_seq, bad_seq_dict = cleanlib(library)
#     print('length of "trimmed" peptide library is: ', len(pep_lib))
#     print("good_seq= ", good_seq,'\n', "bad_seq= ", bad_seq, '\n')
#     trimmed_cp_vs_cp = givencopynum_vs_cpnum(pep_lib)
#     print('copy number statistics for trimmed library (after cleanlib fxn):        ', trimmed_cp_vs_cp, '\n')
#     greater_4_good = find_values_greater_than(pep_lib,4)
#     print('good peptide list peptides found more than 4 times:', '\n', greater_4_good, '\n')
#     greater_4_bad = find_values_greater_than(bad_seq_dict,4)
#     print('bad peptide list peptides found more than 4 times:', '\n', greater_4_bad, '\n')
#     cpnum = givencopynum_vs_cpnum(pep_lib)
#     print ('copy number distribution in good sequence library: ','\n',  cpnum, '\n')
#     cpnum_total = givencopynum_vs_cpnum(library)
#     print ('copy number distribution in whole library: ', '\n', cpnum_total, '\n')
#     cys_p = cys_pattern(pep_lib)
#     print("pattern around cysteines:",'\n', cys_p, '\n')
#     surr_cys = surrounding_cys_longer(pep_lib)
#     smallseq_greater_60 = find_values_greater_than(surr_cys,60)
#     print("small sequences present at greater than 60 copies:",'\n', smallseq_greater_60, '\n')
# ##
# ##    print('####### Levenshtein tests #######', '\n')
# ##    family = dict()
# ##    seed = ['MHQKYKMTKDCFFSFLAHHKKRKLYPMSG']
# ##    levn = levenshtein_test(seed,library,5)
# ##    print(levn)
#
#     C_dict, C_dict_trees, No_C_dict, No_C_dict_trees, badseq_dict, i, j, k = whole_dictionary_tree(pep_lib)
#     print('length of C lib is: ', len(C_dict))
#     print('length of C tree lib is: ', len(C_dict_trees))
#     print('length of No C lib is: ', len(No_C_dict))
#     print('length of No C trees lib is: ', len(No_C_dict_trees))
#     print('length of badseq dict is: ', len(badseq_dict))
#     NoC_cp_vs_cp = givencopynum_vs_cpnum(No_C_dict)
#     print('copy number statistics for No C library:        ', NoC_cp_vs_cp, '\n')
#     C_cp_vs_cp = givencopynum_vs_cpnum(C_dict)
#     print('copy number statistics for Cys containing library:        ', C_cp_vs_cp, '\n')
#
#     print('####### Levenshtein tests 2 #######', '\n')
#     for peptide in greater_4_good:
#         templist = [peptide]
#         levn = rec_levn(templist,library,5)
#         print (levn)
#
# ########################################################################
#
# #in_file = open(input("Peptide sequence file: "))
# #out_file = open(input("File for peptides script can't handle: "), 'w')
#
#
# in_file1 = open('DTT_peps_new.txt')
# out_file1 = open('DTTlibR6_canthandle_new_1.txt','w')
# data_set1 = 'DTT_new_1'
# ##in_file2 = open('EE02-DTTlibR6_S2_R1_pepsequences.txt')
# ##out_file2 = open('EE02-DTTlibR6_S2_R1_pepsequences_canthandle.txt','w')
# ##data_set2 = 'DTT_S2_R1'
# ##in_file3 = open('EE01-NoDTTlibR6_S1_R2_pepsequences.txt')
# ##out_file3 = open('EE01-NoDTTlibR6_S1_R2_pepsequences_canthandle.txt','w')
# ##data_set3 = 'NoDTT_S1_R2'
# ##in_file4 = open('EE01-NoDTTlibR6_S1_R1_pepsequences.txt')
# ##out_file4 = open('EE01-NoDTTlibR6_S1_R1_pepsequences_canthandle.txt','w')
# ##data_set4 = 'NoDTT_S1_R1'
#
#
# Data_analysis(data_set1, in_file1,out_file1)
# ##Data_analysis(data_set2, in_file2,out_file2)
# ##Data_analysis(data_set3, in_file3,out_file3)
# ##Data_analysis(data_set4, in_file4,out_file4)
#







####
####
# THE FOLLOWING IS YET ANOTHER ALTERNATIVE ENDING:
# This can be used to combine two data files to analyze both together
###
####
# def read_file(in_file, out_file, library, diff_lib):
#     counter = 0
#     glob = 0
#     #library = dict()
#     #diff_lib = dict()
#     for line in in_file:
#         temp = line.strip()
#         if len(temp) > 20:
#             glob +=1
#             #if "C" in line:
#             if temp[-12:-1] == 'HHHHHHRL_VA':
#                 rand_region = temp[:-17]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-11:-1] == 'HHHHHHRL_V':
#                 rand_region = temp[:-16]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-5:-1] == 'GSGS':
#                 rand_region = temp[:-4]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-4:-1] == 'GSG':
#                 rand_region = temp[:-3]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] == 'SGSLGH':
#                 rand_region = temp[:-7]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-6:-1] == 'GSGSL':
#                 rand_region = temp[:-5]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='GSGSLG':
#                 rand_region = temp[:-6]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='HHHHHH':
#                 rand_region = temp[:-12]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='GHHHHH':
#                 rand_region = temp[:-11]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='LGHHHH':
#                 rand_region = temp[:-10]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='SLGHHH':
#                 rand_region = temp[:-9]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='GSLGHH':
#                 rand_region = temp[:-8]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='HHHRL_':
#                 rand_region = temp[:-15]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='HHRL_V':
#                 rand_region = temp[:-16]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-7:-1] =='HHHHRL':
#                 rand_region = temp[:-14]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             elif temp[-4:] =="HHHR":
#                 rand_region = temp[:-10]
#                 if rand_region not in library:
#                     library[rand_region] = 1
#                 else:
#                     library[rand_region] += 1
#             else:
#                 out_file.write(line)
#                 counter +=1
#                 if temp not in library:
#                     library[temp] = 1
#                 else:
#                     library[temp] += 1
#                 if temp not in diff_lib:
#                     diff_lib[temp] = 1
#                 else:
#                     diff_lib[temp] += 1
#         else:
#             pass
#     print(counter, glob)
#     return library, diff_lib
#
# def Data_analysis2(data_set, in_file, out_file):
#     print('######## For', data_set, '###############',  "\n", '\n' )
#     library, diff_lib = read_file(in_file, out_file)  #parse the in file and make master dictionaries
#     in_file.close()
#     out_file.close()
#     total_cp_vs_cp = givencopynum_vs_cpnum(library)
#     print('copy number statistics for total library:        ', total_cp_vs_cp, '\n')
#     greater_4 = find_values_greater_than(library,4)
#     print('values in library found more than 4 times: ','\n', greater_4, '\n')
#     print('length of total lib is: ', len(library))
#     pep_lib, good_seq, bad_seq, bad_seq_dict = cleanlib(library)
#     print('length of "trimmed" peptide library is: ', len(pep_lib))
#     print("good_seq= ", good_seq,'\n', "bad_seq= ", bad_seq, '\n')
#     trimmed_cp_vs_cp = givencopynum_vs_cpnum(pep_lib)
#     print('copy number statistics for trimmed library (after cleanlib fxn):        ', trimmed_cp_vs_cp, '\n')
#     greater_4_good = find_values_greater_than(pep_lib,4)
#     print('good peptide list peptides found more than 4 times:', '\n', greater_4_good, '\n')
#     greater_4_bad = find_values_greater_than(bad_seq_dict,4)
#     print('bad peptide list peptides found more than 4 times:', '\n', greater_4_bad, '\n')
#     cpnum = givencopynum_vs_cpnum(pep_lib)
#     print ('copy number distribution in good sequence library: ','\n',  cpnum, '\n')
#     cpnum_total = givencopynum_vs_cpnum(library)
#     print ('copy number distribution in whole library: ', '\n', cpnum_total, '\n')
#     cys_p = cys_pattern(pep_lib)
#     print("pattern around cysteines:",'\n', cys_p, '\n')
#     surr_cys = surrounding_cys_longer(pep_lib)
#     smallseq_greater_60 = find_values_greater_than(surr_cys,60)
#     print("small sequences present at greater than 60 copies:",'\n', smallseq_greater_60, '\n')
# ##
# ##    print('####### Levenshtein tests #######', '\n')
# ##    family = dict()
# ##    seed = ['MHQKYKMTKDCFFSFLAHHKKRKLYPMSG']
# ##    levn = levenshtein_test(seed,library,5)
# ##    print(levn)
#
#     C_dict, C_dict_trees, No_C_dict, No_C_dict_trees, badseq_dict, i, j, k = whole_dictionary_tree(pep_lib)
#     print('length of C lib is: ', len(C_dict))
#     print('length of C tree lib is: ', len(C_dict_trees))
#     print('length of No C lib is: ', len(No_C_dict))
#     print('length of No C trees lib is: ', len(No_C_dict_trees))
#     print('length of badseq dict is: ', len(badseq_dict))
#     NoC_cp_vs_cp = givencopynum_vs_cpnum(No_C_dict)
#     print('copy number statistics for No C library:        ', NoC_cp_vs_cp, '\n')
#     C_cp_vs_cp = givencopynum_vs_cpnum(C_dict)
#     print('copy number statistics for Cys containing library:        ', C_cp_vs_cp, '\n')
#
#     print('####### Levenshtein tests 2 #######', '\n')
#     for peptide in greater_4_good:
#         templist = [peptide]
#         levn = levenshtein_test(templist,library,5)
#         print (levn)
#
# ########################################################################
#
# #in_file = open(input("Peptide sequence file: "))
# #out_file = open(input("File for peptides script can't handle: "), 'w')
#
#
# #in_file1 = open('EE02-DTTlibR6_S2_R2_pepsequences.txt')
# #out_file1 = open('EE02-DTTlibR6_S2_R2_pepsequences_canthandle.txt','w')
# #data_set1 = 'DTT_S2_R2'
# in_file2 = open('EE02-DTTlibR6_S2_R1_pepsequences.txt')
# out_file2 = open('EE02-DTTlibR6_S2_R1_pepsequences_canthandle.txt','w')
# ##data_set2 = 'DTT_S2_R1'
# in_file3 = open('EE01-NoDTTlibR6_S1_R2_pepsequences.txt')
# out_file3 = open('EE01-NoDTTlibR6_S1_R2_pepsequences_canthandle.txt','w')
# ##data_set3 = 'NoDTT_S1_R2'
# ##in_file4 = open('EE01-NoDTTlibR6_S1_R1_pepsequences.txt')
# ##out_file4 = open('EE01-NoDTTlibR6_S1_R1_pepsequences_canthandle.txt','w')
# ##data_set4 = 'NoDTT_S1_R1'
#
# DTT = dict()
# diff_DTT = dict()
# NoDTT = dict()
# diff_NoDTT = dict()
#
# read_file(in_file2, out_file2, DTT, diff_DTT)
# read_file(in_file3, out_file3, NoDTT, diff_NoDTT)
#
# in_file2.close()
# out_file2.close()
# in_file3.close()
# out_file3.close()
#
#
#
# #Data_analysis(data_set1, in_file1,out_file1)
# ##Data_analysis(data_set2, in_file2,out_file2)
# ##Data_analysis(data_set3, in_file3,out_file3)
# ##Data_analysis(data_set4, in_file4,out_file4)


