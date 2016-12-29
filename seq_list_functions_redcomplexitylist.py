in_file = open('EE02-DTTlibR6_S2_L001_R1_001-sequences_forward-2.txt')
out_file2 = open('EE02-DTTlibR6_S2_L001_R1_001-sequences_forward-2-canthandle.txt','w')
#out_file = open('EE02-DTT-R1-histogram.txt')
library = dict()
diff_lib = dict()

def find_key(input_dict, value):
    return {k for k, v in input_dict.items() if v == value}
    
def find_values_greater_than(dictionary,value):
    new_dict = dict()
    for k in dictionary:
        dict_val = dictionary[k]
        if dict_val > value:
            new_dict[k] = dict_val
    return new_dict

def analyze_seq(sequence):
    """
This is used to reduce the complexity of the AA space by group common AAs
    """
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
        elif aa == 'N' or aa == 'Q' or aa == 'S' or aa == 'T':
            temp_seq += 'Q'
        else:
            temp_seq += aa
    return temp_seq
            
def AA_histo_maker(dictionary):
    """
This make a histogram out of all AA frequencies at each position in the library
    """
    AA_histo = dict()
    protein_list = []
    for key in dictionary:
        protein_list.append(key)
    for sequence in protein_list:
        for i in range(0,len(sequence)-1):
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
    AA_list = ['Y', 'C', '-', 'A', 'H', '+', 'M', 'Q', 'P', 'G']
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
    print(C_minus_3, C_minus_2, C_minus_1, C_plus_1, C_plus_2, C_plus_3)

def surrounding_cys(dictionary):
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


def one_mut_seq(dictionary, seq):
    AA_list = ['Y', 'C', '-', 'A', 'H', '+', 'M', 'Q', 'P', 'G']
    check_seq = 'MHQ+A+MQ+-CAAQAYYHH++++YAPMQG'
    for i in range(0,len(seq)-1):
        for amino_acid in AA_list:
            trial_seq = seq[0:i]+amino_acid+seq[i+1:]
            if trial_seq in dictionary and trial_seq != seq:
                print(trial_seq, dictionary[trial_seq])

master_list = ['MHQ+A+MQ+-CAAQAYYHH++++YAPMQG']
def one_mut_seq_large(dictionary, peptide_list2):
    AA_list = ['Y', 'C', '-', 'A', 'H', '+', 'M', 'Q', 'P', 'G']
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

                    
def one_mut_smaller_seq(dictionary):
    AA_list = ['Y', 'C', '-', 'A', 'H', '+', 'M', 'Q', 'P', 'G']
    check_seq = 'MQ+-CAAQA'
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
    small_seq = "ACPA"
    FCPF_counter = 0
    NO_FCPF = 0
    for peptide in dictionary:
        if small_seq in peptide:
            FCPF_counter += 1
        else:
            NO_FCPF +=1
    print('Seq continaining ACPA', FCPF_counter, '\n', 'seq without ACPA', NO_FCPF)

def levenshtein_test(seed, dictionary,value):
    """Intended to be a simple levenshtein sequence clustering, currently only looking at sequences that are the 
    same length as the seed"""
    seed_cluster = dict()
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
    return seed_cluster




counter = 0
glob = 0
for line in in_file:
    temp = line.strip()
    if len(temp) > 20:
        glob +=1
        #if "C" in line:
        if temp[-12:-1] == 'HHHHHHRL_VA':
            seq = temp[:-17]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-11:-1] == 'HHHHHHRL_V':
            seq = temp[:-16]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-5:-1] == 'GSGS':
            seq = temp[:-4]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-4:-1] == 'GSG':
            seq = temp[:-3]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1                
        elif temp[-7:-1] == 'SGSLGH':
            seq = temp[:-7]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-6:-1] == 'GSGSL':
            seq = temp[:-5]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='GSGSLG':
            seq = temp[:-6]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='HHHHHH':
            seq = temp[:-12]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='GHHHHH':
            seq = temp[:-11]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='LGHHHH':
            seq = temp[:-10]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='SLGHHH':
            seq = temp[:-9]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='GSLGHH':
            seq = temp[:-8]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='HHHRL_':
            seq = temp[:-15]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='HHRL_V':
            seq = temp[:-16]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-7:-1] =='HHHHRL':
            seq = temp[:-14]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        elif temp[-4:] =="HHHR":
            seq = temp[:-10]
            rand_region = analyze_seq(seq)
            if rand_region not in library:
                library[rand_region] = 1
            else:
                library[rand_region] += 1
        else:
            out_file2.write(line)
            counter +=1
            seq = analyze_seq(temp)
            if seq not in library:
                library[seq] = 1
            else:
                library[seq] += 1
            if seq not in diff_lib:
                diff_lib[seq] = 1
            else:
                diff_lib[seq] += 1
    else:
        pass
print(counter, glob)
                
def line_print():
    in_file = open('EE02-DTTlibR6_S2_L001_R1_001-sequences_forward-2.txt')
    for line in in_file:
        if line[-12:-1] == 'HHHHHHRL_VA':
            rand_region = line[:-17]
        if line[-5:-1] == 'GSGS':
            ran_region = line[:-4]

    in_file.close()

   
alt = {'MQNSSGSCRTTRRTSARSCRCTLIATYKITGSGSLG': 1}
##histo = AA_histo_maker(library)
##for key in histo:
##    for a, b in key:
##        print("The count of the amino acid: ", a, "at the", b, 'position is', histo[(a,b)]) 

##near_cys = dict()
##for line in in_file:
        #seq = line.strip()
##    for i in range(0,len(seq)):
##            if seq[i] =='C':
##                if i !=len(seq[i-1]) and i !=len(seq[i-2]) and i !=seq(line[i-3]):
##                    realseq = seq[(i-4):(i+5)]
##                    if realseq in near_cys:
##                        near_cys[realseq] += 1
##                    else:
##                        near_cys[realseq] = 1
seq = input("please enter sequence you want to find neighbor sequences for:" )

in_file.close()
#out_file.close() 
out_file2.close()
        
        
