# Author: Ethan D Evans

#Hi! below is the script I used to read a pair of FASTQ files and obtain peptide sequences to be used for
# further data analysis with the other script.

from Bio import pairwise2

def phred_convert(sequence):
    '''
    This will convert a phred value into its corresponding numerical score
    '''
    new_seq = [0]*len(sequence)
    for i in xrange(len(sequence)):
        new_seq[i] = str(ord(sequence[i])-33)
    return ' '.join(new_seq)

def comb_phred(combined_sequences, useavg=False, avgscore=30, usehardcutoff=False, hardcutoff=15, usepercent=False, percentabove=85, percentscore=30):
    '''
    Splits based on whether or not you are using paired end analysis. either way, it will return True (meaning the sequence passes
    the minimum phred score test or False in which case it is not further analyzed
    '''
    if len(combined_sequences) == 4:
        return accept_or_reject(combined_sequences[1], useavg, avgscore, usehardcutoff, hardcutoff, usepercent,
                                percentabove, percentscore) \
               or accept_or_reject(combined_sequences[3],useavg, avgscore, usehardcutoff, hardcutoff, usepercent,
                                   percentabove, percentscore)
    elif len(combined_sequences) == 2:
        return accept_or_reject(combined_sequences[1])
    else: raise ValueError

def accept_or_reject(sequence, useavg, avgscore, usehardcutoff, hardcutoff, usepercent, percentabove, percentscore):
    '''
    This is the phred score analyzer, it looks to see that the sequences have a high quality assignment
    '''
    phred_histo = {}
    #split the sequence string up and make a histogram of scores
    sequence_list = sequence.split(' ')
    for element in sequence_list:
        if element not in phred_histo:
            phred_histo[element] = 1
        else:
            phred_histo[element] += 1
    # if we want to keep sequences only with average scores greater than a certain values we use this code
    if useavg:
        seq_sum = 0
        for score in phred_histo:
            seq_sum += int(score)*phred_histo[score]
        seq_avg = float(seq_sum) / float(len(sequence_list))
        if seq_avg <= avgscore: return False
    # if we want to keep sequences without scores less than a certain amount we use the following
    if usehardcutoff:
        for i in xrange(hardcutoff):
            try:
                if bool(phred_histo[i]):
                    return False
            except:
                continue
    # this is the primary method for removing certain sequences
    if usepercent:
        items_above = 0
        for item in phred_histo.keys():
            if int(item) > percentscore:
                items_above += phred_histo[item]
        true_percent_above = float(items_above)/float(len(sequence_list)) * 100
        if true_percent_above <= float(percentabove):
            return False
    return True

def line_selector(line):
    '''
This function will determine if a sequence is forward (the first two if statements) or needs to be flipped
and reverse complemented
    '''
    if 'TTACAATG' in line:
        return translate_dna(line)
    if line[0:5] == 'TAATA':
        return translate_dna(line)
    if line[0:5] == 'CTAGC':
        new_seq = flip_seq(line)
        return translate_dna(new_seq)
    if "CCGGAGCC" in line:
        new_seq = flip_seq(line)
        return translate_dna(new_seq)
    if "CTATAGCC" in line:
        new_seq = flip_seq(line)
        return translate_dna(new_seq)
    else:
        pass

def align(combined_sequences):
    '''
    This function will align the two paired sequences and return an alignment for them
    '''
    flipped = flip_seq(combined_sequences[2])
    aligments = pairwise2.align.localxx(combined_sequences[0], flipped)
    return aligments[0]

def Align_and_pick(combined_sequences):
    '''
    An alternative aligning function for two paired sequences
    '''
    flipped = flip_seq(combined_sequences[2])
    phred_0 = combined_sequences[1].split(' ')
    phred_1 = combined_sequences[3].split(' ')
    alignments = pairwise2.align.localxs(combined_sequences[0], flipped, -0.5, -0.5)
    final_sequence = [_ for _ in range(len(alignments[0][0]))]
    base_count_0 = 0 #used to determine the position in the original sequence to find Q values
    base_count_1 = 0
    for i in xrange(len(alignments[0][0])):
        if alignments[0][0][i] == alignments[0][1][i]:
            final_sequence[i] = alignments[0][0][i]
            base_count_0 += 1
            base_count_1 += 1
        else:
            if alignments[0][0][i] == '-':
                final_sequence[i] = alignments[0][1][i]
                base_count_1 += 1
            elif alignments[0][1][i] == '-':
                final_sequence[i] = alignments[0][0][i]
                base_count_0 += 1
            else:
                # take the base with the higher q score
                if int(phred_0[base_count_0]) > int(phred_1[len(combined_sequences[2])-1-base_count_1]):
                    final_sequence[i] = alignments[0][0][i]
                    #print alignments[0][0][i], phred_0[base_count_0]
                    base_count_0 += 1
                    base_count_1 += 1
                else:
                    final_sequence[i] = alignments[0][1][i]
                    #print alignments[0][1][i], phred_1[len(combined_sequences[2])-1-base_count_1]
                    base_count_1 += 1
                    base_count_0 += 1
    return ''.join(final_sequence)


def switch(letter):
    '''
    Function used to help determine the reverse complement of a sequence
    '''
    if letter == 'A':
        return 'T'
    if letter == 'T':
        return 'A'
    if letter == 'G':
        return 'C'
    if letter == 'C':
        return 'G'
def complement_base(string):
    '''
Reads through the string and adds the complement base of the new_list then combines new_list into a string and return this
    '''
    new_list = []
    for letter in string:
        new_list.append(switch(letter))
    break_point = ''
    seq = break_point.join(new_list)
    return seq
def flip_seq(string):
    '''
Returns the reverse of the string and then returns the result of the complement_base function on the reversed sequence
    '''
    temp = string[::-1]
    return complement_base(temp)
def translate_dna(sequence):
    '''
    As the name implies, this will translate a seqeunce into its corresponding peptide
    '''
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    proteinsequence = ''
    start = sequence.find('ATG')
    sequencestart = sequence[int(start):]
    for n in range(0,len(sequencestart),3):
        if sequencestart[n:n+3] in codontable:
            proteinsequence += codontable[sequencestart[n:n+3]]
    return proteinsequence

def main():
    '''
    Reads in from initial FASTQ files, combines paired reads, removes low quality reads
    '''
    #in_file_1 = open(raw_input("please enter a FASTQ file name to parse: "))
    in_file_1 = open('EE01-NoDTTLibR6_S1_L001_R1_001.fastq')
    paired = 'yes'
    in_file_2 = open('EE01-NoDTTLibR6_S1_L001_R2_001.fastq')
    not_valid = True
    while not_valid: #this part would just be used in a more interactive environment with someone not familiar with programming - it wasn't used
        try:
            #paired = raw_input("Paired end? (please type 'yes' or 'no', case sensitive): ")
            if paired == 'yes':
                #in_file_2 = open(raw_input("please enter a second FASTQ file name to parse: "))
                not_valid = False
            elif paired != 'no':
                print 'You did not input a valid answer'
            else:
                paired = False
                not_valid = False
        except:
            print 'did you forget to use a quotation marks or input a nonexistant file?'
    #out_file = open(raw_input("please enter a .txt file name to write sequences to: "),'w')
    out_file = open('NoDTT_peps_final.txt','w')
    counter_seq = 0
    accept = 0
    reject = 0
    while in_file_1: # reads in the FASTQ and removes unneeded lines and gets the total phred score for the sequence
        in_file_1.readline()
        line2_f1 = in_file_1.readline().strip('\n')
        counter_seq += 1
        in_file_1.readline()
        line4_f1 = phred_convert(in_file_1.readline().strip('\n'))
        if bool(paired): # read in the sequence pair if this value is True and then combine with the first sequence
            in_file_2.readline()
            line2_f2 = in_file_2.readline().strip('\n')
            in_file_2.readline()
            line4_f2 = phred_convert(in_file_2.readline().strip('\n'))
            combined = [line2_f1, line4_f1, line2_f2, line4_f2]
        else:
            combined = [line2_f1, line4_f1]
        # now determine what sequences will actually be used. useavg will turn on flag for using a avg phred score and
        # the hard cutoff is is you want all positions to be above a certain value
        # usepercent is probably the best choice to use with 85% above Q30 ... this seems to be pretty standard
        write_seq = comb_phred(combined, useavg=False, avgscore=30, usehardcutoff=False, hardcutoff=15,
                             usepercent=True, percentabove=85, percentscore=30)
        if write_seq:
            # *****OPTION 1****** align all sequences and write -- takes many days
            #sequence = Align_and_pick(combined)
            # out_file.write(str(sequence))

            accept += 1 # just used to monitor the progress since depending on which option you choose, it can take a while!
            if accept % 1000 == 0:
                print accept

            #   *****OPTION 2**** just write all the DNA sequences
            # out_file.write(combined[0])
            # out_file.write('\n')
            # out_file.write(combined[1])
            # out_file.write('\n')
            # out_file.write(combined[2])
            # out_file.write('\n')
            # out_file.write(combined[3])
            # out_file.write('\n')

            # *****OPTION 3****** translate both pairs and if not the same then align, pick seq and translate
            seq1 = line_selector(combined[0])
            seq2 = line_selector(combined[2])
            try:
                if seq1[:30] != seq2[:30]:
                    # align = pairwise2.align.localxs(seq1,seq2,-.5,-.5)
                    # print align[0]
                    sequence = Align_and_pick(combined)
                    sequence = line_selector(sequence)
                    if sequence == None:
                        if 'GSGSLG' in seq1:
                            sequence = seq1
                        if 'GSGSLG' in seq2:
                            sequence = seq2
                    out_file.write(sequence)
                    out_file.write('\n')
                else:
                    if len(seq1) > len(seq2):
                        out_file.write(seq1)
                        out_file.write('\n')
                    else:
                        out_file.write(seq2)
                        out_file.write('\n')
            except:
                pass
        else:
            reject += 1
    print accept
    print reject
    print('total sequences: ', counter_seq)

    in_file_1.close()
    in_file_2.close()
    out_file.close()



if __name__ == '__main__':
    main()
