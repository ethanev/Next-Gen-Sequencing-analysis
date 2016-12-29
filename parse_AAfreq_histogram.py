#in_file = open('EE01-NoDTTlibR6-S1_R1-seq_forward_histo.txt')
in_file = open(input('please enter a histogram file: ' ))

histo_list = []

for line in in_file:
    delimiter = '('
    histo_list = line.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = ')'
    histo_list = histo_str.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = '{'
    histo_list = histo_str.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = '}'
    histo_list = histo_str.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = ':'
    histo_list = histo_str.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = "'"
    histo_list = histo_str.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = ','
    histo_list = histo_str.split(delimiter)
    delimiter = ''
    histo_str = delimiter.join(histo_list)
    delimiter = ' '
    histo_list = histo_str.split(delimiter)
    
AA_list = ['A', 'N', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
final_list = []
for i in range(0,len(histo_list)):
    if histo_list[i] in AA_list:
        temp1 = histo_list[i]
        temp2 = histo_list[i+1]
        temp3 = histo_list[i+2]
        print(temp1, temp2, temp3)
        #final_list.append(te)
    else:
        pass
delimiter = '\n'
out = delimiter.join(final_list)
print(out)    
#for each ', (' in in_file:
    
