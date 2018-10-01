import sys
import re
import pandas as pd
import numpy as np
import linecache

from coding import coding
from coding import decoding

def main():
    args = sys.argv
    ## args[1] : threshold = 400?
    threshold = int(args[1])
    cancer_type_list = ['adrenal_gland', 'biliary_tract', 'bone', 'breast', 'central_nervous_system', 'cervix', 'endometrium',\
            'haematopoietic_and_lymphoid_tissue', 'kidney', 'large_intestine', 'liver', 'lung', 'oesophagus', 'ovary', 'pancreas', 'parathyroid',\
            'prostate', 'skin', 'small_intestine', 'soft_tissue', 'stomach', 'thyroid', 'upper_aerodigestive_tract',\
            'urinary_tract', 'vulva']

    Mutation_file = '~/Data/Mutation_Signature/CosmicMutantExport.tsv'
    pre_data = pd.read_csv(Mutation_file, delimiter = '\t')
    pre_data = select_data(pre_data)

    for cancer_type in cancer_type_list:
        data = select_cancer_type(pre_data, threshold, cancer_type)
        data.reset_index(drop=True, inplace=True)
        name_list = data['Sample name']
        last_document = 'xxxxxxxx'
        number_of_document = 0
        for name in name_list:
            if(name != str(last_document)):
                last_document = name
                number_of_document += 1
        if (number_of_document >= 10):
            output_to_file(data, cancer_type, number_of_document)

def select_data(data):
    ## Extract substitution
    data = data[data['Mutation Description'].str.contains('Substitution')]

    ## Select which use GRCh38
    data = data[data['GRCh'] == 38]
    data = data.loc[:,['Sample name', 'Mutation CDS', 'Mutation genome position', 'Mutation Description', 'Primary site', 'Primary histology', 'Mutation strand']]
    data = data.sort_values(by = 'Sample name')
    data.reset_index(drop = True, inplace = True)

    ## Select Mutation which is SNV
    count = 0
    CDS = data['Mutation CDS']
    drop_list = list()
    for x in CDS:
        char_list = list(x)
        if((char_list[len(char_list)-4] in ['A','C','G','T','>']) or char_list[len(char_list)-2] != '>'):
            drop_list.append(count)
        count += 1
    data.drop(drop_list, inplace = True)
    return data

def select_cancer_type(pre_data, threshold, cancer_type):
    if(cancer_type != 'all'):
        pre_data = pre_data.loc[pre_data['Primary site'] == cancer_type]
    pre_data.sort_values(by='Sample name', inplace=True)
    pre_data.reset_index(drop=True, inplace=True)
    name_list = pre_data['Sample name']
    last_document = name_list[0]
    drop_list = list()
    sum_of_words = 0
    temp_index = 0
    index_list = list()
    for name in name_list:
        if(name != last_document):
            if(sum_of_words < threshold):
                drop_list.extend(index_list)
            sum_of_words = 0
            last_document = name
            index_list = list()
        sum_of_words += 1
        index_list.append(temp_index)
        temp_index += 1
    if(sum_of_words < threshold):
        drop_list.extend(index_list)
    pre_data.drop(drop_list, inplace=True)
    return pre_data

def output_to_file(pre_data, cancer_type, number_of_document):
    output_file = 'data/data_' + cancer_type + '.txt'
    output = open(output_file, 'w')

    name_list = pre_data['Sample name']
    last_document = name_list[0]
    data_mat = [[] for i in range(number_of_document)]
    index = 0
    error = 0
    document = 0
    for name in name_list:
        if(name != str(last_document)):
            last_document = name
            document += 1
        Mutation = calc_mutation(pre_data['Mutation CDS'][index],pre_data['Mutation genome position'][index],pre_data['Mutation strand'][index])
        if(Mutation == -1):
            error += 1
        else:
            data_mat[document].append(Mutation)
        index += 1
        print(index)
 
    output.write(str(number_of_document) + '\n')
    for i in range(number_of_document):
        for j in range(len(data_mat[i])):
            if(j != len(data_mat[i]) - 1):
                output.write(data_mat[i][j] + ' ')
            else:
                output.write(data_mat[i][j])
        output.write('\n')

def calc_mutation(mutation, position, strand):
    before = mutation[len(mutation)-3]
    after = mutation[len(mutation)-1]
    position_list = re.split(r'[:-]',position)
    if(len(position_list) != 3):
        print('position error')
        return -1
    if(int(position_list[0]) == 23):
        chromosome = 'X'
    elif(int(position_list[0]) == 24):
        chromosome = 'Y'
    elif(int(position_list[0]) == 25):
        chromosome = 'M'
    else:
        chromosome = int(position_list[0])
    start = int(position_list[1])
    num = int(position_list[2]) - int(position_list[1]) + 1
    GRCh_file = '/Users/tarom/Data/Mutation_Signature/chr' + str(chromosome) + '.fa'
    quotient = start // 50
    surplus = start % 50

    if(surplus != 0):
        target_index = int(surplus) - 1
    else:
        quotient -= 1
        target_index = 49
    target_line = linecache.getline(GRCh_file, int(quotient)+1)
    pre_line = linecache.getline(GRCh_file, int(quotient))
    post_line = linecache.getline(GRCh_file, int(quotient)+2)

    if(((target_line[target_index] != before) and (strand == '+')) or ((target_line[target_index] != swap(before))and(strand == '-'))):
        print('error: ' + mutation)
        print('target: ' + target_line[target_index])
        print('strand: ' + strand)
        strand = swap(strand)
        if(((target_line[target_index] != before) and (strand == '+')) or ((target_line[target_index] != swap(before))and(strand == '-'))):
            print('still error')
            return -1

    bases = pre_line[:-1] + target_line[:-1] + post_line[:-1]
    target_index += 50
    print('bases:',bases)

    Upstream = bases[target_index-10:target_index]
    Downstream = bases[target_index+1:target_index+11]
    print('Upstream:',Upstream)
    print('Downstream:',Downstream)

    if(((strand == '+') and (before in ['A', 'G'])) or ((strand == '-') and (before in ['C', 'T']))):
        buf_Upstream = Upstream
        Upstream = ''
        for i in range(10):
            Upstream += swap(Downstream[9-i])
        Downstream = ''
        for i in range(10):
            Downstream += swap(buf_Upstream[9-i])
    if(before in ['A', 'G']):
        before = swap(before)
        after = swap(after)

    Upstream = list(Upstream)
    Downstream = list(Downstream)
    for i in range(10):
        Upstream[i] = coding(Upstream[i])
        Downstream[i] = coding(Downstream[i])
    Word = coding(str('[' + before + '>' + after + ']'))
    answer = ''
    for i in range(10):
        answer += str(Upstream[i])
    answer += str(Word)
    for i in range(10):
        answer += str(Downstream[i])
    answer += ',' + position
    return(answer)

def swap(base):
    if(base == 'A'):
        return('T')
    elif(base == 'C'):
        return('G')
    elif(base == 'G'):
        return('C')
    elif(base == 'T'):
        return('A')
    elif(base == '+'):
        return('-')
    elif(base == '-'):
        return('+')
    else:
        return(base)

if __name__ == '__main__':
    main()
