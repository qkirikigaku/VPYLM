import os

def main():
    PATH = 'data/'
    cancer_type_list = []
    for x in os.listdir(PATH):
        if x.startswith('data_') and x.endswith('.txt'):
            index = x.index('.txt')
            cancer_type_list.append(x[5:index])
    for cancer_type in cancer_type_list:
        input_file = 'data/data_' + cancer_type + '.txt'
        lines = open(input_file, 'r').readlines()
        corpus = []
        for i,line in enumerate(lines):
            if(i != 0):
                words = line.split()
                for word in words:
                    mutation = word.split(',')[0]
                    corpus.append(mutation)
        output_file_name = 'data/corpus_' + cancer_type + '.txt'
        output = open(output_file_name, 'w')
        for mutation in corpus:
            output.write(str(int(mutation[10])-4) +  ' ')
            for i in range(10):
                output.write(mutation[9-i])
            output.write(' ' + mutation[11:])
            output.write('\n')
        output.close()

if __name__ == '__main__':
    main()
