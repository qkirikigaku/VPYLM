import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
from draw_ngram import load_data

MAX_INDEX = 100
NUM_VOCAB = 6
MAX_NUM_CONTEXT = 10

def main():
    args = sys.argv
    cancer_type = args[1]
    mutation, upstream, downstream = load_data(cancer_type)
    up_LL, up_error, up_num_mutation \
    = load_exp('Up', cancer_type, mutation)
    down_LL, down_error, down_num_mutation \
    = load_exp('Down', cancer_type, mutation)
    PATH = 'result/' + cancer_type
    if(os.path.exists(PATH) == False): os.mkdir(PATH)
    draw_scatter("Up", cancer_type, up_LL, up_error, up_num_mutation)
    draw_scatter("Down", cancer_type, down_LL, down_error, down_num_mutation)

def load_exp(direction, cancer_type, mutation):
    likelihood_list = list()
    ngram_list = list()
    index = 0
    for i in range(MAX_INDEX):
        File_name = 'result/' + cancer_type + '_' + str(i+1) + '/corpus_' \
                  + cancer_type + '_' + direction + '.txt'
        if(os.path.isfile(File_name)):
            File = open(File_name, 'r')
            temp_likelihood = float(File.readline())
            if(temp_likelihood < -100):
                likelihood_list.append(temp_likelihood)
                line = File.readline().split()
                ngram_list.append([])
                for word in line:
                    ngram_list[index].append(int(word))
                index += 1
    true_ngram = ngram_list[likelihood_list.index(max(likelihood_list))]
    num_mutation = len(ngram_list[0])
    shaped_true_ngram = [[]for i in range(NUM_VOCAB)]
    for i,x in enumerate(true_ngram):
        shaped_true_ngram[mutation[i]].append(x)
    for i in range(NUM_VOCAB):
        shaped_true_ngram[i] = shape(shaped_true_ngram[i])

    error_list = np.zeros([len(ngram_list)])
    for i,x in enumerate(ngram_list):
        shaped_ngram = [[]for i in range(NUM_VOCAB)]
        for j,y in enumerate(x):
            shaped_ngram[mutation[j]].append(y)
        for j in range(NUM_VOCAB):
            shaped_ngram[j] = shape(shaped_ngram[j])
        for j,y in enumerate(shaped_ngram):
            for k,z in enumerate(y):
                error_list[i] += abs(z - shaped_true_ngram[j][k])
    return likelihood_list, error_list, num_mutation

def shape(ngram):
    shaped_list = [0 for i in range(MAX_NUM_CONTEXT + 1)]
    for x in ngram:
        shaped_list[x] += 1
    return shaped_list

def draw_scatter(direction, cancer_type, LL, error, num_mutation):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(LL, error)
    ax.set_xlabel('log-likelihood')
    ax.set_ylabel('error')
    ax.set_title(cancer_type + ' ' + direction \
                 + ' (' + str(num_mutation) + ' mutations)')
    fig.tight_layout()
    result_name = 'result/' + cancer_type + '/@LL-error_' + direction + '.png'
    fig.savefig(result_name, dpi=300)
    plt.close(1)

if __name__ == '__main__':
    main()
