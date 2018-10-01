import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from coding import coding
from coding import decoding

num_vocabulary = 6
num_base = 4
num_max_context = 10

def main():
    args = sys.argv
    cancer_type = args[1]
    mutation,upstream,downstream = load_data(cancer_type)
    up_ngram, up_posterior_sampling = load_result(cancer_type, 'Up')
    down_ngram, down_posterior_sampling = load_result(cancer_type, 'Down')
    draw_fig_ngram('Up', up_ngram, cancer_type)
    draw_fig_ngram('Down', down_ngram, cancer_type)
    draw_fig_ngram_by_mutation('Up', mutation, up_ngram, cancer_type, upstream)
    draw_fig_ngram_by_mutation('Down', mutation, down_ngram, cancer_type, downstream)
    output_context('Up', mutation, up_ngram, cancer_type, upstream)
    output_context('Down', mutation, down_ngram, cancer_type, downstream)
    output_posterior_sampling('Up', up_posterior_sampling, cancer_type)
    output_posterior_sampling('Down', down_posterior_sampling, cancer_type)

def load_data(cancer_type):
    file_name = 'data/corpus_' + cancer_type + '.txt'
    lines = open(file_name, 'r').readlines()
    mutation = list()
    upstream = list()
    downstream = list()
    for line in lines:
        temp_list = line.split()
        mutation.append(int(temp_list[0]))
        temp_upstream = list()
        for i in range(len(temp_list[1])):
            temp_upstream.append(int(temp_list[1][i]))
        upstream.append(temp_upstream)
        temp_downstream = list()
        for i in range(len(temp_list[2])):
            temp_downstream.append(int(temp_list[2][i]))
        downstream.append(temp_downstream)
    return mutation,upstream,downstream

def load_result(cancer_type, direction):
    file_name = 'result/' + cancer_type + '/corpus_' + direction + '.txt'
    File = open(file_name, 'r')
    File.readline()
    temp_list = File.readline().split()
    ngram = list()
    for x in temp_list:
        ngram.append(int(x))
    posterior_sampling = list()
    for i in range(num_vocabulary):
        temp_list = File.readline().split()
        posterior_sampling.append([])
        for x in temp_list:
            posterior_sampling[i].append(x)
    return ngram, posterior_sampling

def draw_fig_ngram(direction, ngram, cancer_type):
    ngram_count = np.zeros(num_max_context + 1)
    for x in ngram:
        ngram_count[x] += 1
    left = range(num_max_context + 1)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.bar(left, ngram_count, width=1, align='center')
    ax1.set_xlim(-1, num_max_context + 1)
    ax1.set_xticks(left)
    ax1.set_xlabel('predicted n-gram')
    ax1.set_ylabel('number of mutations')
    fig.tight_layout()
    fig.savefig('result/' + cancer_type \
                + '/ngram_count_' + direction + '.png', dpi=300)
    plt.close(1)

def draw_fig_ngram_by_mutation(direction, mutation, ngram, cancer_type, context):
    new_ngram = list()
    vocab_list = ["[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]"]
    for i in range(num_vocabulary):
        new_ngram.append([])
    for i,x in enumerate(mutation):
        new_ngram[x].append(ngram[i])
    
    max_ngram = list()
    second_ngram = list()
    ngram_count = np.zeros([num_vocabulary, num_max_context+1])
    for i in range(num_vocabulary):
        for x in new_ngram[i]:
            ngram_count[i,x] += 1
        max_ngram.append(np.argmax(ngram_count[i]))
    for i,x in enumerate(max_ngram): # max_ngramがユニグラムでなければ該当文脈を参照するようにする。
        max_ngram[i] -= 1
    
    fig = plt.figure()
    ax1 = fig.add_subplot(2,3,1); ax2 = fig.add_subplot(2,3,2); ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4); ax5 = fig.add_subplot(2,3,5); ax6 = fig.add_subplot(2,3,6)
    axs = plt.gcf().get_axes()
    for i,ax in enumerate(axs):
        left = range(num_max_context + 1)
        ax.bar(left, ngram_count[i], width=1, align='center')
        ax.set_xlim(-1, num_max_context + 1)
        ax.set_xticks(left)
        ax.set_xlabel('predicted n-gram order')
        ax.set_ylabel('number of mutations')
        ax.set_title(vocab_list[i] + ' (' + str(int(sum(ngram_count[i]))) + ')')
    fig.tight_layout()
    fig.savefig('result/' + cancer_type + '/ngram_count_' + direction + '_by_mutation.png', dpi=300)
    plt.close(1)

    # 推定されたn-gramの中でもっとも多いnオーダーに該当する塩基の割合を出力
    height = np.zeros([num_vocabulary, num_base])
    base_list = ['A', 'C', 'G', 'T']
    for i,x in enumerate(mutation):
        if(max_ngram[x] != -1): height[x, context[i][max_ngram[x]]] += 1
    fig = plt.figure()
    ax1 = fig.add_subplot(2,3,1); ax2 = fig.add_subplot(2,3,2); ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4); ax5 = fig.add_subplot(2,3,5); ax6 = fig.add_subplot(2,3,6)
    axs = plt.gcf().get_axes()
    for i,ax in enumerate(axs):
        left = range(num_base)
        if(max_ngram[i] != -1): 
            ax.bar(left, height[i], width=1, align='center')
            ax.set_xlim(-1, num_base)
            ax.set_xticks(left)
            ax.set_xticklabels(base_list)
            ax.set_xlabel('Base type')
            ax.set_ylabel('Proportion of the base apart ' + str(max_ngram[i] + 1) + ' bases')
            ax.set_title(vocab_list[i] + ' (' + str(int(sum(height[i]))) + ')')
        else:
            ax.bar(0,0, width=1, align='center')
            ax.set_title(vocab_list[i] + ' is estimated with unigram.')
    fig.tight_layout()
    fig.savefig('result/' + cancer_type + '/Nth_context_' + direction + '.png', dpi=300)
    plt.close(1)

def output_context(direction, mutation, ngram, cancer_type, context):
    dic_list = [{} for i in range(num_vocabulary)]
    for i,x in enumerate(ngram):
        temp_mutation = mutation[i]
        temp_context_list = context[i][:x]
        temp_context = ''
        for y in temp_context_list:
            temp_context += str(y)
        if(temp_context in dic_list[temp_mutation]): 
            dic_list[temp_mutation][temp_context] += 1
        else: dic_list[temp_mutation].update({temp_context:1})
    fig = plt.figure()
    for i in range(num_vocabulary):
        fig.add_subplot(2,3,i+1)
    axs = plt.gcf().get_axes()
    for i,ax in enumerate(axs):
        left_label = list()
        height = list()
        for k,v in sorted(dic_list[i].items(), key=lambda x: -x[1]):
            if(direction == 'Up'):
                name = '5\'-'
                for j,y in enumerate(list(k)[::-1]):
                    name += decoding(y)
                name += decoding(str(i+4))
                left_label.append(name)
                height.append(v)
            elif(direction == 'Down'):
                name = decoding(str(i+4))
                for j,y in enumerate(list(k)):
                    name += decoding(y)
                name += '-3\''
                left_label.append(name)
                height.append(v)
        left_label = left_label[:10]
        shaped_height = height[:10]
        left = range(len(left_label))
        ax.bar(left, shaped_height, width=1, align='center')
        ax.set_xticks(left)
        ax.set_xticklabels(left_label)
        for tick in ax.get_xticklabels():
            tick.set_rotation(85)
        ax.set_xlabel('Predicted Context')
        ax.set_ylabel('Count')
        ax.set_title(decoding(str(i+4)) + ' (' + str(int(sum(height))) + ')')
    fig.tight_layout()
    fig.savefig('result/' + cancer_type + '/predicted_context_' + direction + '.png', dpi=300)
    plt.close(1)

def output_posterior_sampling(direction, sampling, cancer_type):
    dic_list = list()
    fig = plt.figure()
    ax1 = fig.add_subplot(2,3,1); ax2 = fig.add_subplot(2,3,2); ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4); ax5 = fig.add_subplot(2,3,5); ax6 = fig.add_subplot(2,3,6)
    axs = plt.gcf().get_axes()
    for i,ax in enumerate(axs):
        dic_list.append({})
        num_all = len(sampling[i])
        for x in sampling[i]:
            if(x in dic_list[i]):
                dic_list[i][x] += 1
            else:
                dic_list[i].update({x:1})
        left_label = list()
        height = list()
        for k,v in sorted(dic_list[i].items(), key=lambda x: -x[1]):
            if(direction == 'Up'): 
                name = '5\'-'
                for j,y in enumerate(list(k)[::-1]):
                    if(j != len(list(k))-1): name += decoding(y)
                name += decoding(str(i+4))
                left_label.append(name)
                height.append(v)
            elif(direction == 'Down'): 
                name = decoding(str(i+4))
                for j,y in enumerate(list(k)):
                    if(j != 0): name += decoding(y)
                name += '-3\''
                left_label.append(name)
                height.append(v)
        left_label = left_label[:10]
        height = height[:10]
        left = range(len(left_label))
        ax.bar(left, height, width=1, align='center')
        ax.set_xlim(-1, len(left_label)+1)
        ax.set_xticks(left)
        ax.set_xticklabels(left_label)
        for tick in ax.get_xticklabels():
            tick.set_rotation(85)
        ax.set_xlabel('Probabilistic phrase')
        ax.set_ylabel('Count')
        ax.set_title(decoding(str(i+4)))
    fig.tight_layout()
    fig.savefig('result/' + cancer_type + '/posterior_sampling_' + direction + '.png', dpi=300)
    plt.close(1)

if __name__ == '__main__':
    main()
