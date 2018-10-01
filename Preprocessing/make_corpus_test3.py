import numpy as np
from coding import coding
from coding import decoding

N = 5000

mutation_prob = [0.2, 0.2, 0.2, 0.2, 0.1, 0.1]

out_prob = [[0.90, 0.09, 0.005, 0.005], \
            [0.005, 0.90, 0.09, 0.005], \
            [0.005, 0.005, 0.90, 0.09], \
            [0.09, 0.005, 0.005, 0.90], \
            [0.09, 0.90, 0.005, 0.005], \
            [0.90, 0.005, 0.005, 0.09]]

out_uniform = [0.25, 0.25, 0.25, 0.25]

num_context = 10

stop_prob = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]

def main():
    output_file_name = 'data/corpus_test3.txt'
    output = open(output_file_name, 'w')
    for i in range(N):
        output.write(make_mutation())
    output.close()

def make_mutation():
    mutation = ''
    substitution = choose_substitution()
    mutation = str(substitution)
    for i in range(2):
        mutation += ' '
        stop_flag = 0
        for j in range(num_context):
            if(stop_flag == 0): stop_flag = np.random.binomial(1, stop_prob[j])
            if(stop_flag == 0): temp_context = choose_base(substitution)
            else: temp_context = choose_base(-1)
            mutation += str(temp_context)
    mutation += '\n'
    return mutation

def choose_substitution():
    array = np.random.multinomial(1, mutation_prob)
    for i,x in enumerate(array):
        if(x == 1): return i

def choose_base(substitution):
    prob = list()
    if(substitution == -1):
        prob = out_uniform
    else: prob = out_prob[substitution]
    array = np.random.multinomial(1, prob)
    for i,x in enumerate(array):
        if(x == 1): return i

if __name__ == '__main__':
    main()
