import sys
import os
import shutil
import glob

Max_index = 100

def main():
    args = sys.argv
    cancer_type = args[1]
    up_index = find_index('Up', cancer_type)
    down_index = find_index('Down', cancer_type)
    PATH = 'result/' + cancer_type
    if(os.path.exists(PATH) == False): os.mkdir(PATH)
    up_file = 'result/' + cancer_type + '_' + str(up_index) \
            + '/corpus_' + cancer_type + '_Up.txt'
    down_file = 'result/' + cancer_type + '_' + str(down_index) \
            + '/corpus_' + cancer_type + '_Down.txt'
    up_likelihood = 'result/' + cancer_type + '_' + str(up_index) \
                  + '/corpus_' + cancer_type + '_Up_log_likelihood.txt'
    down_likelihood = 'result/' + cancer_type + '_' + str(down_index) \
                  + '/corpus_' + cancer_type + '_Down_log_likelihood.txt'
    shutil.copyfile(up_file, PATH + '/corpus_Up.txt')
    shutil.copyfile(down_file, PATH + '/corpus_Down.txt')
    shutil.copyfile(up_likelihood, PATH + '/corpus_Up_log_likelihood.txt')
    shutil.copyfile(down_likelihood, PATH + '/corpus_Down_log_likelihood.txt')
    rm_list = glob.glob('result/' + cancer_type + '_*')
    for x in rm_list:
        shutil.rmtree(x)

def find_index(direction, cancer_type):
    likelihood_list = list()
    index_list = list()
    for i in range(Max_index):
        File_name = 'result/' + cancer_type + '_' + str(i+1) \
                  + '/corpus_' + cancer_type + '_' + direction + '.txt'
        if(os.path.isfile(File_name)):
            File = open(File_name, 'r')
            likelihood = float(File.readline())
            if(likelihood < -100):
                likelihood_list.append(likelihood)
                index_list.append(i)
    print(direction)
    print(likelihood_list)
    return index_list[likelihood_list.index(max(likelihood_list))]+1

if __name__ == '__main__':
    main()
