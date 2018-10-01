import sys
import os
from datetime import datetime
from time import sleep

MAX_EXP_NUM = 99

def main():
    cancer_types = load_cancer_type()
    while(1):
        complete_list = list()
        for cancer_type in cancer_types:
            end_num = 0
            for i in range(MAX_EXP_NUM):
                start_unix_time = os.stat(error_log_file \
                                + '.' + str(i+1)).st_mtime
                base_file_name = 'result/' + cancer_type + '_' + str(i+1) \
                               + '/corpus_' + cancer_type + '_'
                for direction in ['Up', 'Down']:
                    file_name = base_file_name + direction + '.txt'
                    end_num += check_complete(file_name, start_unix_time)
            print(cancer_type + ' : ' + str(end_num) + ' / ' \
                  + str(MAX_EXP_NUM*2) + ' jobs completed.\n')
            if(end_num == MAX_EXP_NUM*2): complete_list.append(cancer_type)
        print('\n')
        print('Completed : ', complete_list)
        print('\n')
        sleep(60)

def load_cancer_type():
    file_name = 'data/primary_lesion_list.txt'
    cancer_types = list()
    for line in open(file_name, 'r').readlines():
        cancer_types.append(line[:-1])
    return cancer_types

def check_complete(file_name, start_unix_time):
    if(os.path.exists(file_name) == False): return 0
    else:
        update_unix_time = os.stat(file_name).st_mtime
        if(update_unix_time > start_unix_time): return 1
        else: return 0

if __name__ == '__main__':
    args = sys.argv
    error_log_file = args[1]
    main()
