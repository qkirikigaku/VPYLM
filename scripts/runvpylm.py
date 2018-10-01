import sys
import os
from multiprocessing import Pool
import subprocess
from datetime import datetime

def main():
    primary_lesions = load_lesions()
    arguments = list()
    for x in primary_lesions:
        arguments.append((x, 'Upstream'))
        arguments.append((x, 'Downstream'))
    pool = Pool()
    _ = pool.starmap(execute, arguments)

def load_lesions():
    file_name = 'data/primary_lesion_list.txt'
    primary_lesions = list()
    for line in open(file_name, 'r').readlines():
        primary_lesions.append(line[:-1])
    return primary_lesions

def execute(primary_lesion, direction):
    result_path = 'result/' + primary_lesion + '_' + experiment_index
    if(os.path.exists(result_path) == False): os.mkdir(result_path)
    cmd = './VPYLM ' + primary_lesion + ' ' + direction + ' ' + experiment_index
    subprocess.call(cmd.split())
    finstr  = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    finstr += 'arg: ' + primary_lesion + ', ' + direction
    print(finstr, file=sys.stderr)

if __name__ == '__main__':
    args = sys.argv
    experiment_index = args[1]
    main()
