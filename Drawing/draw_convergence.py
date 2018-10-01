import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    test1 = load_data('VPYLM.sh.o86846841.')
    test2 = load_data('VPYLM.sh.o86846843.')
    test3 = load_data('VPYLM.sh.o86860820.')
    height1 = np.zeros([499])
    height2 = np.zeros([499])
    height3 = np.zeros([499])
    for i in range(499):
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for j in range(20):
            for k in range(11):
                sum1 += abs(test1[j][i+1][k] - test1[j][i][k])
                sum2 += abs(test2[j][i+1][k] - test2[j][i][k])
                sum3 += abs(test3[j][i+1][k] - test3[j][i][k])
        sum1 /= 20
        sum2 /= 20
        sum3 /= 20
        height1[i] = sum1
        height2[i] = sum2
        height3[i] = sum3
    left = range(1,500)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(left, height1, label='test1')
    ax1.plot(left, height2, label='test2')
    ax1.plot(left, height3, label='test3')
    ax1.set_xlim(0, 500)
    ax1.legend()
    fig.tight_layout()
    fig.savefig('result/@other/@figure/convergence.png', dpi=300)
    plt.close(1)

def load_data(File_head):
    transition_list = list()
    for i in range(20):
        transition_list.append(list())
        File_name = 'result/@other/' + File_head + str(i+1)
        File = open(File_name, 'r')
        lines = File.readlines()
        for j in range(500):
            transition_list[i].append(list())
            for k in range(11):
                if(k != 10):
                    transition_list[i][j].append(int(lines[15 * j + k + 8][4:]))
                else:
                    transition_list[i][j].append(int(lines[15 * j + k + 8][5:]))
    return transition_list

if __name__ == '__main__':
    main()
