import sys
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    args = sys.argv
    cancer_type = args[1]
    direction = args[2]
    File_name = 'result/' + cancer_type + '/corpus_' + direction + '_log_likelihood.txt'
    height = list()
    for line in open(File_name, 'r').readlines():
        height.append(float(line))
    left = range(0, len(height))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(left, height)
    ax1.set_xlim(0, len(height))
    ax1.set_xlabel('Sampling iteration')
    ax1.set_ylabel('Log-likelihood')
    ax1.set_title('Transition of log-likelihood')
    fig.tight_layout()
    name = 'result/' + cancer_type + '/@transition_likelihood_' + direction + '.png'
    fig.savefig(name, dpi=300)
    plt.close(1)

if __name__ == '__main__':
    main()
