import os

def main():
   file_list = os.listdir('data/')
   files = list()
   for File in file_list:
       if(File.startswith('corpus_') and File.endswith('.txt') and File.find('test') == -1):
           files.append(File)
   output_file = "data/corpus_all.txt"
   output = open(output_file, 'w')
   for File in files:
       input_file = open("data/" + File, 'r')
       for line in input_file.readlines():
           output.write(line)

if __name__ == '__main__':
    main()
