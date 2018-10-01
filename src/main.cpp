//
// Created by T. Matsutani.
//

#include "VPYLM.h"

void run(string cancer_type, string direction, int exp_index);
int ctoi(const char c);

class TREES{
    int num_mutation;
    vector<int> mutation;
    string cancer_type;
    vector<vector<int> > upstream;
    vector<vector<int> > downstream;
public:
    int exp_index;
    string direction;
    VPYLM* direction_tree;

    TREES(string x, string y, int z);
    ~TREES();
    bool load_data();
    void set_vpylm();
};

TREES::TREES(string x, string y, int z){
    cancer_type = x;
    direction = y;
    exp_index = z;
}

TREES::~TREES(){
    delete direction_tree;
}

bool TREES::load_data(){
    ifstream ifs;
    string input_file_name = "data/corpus_" + cancer_type + ".txt";
    ifs.open(input_file_name.c_str(), ios::in);
    if (!ifs){
        cout << "Cannot open " + input_file_name << endl;
        return false;
    }
    char buf[256];
    char *temp;
    string temp_context;
    while (ifs.getline(buf,256)){
        temp = strtok(buf, " ");
        mutation.push_back(atoi(temp));
        temp_context = strtok(NULL, " ");
        vector<int> temp_upstream;
        for (int i=0; i < MAX_CONTEXT_NUM; i++){
            int temp_num = ctoi(temp_context[i]);
            temp_upstream.push_back(temp_num);
        }
        upstream.push_back(temp_upstream);
        temp_context = strtok(NULL, " ");
        vector<int> temp_downstream;
        for (int i=0; i < MAX_CONTEXT_NUM; i++){
            int temp_num = ctoi(temp_context[i]);
            temp_downstream.push_back(temp_num);
        }
        downstream.push_back(temp_downstream);
    }
    num_mutation = mutation.size();
    return true;
}

void TREES::set_vpylm(){
    if(direction == "Upstream"){
        direction_tree = new VPYLM(cancer_type, "Up", mutation, upstream, exp_index);
    }
    else if(direction == "Downstream"){
        direction_tree = new VPYLM(cancer_type, "Down", mutation, downstream, exp_index);
    }
}

int main(int argc, char * argv[]){
    if (argc != 4){
        cout << "The number of argument is invalid." << endl;
        return(1);
    }
    string cancer_type = argv[1]; // = (string) primary lesion (e.g. lung).
    string direction = argv[2];
    int exp_index = atoi(argv[3]);
    run(cancer_type, direction, exp_index);
}

void run(string cancer_type, string direction, int exp_index){
    TREES trees(cancer_type, direction, exp_index);
    trees.load_data();
    trees.set_vpylm();
    cout << "-------------------" << endl;
    cout << "Setting " + direction + "..." << endl;
    cout << "-------------------" << endl;
    trees.direction_tree -> run_learning();
    trees.direction_tree -> write_data();
}

int ctoi(const char c){
    if('0' <= c && c<= '9') return(c-'0');
    return -1;
}

