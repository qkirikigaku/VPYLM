//
// Created by T. Matsutani.
//

#include <string.h>
#include "sampler.h"

#define MAX_CONTEXT_NUM 10

using namespace std;

class Node{
public:
    string direction; // = {Up, Down};
    int base; // = {0,1,2,3} : {A,C,G,T};
    int depth; // the depth of n-gram node (root -> depth = 0).
    
    Node *parent; // = parent ID;
    vector<Node *> children; // = children IDs;
    
    int num_vocabulary; // = 6
    double alpha; // hyper parameter 1;
    double beta; // hyper parameter 2;
    double a; // # mutations stopped at this node;
    double b;// # mutations that passed this node;
    double initial_stop_probability; // q_i

    void set_a(string x, int y, int z);
    void set_parent(Node *x);
    bool add_child(Node *x);
    bool exist_next_node(int x); // x : the type of mutation;
    Node* find_next_node(int x); // x : the type of mutation;

    vector<double> p_m; // = p(m|seat_arrangement);
    vector<vector<int> > count; // = C_uwk;
    vector<int> table; // = t_uw;
    int num_all_table; // = t_u;

    double rho; // rho, sigma, and tau are auxiliary variables to calculate hyperparameter.
    vector<double> sigma;
    vector<vector<vector<double> > > tau;

    void add_table(int mutation);
    void remove_table();
    void update_table();
    void calc_auxiliary_variables(vector<double> &dis, vector<double> &theta);
};

void Node::set_a(string x, int y, int z){
    direction = x;
    base = y;
    depth = z;
    children.reserve(4);
    num_vocabulary = 6;
    alpha = 1.0;
    beta = 5.0;
    a = 0.0;
    b = 0.0;
    initial_stop_probability = sampler::beta(alpha, beta);
    if (depth == MAX_CONTEXT_NUM) initial_stop_probability = 1.0;
    p_m.resize(num_vocabulary, 1.0/num_vocabulary);
    count.resize(num_vocabulary);
    table.resize(num_vocabulary, 0);
    num_all_table = 0;
    tau.resize(num_vocabulary);
}

void Node::set_parent(Node *x){
    parent = x;
}

bool Node::add_child(Node *x){
    children.push_back(x);
    if(children.size() > 4) return false;
    else return true;
}

bool Node::exist_next_node(int x){
    for (int i=0; i < children.size(); i++){
        if(children[i] -> base == x) return true;
    }
    return false;
}

Node* Node::find_next_node(int x){
    Node* index;
    for (int i=0; i < children.size(); i++){
        if(children[i] -> base == x) index = children[i];
    }
    return index;
}

void Node::add_table(int x){
    for (int i=0; i < num_vocabulary; i++){
        if(i == x) count[i].push_back(1);
        else count[i].push_back(0);
    }
    update_table();
}

void Node::remove_table(){
    int erase_index = -1;
    for (int i=0; i < num_all_table; i++){
        int erase_flag = 1;
        for (int j=0; j < num_vocabulary; j++){
            if(count[j][i] != 0) erase_flag = 0;
        }
        if(erase_flag == 1) erase_index = i;
    }
    if (erase_index != -1){
        for (int i=0; i < num_vocabulary; i++){
            count[i].erase(count[i].begin() + erase_index);
        }
    }
    update_table();
}

void Node::update_table(){
    num_all_table = count[0].size();
    for (int i=0; i < num_vocabulary; i++){
        table[i] = 0;
        for (int j=0; j < num_all_table; j++){
            if(count[i][j] != 0) table[i]++;
        }
    }
}

void Node::calc_auxiliary_variables(vector<double> &dis, vector<double> &theta){
    int num_all_customers = 0;
    for (int i=0; i < num_vocabulary; i++){
        for (int j=0; j < num_all_table; j++){
            num_all_customers += count[i][j];
        }
    }
    rho = sampler::beta(theta[depth]+1, num_all_customers-1);
    
    if(sigma.size() != 0) sigma.clear();
    for (int i=0; i < num_all_table-1; i++){
        sigma.push_back(sampler::bernoulli((theta[depth])/(theta[depth] + dis[depth] * (i+1))));
    }
    
    for (int i=0; i < num_vocabulary; i++){
        if(tau[i].size() != 0) tau[i].clear();
        for (int j=0; j < num_all_table; j++){
            vector<double> vec;
            tau[i].push_back(vec);
            for (int k=0; k < count[i][j]-1; k++){
                tau[i][j].push_back(sampler::bernoulli(k/(k+1-dis[depth])));
            }
        }
    }
}
