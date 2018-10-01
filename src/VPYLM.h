//
// Created by T. Matsutani.
//

#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <numeric>
#include "node.h"
#include "sampler.h"

using namespace std;

class VPYLM{
    string cancer_type;
public:
    Node* root;
    vector<Node*> nodes;
    int num_mutation;
    string direction; // = {Up, Down}
    vector<int> mutation; // mutation[i]; (1 <= i <= num_mutation)
    vector<vector<int> > context; // context[i][MAX_CONTEXT_NUM];
    vector<Node*> real_customer_arrangement; // 1 <= i <= num_mutation;

    vector<double> dis; // discount parameter[depth]; 0 <= d[depth] < 1, 
    vector<double> theta; // strength parameter[depth]; theta[depth] > -d[depth]

    int exp_index;

    VPYLM(string x, string y, vector<int> &z, vector<vector<int> > &w, int v);
    ~VPYLM();
    void run_learning();
    void make_suffix_tree();
    void initialize();

    void add_customer(int index, Node* position);
    
    void calc_p_m(Node* start_node);
    
    void remove_customer(int index, Node* position);
    
    void reset_ab(int index, Node* position);
    void arrange_real_customer(int index);
    
    void update_hyperparameter();
    void update_auxiliary_variables(Node* start_node);
    void calc_argument(Node* start_node, vector<double> &dis1,
                       vector<double> &dis2, vector<double> &theta1, 
                       vector<double> &theta2);
    
    double SumVec(vector<double> &Vec);
    void Normalize(vector<double> &Vec);
    
    double temp_log_likelihood;
    double old_log_likelihood;
    void calc_log_likelihood();
    vector<double> vec_log_likelihood;

    void write_data();
    string posterior_sampling(int mutation);

    void print_depth_distribution();
    void print_root_count(int iter, int index, string status);
};

VPYLM::VPYLM(string x, string y, vector<int> &z,
             vector<vector<int> > &w, int v){
    cancer_type = x;
    direction = y;
    mutation = z;
    context = w;
    exp_index = v;
    num_mutation = context.size();
    real_customer_arrangement.resize(num_mutation);
    root = new Node(); root -> set_a(direction, -1, 0);
    for (int i=0; i < MAX_CONTEXT_NUM+1; i++){
        dis.push_back(sampler::uniform(0.0, 0.999999));
        theta.push_back(sampler::uniform(-dis[i] + 1.0000, -dis[i] + 2.0000));
    }
    old_log_likelihood = -10000000000;
}

VPYLM::~VPYLM(){
    delete root;
    for (int i=0; i < nodes.size(); i++){
        delete nodes[i];
    }
}

void VPYLM::run_learning(){
    make_suffix_tree();
    initialize();
    cout << "---------------" << endl;
    cout << "Sampling Start." << endl;
    cout << "---------------" << endl;
    for (int iter=0; iter < 100; iter++){
        cout << "iteration : "  << to_string(iter) << endl;
        print_depth_distribution();
        calc_p_m(root);
        for(int index=0; index < num_mutation; index++){
            remove_customer(index, real_customer_arrangement[index]);
            reset_ab(index, real_customer_arrangement[index]);
            arrange_real_customer(index);
            add_customer(index, real_customer_arrangement[index]);
        }
        update_hyperparameter();
        calc_log_likelihood();
        cout << "temp log_likelihood : " 
             << to_string(temp_log_likelihood) << endl;
        cout << "Improved : " 
             << to_string(temp_log_likelihood - old_log_likelihood) << endl;
        old_log_likelihood = temp_log_likelihood;
        cout << "completed" << endl;
        cout << endl;
    }
}

void VPYLM::make_suffix_tree(){
    for (int i=0; i < num_mutation; i++){
        Node *parent = root;
        for (int j=0; j < MAX_CONTEXT_NUM; j++){
            int base = context[i][j];
            if(parent -> exist_next_node(base)){
                parent = parent -> find_next_node(base);
            }
            else{
                Node *node = new Node();
                nodes.push_back(node);
                node -> set_a(direction, base, j+1); node -> set_parent(parent);
                parent -> add_child(node); parent = node;
            }
        }
    }
}

void VPYLM::initialize(){
    for (int i=0; i < num_mutation; i++){
        Node *parent = root;
        int temp_mutation = mutation[i];
        int j=0;
        while(j < MAX_CONTEXT_NUM+1){
            int stop_flag 
                = sampler::bernoulli(parent -> initial_stop_probability);
            if(stop_flag == 1){
                real_customer_arrangement[i] = parent;
                parent -> a[temp_mutation] += 1.0;
                break;
            }
            else{
                int base = context[i][j];
                Node *child = parent -> find_next_node(base);
                parent -> b[temp_mutation] += 1.0;
                parent = child;
            }
            j++;
        }
        add_customer(i, real_customer_arrangement[i]);
    }
}

void VPYLM::add_customer(int index, Node* node){
    int temp_mutation = mutation[index];
    //calc_p_m(root);
    vector<double> p;
    if(node != root){
        int dep = node -> depth;
        for (int i=0; i < node -> num_all_table; i++){
            double candidate = node -> count[temp_mutation][i] 
                             - dis[dep] * node -> table[temp_mutation];
            p.push_back(max(0.0, candidate));
        }
        double component = theta[dep] + dis[dep] * node -> num_all_table;
        p.push_back(component * node -> parent -> p_m[temp_mutation]);
        Normalize(p);
    }
    else{
        int flag=0;
        for (int i=0; i < node -> num_all_table; i++){
            if(node -> count[temp_mutation][i] != 0){
                p.push_back(1.0); flag=1;
            }
            else p.push_back(0.0);
        }
        if(flag == 0) p.push_back(1.0);
        else p.push_back(0.0);
    }
    int temp_table = sampler::categorical(node -> num_all_table + 1, p);
    if(temp_table != node -> num_all_table){
        node -> count[temp_mutation][temp_table]++;
        calc_p_m(node);
    }
    else{
        node -> add_table(temp_mutation);
        if(node != root) add_customer(index, node -> parent);
    }
};

void VPYLM::calc_p_m(Node* start_node){
    vector<double> parent_p_m;
    if(start_node != root){
        parent_p_m = start_node -> parent -> p_m;
        int dep = start_node -> depth;
        vector<int> count_w;
        count_w.resize(start_node -> num_vocabulary, 0);
        for (int i=0; i < start_node -> num_vocabulary; i++){
            for (int j=0; j < start_node -> num_all_table; j++){
                count_w[i] += start_node -> count[i][j];
            }
        }
        for (int i=0; i < start_node -> num_vocabulary; i++){
            start_node -> p_m[i] = count_w[i] - dis[dep] 
                                 * start_node -> table[i] 
                                 + (theta[dep] + dis[dep] 
                                 * start_node -> num_all_table) * parent_p_m[i];
        }
        Normalize(start_node -> p_m);
    }
    else{
        for (int i=0; i < root -> num_vocabulary; i++){
            start_node -> p_m[i] = 1.0 / start_node -> num_vocabulary;
        }
    }
    for (int i=0; i < start_node -> children.size(); i++){
        Node *node = start_node -> children[i];
        if((SumVec(node -> a) != 0.0) || (SumVec(node -> b) != 0.0)){
            calc_p_m(node);
        }
    }
}

void VPYLM::remove_customer(int index, Node* node){
    int temp_mutation = mutation[index];
    if(node == root){
        int temp_table;
        for (int i=0; i < node -> num_all_table; i++){
            if(node -> count[temp_mutation][i] != 0){
                node -> count[temp_mutation][i]--;
                temp_table = i;
            }
        }
        if(node -> count[temp_mutation][temp_table] == 0){
            node -> remove_table();
        }
        calc_p_m(node);
    }
    else{
        vector<double> p;
        for (int i=0; i < node -> num_all_table; i++){
            p.push_back(node -> count[temp_mutation][i]);
        }
        Normalize(p);
        int temp_table = sampler::categorical(node -> num_all_table, p);
        node -> count[temp_mutation][temp_table]--;
        if(node -> count[temp_mutation][temp_table] == 0){
            node -> remove_table();
            Node* next_node = node -> parent;
            remove_customer(index, next_node);
        }
        else{
            calc_p_m(node);
        }
    }
}

void VPYLM::reset_ab(int index, Node* node){
    int temp_mutation = mutation[index];
    if (node == real_customer_arrangement[index]){
        node -> a[temp_mutation] -= 1.0;
    }
    else node -> b[temp_mutation] -= 1.0;
    if(node != root) reset_ab(index, node -> parent);
}

void VPYLM::arrange_real_customer(int index){
    int temp_mutation = mutation[index];
    vector<double> p_l;
    for(int i=0; i < MAX_CONTEXT_NUM+1; i++){
        Node* node = root;
        double prod_pass = 1.0;
        for (int j=0; j < i; j++){
            prod_pass *= (node -> b[temp_mutation] + node -> beta) 
                        / (node -> a[temp_mutation] + node -> b[temp_mutation]
                            + node -> alpha + node -> beta);
            node = node -> find_next_node(context[index][j]);
        }
        double prob = (node -> a[temp_mutation] + node -> alpha) 
                    / (node -> a[temp_mutation] + node -> b[temp_mutation] 
                        + node -> alpha + node -> beta) * prod_pass;
        prob *= node -> p_m[temp_mutation];
        p_l.push_back(prob);
    }
    Normalize(p_l);
    int stop_node_depth = sampler::categorical(MAX_CONTEXT_NUM+1, p_l);
    Node* node = root;
    int i=0;
    while(1){
        if (node -> depth == stop_node_depth) break;
        node -> b[temp_mutation] += 1.0;
        node = node -> find_next_node(context[index][i]);
        i++;
    }
    node -> a[temp_mutation] += 1.0;
    real_customer_arrangement[index] = node;
}

void VPYLM::update_hyperparameter(){
    update_auxiliary_variables(root);
    const double Kappa = 1.0, Lambda = 1.0, Mu = 1.0, Nu = 1.0;
    vector<double> dis_argument1, dis_argument2,
                   theta_argument1, theta_argument2;
    dis_argument1.resize(MAX_CONTEXT_NUM+1, Kappa); 
    dis_argument2.resize(MAX_CONTEXT_NUM+1, Lambda);
    theta_argument1.resize(MAX_CONTEXT_NUM+1, Mu); 
    theta_argument2.resize(MAX_CONTEXT_NUM+1, Nu);
    calc_argument(root, dis_argument1, dis_argument2, 
                    theta_argument1, theta_argument2);
    for (int i=0; i < MAX_CONTEXT_NUM+1; i++){
        dis[i] = sampler::beta(dis_argument1[i], dis_argument2[i]);
        theta[i] = sampler::gamma(theta_argument1[i], theta_argument2[i]);
    }
}

void VPYLM::update_auxiliary_variables(Node* start_node){
    start_node -> calc_auxiliary_variables(dis, theta);
    for (int i=0; i < start_node -> children.size(); i++){
        Node *next_node = start_node -> children[i];
        if((SumVec(next_node -> a) != 0.0) || (SumVec(next_node -> b) != 0.0)){
            update_auxiliary_variables(next_node);
        }
    }
}

void VPYLM::calc_argument(Node* start_node, vector<double> &dis1,
                          vector<double> &dis2, vector<double> &theta1, 
                          vector<double> &theta2){
    int temp_depth = start_node -> depth;
    if(start_node -> num_all_table >= 2){
        for (int i=0; i < start_node -> num_all_table-1; i++){
            dis1[temp_depth] += 1 - start_node -> sigma[i];
            theta1[temp_depth] += start_node -> sigma[i];
        }
        theta2[temp_depth] -= log(start_node -> rho);
    }
    for (int i=0; i < start_node -> num_vocabulary; i++){
        for (int j=0; j < start_node -> num_all_table; j++){
            if(start_node -> count[i][j] >= 2){
                for (int k=0; k < start_node -> count[i][j]-1; k++){
                    dis2[temp_depth] += 1 - start_node -> tau[i][j][k];
                }
            }
        }
    }
    for (int i=0; i < start_node -> children.size(); i++){
        Node *next_node = start_node -> children[i];
        if((SumVec(next_node -> a) != 0.0) || (SumVec(next_node -> b) != 0.0)){
            calc_argument(next_node, dis1, dis2, theta1, theta2);
        }
    }
}

double VPYLM::SumVec(vector<double> &Vec){
    double sum = 0;
    for (int i=0; i < Vec.size(); i++){
        sum += Vec[i];
    }
    return sum;
}

void VPYLM::Normalize(vector<double> &Vec){
    double sum = 0;
    for (int i=0; i < Vec.size(); i++){
        sum += Vec[i];
    }
    for (int i=0; i < Vec.size(); i++){
        Vec[i] /= sum;
    }
}

void VPYLM::calc_log_likelihood(){
    calc_p_m(root);
    temp_log_likelihood = 0;
    vector<double> temp_ll;
    vector<int> temp_count;
    temp_ll.resize(6,0); temp_count.resize(6,0);
    for (int i=0; i < num_mutation; i++){
        int temp_mutation = mutation[i];
    	temp_ll[temp_mutation] 
            += log(real_customer_arrangement[i] -> p_m[temp_mutation]);
        temp_count[temp_mutation]++;
    }
    for (int i=0; i < 6; i++){
        temp_log_likelihood += temp_ll[i] / temp_count[i];
    }
    vec_log_likelihood.push_back(temp_log_likelihood);
}

void VPYLM::write_data(){
    cout << "------------------" << endl;
    cout << "Writing results..." << endl;
    cout << "------------------" << endl;
    ofstream ofs;
    string output_file_name = "result/" + cancer_type + "_" 
                            +  to_string(exp_index)+ "/corpus_" + cancer_type 
                            + "_" + direction + ".txt";
    ofs.open(output_file_name, ios::out);
    ofs << to_string(temp_log_likelihood) << endl;
    for (int i=0; i < num_mutation; i++){
        ofs << to_string(real_customer_arrangement[i] -> depth) << " ";
    } ofs << endl;
    for (int i=0; i < root -> num_vocabulary; i++){
        for (int j=0; j < 10000; j++){
            string context = posterior_sampling(i);
            ofs << context << " ";
        } ofs << endl;
    }
    for (int i=0; i < MAX_CONTEXT_NUM+1; i++){
        ofs << to_string(dis[i]) << " ";
    } ofs << endl;
    for (int i=0; i < MAX_CONTEXT_NUM+1; i++){
        ofs << to_string(theta[i]) << " ";
    } ofs << endl;
    ofs.close();
    output_file_name = "result/" + cancer_type + "_" + to_string(exp_index)
                     + "/corpus_" + cancer_type + "_" + direction
                     + "_log_likelihood.txt";
    ofs.open(output_file_name, ios::out);
    for (int i=0; i < vec_log_likelihood.size(); i++){
        ofs << to_string(vec_log_likelihood[i]) << endl;
    }
    ofs.close();
}

string VPYLM::posterior_sampling(int mutation){
    Node* temp_node = root;
    string str = to_string(mutation);
    for (int i=0; i < MAX_CONTEXT_NUM+1; i++){
        double stop_probability 
            = (temp_node -> a[mutation] + temp_node -> alpha) 
            / (temp_node -> a[mutation] + temp_node -> alpha 
               + temp_node -> b[mutation] + temp_node -> beta);
        if(sampler::bernoulli(stop_probability) || (i == MAX_CONTEXT_NUM)){ 
            break;
        }
        vector<double> p; p.resize(temp_node -> children.size(), 0);
        int num_temp_mutation = 0;
        for (int j=0; j < temp_node -> children.size(); j++){
            Node* child = temp_node -> children[j];
            for (int k=0; k < child -> num_all_table; k++){
                p[j] += child -> count[mutation][k];
                num_temp_mutation += child -> count[mutation][k];
            }
        }
        if(num_temp_mutation == 0) break;
        Normalize(p);
        int index_child = sampler::categorical(temp_node -> children.size(), p);
        Node* next_node = temp_node -> children[index_child];
        if((SumVec(next_node -> a) == 0.0) && (SumVec(next_node -> b) == 0.0)){
            break;
        }
        int next_base = next_node -> base;
        str += to_string(next_base);
        temp_node = next_node;
    }
    return str;
}

void VPYLM::print_depth_distribution(){
    cout << "depth distribution" << endl;
    vector<int> num_depth;
    num_depth.resize(MAX_CONTEXT_NUM+1, 0);
    for (int i=0; i < num_mutation; i++){
        num_depth[real_customer_arrangement[i] -> depth]++;
    }
    for (int i=0; i < MAX_CONTEXT_NUM+1; i++){
        cout << to_string(i) + " : " + to_string(num_depth[i]) << endl;
    }
}

void VPYLM::print_root_count(int iter, int index, string status){
    cout << "iter : " << to_string(iter) << ", index : " << to_string(index) 
         << ", mutation : " << to_string(mutation[index]) << endl;
    cout << "status : " << status << endl;
    for (int i=0; i < root -> num_vocabulary; i++){
        for (int j=0; j < root -> num_all_table; j++){
            cout << to_string(root -> count[i][j]) << " ";
        }
        cout << endl;
    }
}

