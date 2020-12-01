//
//  networks.hpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef networks_hpp
#define networks_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <math.h>
#include <algorithm>
#include <tuple>
#include <map>

#include "mole_network.hpp"
#include "snp_network.hpp"

using namespace std;

#define hash0 98
#define hash1 99

class covar_networks
{
    
public:
    
    covar_networks(): mole_network(), snp_network() {}; 
    
    int read_line(string StrLine, vector<int> &p_col, vector<int> &n_col);
    
    void load_data(string datafile);
        
    void connect_snp();
    
    void connect_mole();
    
    void init(string datafile);
    
    void run(vector<float> cutoffs);
    
    void pre_filtration();
    
    void output(string);
    
    map<int, std::string> labels;
    
    Mole_network mole_network;
    Snp_network snp_network;

private:
    
    bool hash_bool[200];
    
    void cycle(float cutoff);
    
    float factorial_log(int end);
    
    double bionomial(int total, int obs);
    
    float cul_bionom(int total, int obs);
    
    float bionomial_error(int total, int obs);
    
    float cul_bionom_error(int total, int obs);
    
    float confidence(int total, int obs);
    
    void sv_innerfilter(vector<tuple<int, vector<int>, vector<int>>>& );
};



#endif /* networks_hpp */
