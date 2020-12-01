//
//  snps_network.hpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef snp_network_hpp
#define snp_network_hpp

#include <iostream>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <utility>
#include <algorithm>
#include <tuple>

using namespace std;

#define hash0 98
#define hash1 99
#define maxvalue 20
#define combine_dis 5


class Snp_network
{
    
public:
    
    Snp_network();
    
    
    int findindex(int coordi);
        
    int load_snp(int coordi, vector<int> &bar_p, vector<int> &bar_n);
    
    void set_hash(vector<int> &vecx_p, vector<int> &vecx_n, vector<bool> &badbar);
    
    void itsect(vector<int> &vecy_p, vector<int> &vecy_n, vector<bool> &badbar, int *posi_counter, int *nega_counter);
    
    void covar(float covar_cutoff);
    
    vector<int> snp_coordi;
    
    vector<vector<int>> snp_matrix_p; 
    
    vector<vector<int>> snp_matrix_n;
    
    vector<vector<int>> snp_snp;
    
    vector<bool> *badbar;
    
    vector<bool> badindex;
    
    vector<int> covar_num;
    vector<int> covar_dom;
    
    int snp_counter;
    
    
private:
    
    vector<int> hash_index[100][100];
    vector<int> hash_coordi[100][100];
    
    vector<int> hash_dic[40];
    vector<bool> hash_dic_sign[40];
    
    void find_uniq(vector<int>& list);
};



#endif /* snp_network_hpp */
