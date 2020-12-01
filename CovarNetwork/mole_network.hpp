//
//  mole_network.hpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef mole_network_hpp
#define mole_network_hpp

#include <iostream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <tuple>

using namespace std;

#define hash0 98
#define hash1 99
#define maxvalue 20
#define mole_combine 0

class Mole_network
{
    
public:
    
    Mole_network();
    
    int findindex(int barcode);
    
    void load_mole(int barcode, vector<int> &coordi_p, vector<int> &coordi_n);
    
    void set_hash(vector<int> &vecx_p, vector<int> &vecx_n, vector<bool> &badbar);
    
    void itsect(vector<int> &vecy_p, vector<int> &vecy_n, vector<bool> &badbar, int *posi_counter, int *nega_counter);
    
    void covar(float covar_cutoff);
    
    vector<int> mole_barcode;
    
    vector<vector<int>> mole_matrix_p;
    
    vector<vector<int>> mole_matrix_n;
    
    vector<vector<int>> mole_mole;
    
    vector<bool> *badbar;
    
    vector<bool> badindex;
    
    vector<int> covar_num;
    vector<int> covar_dom;
    
    int mole_counter, current_mole_index, current_mole_coordi;
    
    
    
private:
    
    vector<int> hash_index[100][100];
    vector<int> hash_coordi[100][100];
    
    vector<int> hash_dic[40];
    vector<bool> hash_dic_sign[40];
    
    
};



#endif /* mole_network_hpp */
