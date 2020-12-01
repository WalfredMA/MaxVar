//
//  probability_calculator.hpp
//  snp_caller
//
//  Created by Wangfei MA on 7/8/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef probability_calculator_hpp
#define probability_calculator_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <cstring>
#include "allele.hpp"
#include "mutation.hpp"

using namespace std;


class prob_calculator{
    
public:
    
    prob_calculator(){};
    
    void loadStat(string statfile);
    
    void loadMole(string molefile);
    
    void loadSNP(string snpfile);
    
    float systemfilter(int first_alleleindex, int second_alleleindex, int first_obs, int second_obs, int &f_freq_cut, int &r_freq_cut, int &f_hcutoff, int &r_hcutoff);
    
    int BaseQualfilter(int first_alleleindex, int second_alleleindex, float first_qual, float second_qual);
        
    float AlleleHaplo(int allele_index, int *allele_count, float *posi_label, float *nega_label);
    
    int BaseHaplo(vector<pair<vector<int>&,vector<int>&>> &all_moles, float *false_n_1, float *false_n_2);
    
    float cul_bionom(int total, int obs);
    
    float RandomError(int total, int obs);
    
    float cul_bionom_error(int total, int obs);
    
    int valide_start, valide_end, valide_size;
    
    std::map<int, std::vector<int>> mole_phase;
    
    std::set<int> refsnp, lowsnp, used_group;
    
    std::map<int, std::vector<char>> labels;
    
    mutation* mutations[1024];
    
    allele* alleles[8];  
    
private:
    
    float permutation_log(int start,int end);
    
    float factorial_log(int end);
    
    double bionomial(int total, int obs);
    
    double bionomial_error(int total, int obs);
    
    float log_add(float x, float y);
    
    bool** save_results[64];
    bool** ifsave_results[64];
    
};


#endif /* probability_calculator_hpp */
