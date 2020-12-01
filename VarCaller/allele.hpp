//
//  alleles.hpp
//  snp_caller
//
//  Created by Wangfei MA on 7/22/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef allele_hpp
#define allele_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

class allele{

public:
    
    allele(){};
    
    float coverage_confi(int obs, int ifhomo);
    
    const char* type="default";
    
    string name;
    
    int index;
    
    bool ifpoisson=1;
    
    float coverage_exp=0.0,coverage_errorlog=0.0;
    
    float qualcut=0.0, lowqual=0.0;

private:
    
    float obs_prob(int obs, float assumcopy);
    
    float allcopy_prob(int obs);
    
    vector<float> count_prob[2];
    vector<bool> ifsave[2];
    
    float permutation_log(int start,int end);
    
    float factorial_log(int end);
    
    float log_add(float x, float y);

};

#endif /* allele_hpp */
