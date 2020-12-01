//
//  mutation.hpp
//  snp_caller
//
//  Created by Wangfei MA on 7/22/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef mutation_hpp
#define mutation_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>

extern int confidence;

using namespace std;

class mutation{
    
public:
    
    mutation(){};
    
    float error_confi(int corr,int error);
    
    const char* type="default";
    
    string name;
    
    int index;
    
    int errorcut=19/confidence;
    
    int highcut=29/confidence;
    
    float error_curve[2]={0};
    
    float error_coefs[3]={0};
    
    float corr_coefs[3]={0};
    
private:
    
    vector<vector<float>*> count_prob;
    
    vector<vector<bool>*> ifsave;
    
    float permutation_log(int start,int end);
    
    float factorial_log(int end);
    
};

#endif /* mutation_hpp */
