//
//  MLE.hpp
//  sample_analyzer
//
//  Created by Wangfei MA on 7/17/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef MLE_hpp
#define MLE_hpp


#include <iostream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

class MLE{
    
public:
    
    MLE(){};
        
    bool log_Gaussian(vector<double> &frequencies, int binstart, float *coef_a, float *coef_b, float *coef_c);
    
    bool Gaussian(vector<double> &frequencies, int binstart, float *err, float *mean);
    
    bool Poisson(vector<double> &frequencies, int binstart, float *lambda);
    
    float quadratic_fitting(const vector<float>& x, const vector<float>& y, float *coef_a, float *coef_b, float *coef_c);
    
    float linear_fitting(const vector<float>& x, const vector<float>& y, float *slope, float *intercept);
        
    float slope(const vector<float>& x, const vector<float>& y, float *slope );

private:
    
    float permutation_log(int start,int end);
    
    float factorial_log(int end);
        
};

#endif /* MLE_hpp */
