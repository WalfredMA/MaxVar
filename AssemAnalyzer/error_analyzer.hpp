//
//  stochastic_analyzer.hpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/27/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef stochastic_analyzer_hpp
#define stochastic_analyzer_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <map>

#include "MLE.hpp"

#define  Vprint( n )  for ( auto i: n ) { cout<<i<<" ";} cout<<endl;


typedef struct
{
    
    vector<float> curve_quadratic{0.0,0.0,0.0};
    
    double errorrate=0.0;
    
    int errorcut=-1;
    
    int lowconfi_cut=-1;
    
}error;


class Error_anal: MLE{
    
public:
    
    Error_anal(){};
    
    bool modeling(vector<double> &frequencies, int binstart, error& type_error);
    
    float estimate_error(vector<double> &frequencies, int input_coverage);
    
    
private:
    
    bool TailFitting(vector<double> &frequencies, int binstart, error& type_error);
    
    int find_Inflection(vector<double> &frequencies,int *first_decrea, int* first_increa);
    
};


#endif /* stochastic_analyzer_hpp */
