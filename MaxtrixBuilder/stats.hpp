//
//  stats.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef stats_hpp
#define stats_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <math.h>

#define exp_pro 0.90
#define unex_pro 0.10

//log2(0.9)
#define exp_pro_log -0.152

//log2(0.1)
#define unex_pro_log -3.322

class stats{
    
public:
    
    float Cul_Bionom(int total, int obs);
    
    double Bionomial(int total, int obs);
    
    float Factorial_log(int start,int end);
        

private:
    
    float results_storage[100][100];
    float results_history[100][100];
    
};


#endif /* stats_hpp */
