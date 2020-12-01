//
//  stats.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "stats.hpp"

using namespace std;


float stats::Factorial_log(int start,int end){
    
    float factorial=0.0;
    
    for (int i= start;i<=end;++i){
        factorial+=log2f(i);
    }
    return factorial;
}

double stats::Bionomial(int total, int obs){
    
    int notobs=total-obs;
    
    float result= obs* exp_pro_log+notobs* unex_pro_log +Factorial_log(total-obs+1, total) - Factorial_log(2, obs);
    
    return result;
}


float stats::Cul_Bionom(int total, int obs){
    
    
    if (total-obs>=10 && total<100){
        return 10.0;
    }
    
    else if (!total || !obs || !(total-obs)  || obs>total || total>=1000) {
        return 0.0;
    }
    
    if (total<100 && obs<100 &&results_history[total][obs]){
        
        return results_storage[total][obs];
        
    }
    
    double culmi=0.0;
    for (int i=total-obs;i<=total;++i){
        
        culmi+=pow(2,Bionomial(total, i));
    }
    
    float result=-log2f(culmi);
    
    if (total<100 && obs<100){
        results_storage[total][obs]=result;
        results_history[total][obs]=1;
    }
    return result;
}
