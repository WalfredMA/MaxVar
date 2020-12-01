//
//  alleles.cpp
//  snp_caller
//
//  Created by Wangfei MA on 7/22/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "allele.hpp"


float allele::log_add(float x, float y){
    
    
    float higher=(x>y)?x:y;
    
    float diff=fabs(x-y);
    
    if (diff>3){
        
        return higher;
        
    }
    
    return higher+log(1+exp(-diff));
    
}

float allele::factorial_log(int end){
    
    static float results_storage[100];
    static bool results_index[100]={0};
    
    if (end<100 && results_index[end]){
        
        return results_storage[end];
        
    }
    
    float factorial=0.0;
    
    for (int i=1;i<=end;i++){
        factorial+=log(i);
    }
    
    
    if (end<100){
        
        results_index[end]=1;
        results_storage[end]= factorial;
        
    }
    
    return  factorial;
}

float allele::permutation_log(int start,int end){
    
    return factorial_log(end)-factorial_log(start);
}



float allele::coverage_confi(int obs, int ifhomo){
    
    if (obs<ifsave[ifhomo].size() && ifsave[ifhomo][obs]){
        
        return count_prob[ifhomo][obs];
    }
    
    int push_size=(int)(obs-count_prob[ifhomo].size()+1);
    if (push_size){
        
        for (int i=0;i<push_size;i++){
            
            ifsave[ifhomo].push_back(0);
            count_prob[ifhomo].push_back(0.0);
        }
    }
    
    float confi=0.0;
    if (ifhomo){
        
        confi=log_add(obs_prob(obs, 2.0),obs_prob(obs, 3.0));
    
    }
    
    else {
        
        confi=obs_prob(obs, 1.0);
        
    }
    
    
    count_prob[ifhomo][obs]=confi;
    ifsave[ifhomo][obs]=1;
    
    return confi;
    
}


float allele::obs_prob(int obs, float assumcopy){
    
    float result=0.0;
    
    float expected=coverage_exp*assumcopy/2;
    float errorlog=coverage_errorlog+log(assumcopy/2)/2;
    
    if (ifpoisson==1) {
        
        for (int i=0; i<obs; i++){
            
            result+=exp(i*2*errorlog-expected-permutation_log(1,i));
        }
        result=log(result);
        
    }
    
    else {
        
        float error=exp(errorlog);
        float discre=(expected<obs)?0:(expected-obs);
        result= log(erfc(discre/error)/2);
    }
    
    return result;
}

