//
//  mutations.cpp
//  snp_caller
//
//  Created by Wangfei MA on 7/22/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "mutation.hpp"

extern float median_coverage;
extern float mean_coverage;

float mutation::factorial_log(int end)
{
    
    static float results_storage[1000];
    static bool results_index[1000]={0};
    
    if (end<1000 && results_index[end])
    {
        
        return results_storage[end];
        
    }
    
    float factorial=0.0;
    
    for (int i=1;i<=end;i++)
    {
        factorial+=log(i);
    }
    
    
    if (end<1000)
    {
        
        results_index[end]=1;
        results_storage[end]= factorial;
        
    }
    
    return  factorial;
}

float mutation::permutation_log(int start,int end)
{
    
    return factorial_log(end)-factorial_log(start);
}



float mutation::error_confi(int corr,int error)
{
    
    if (error<((int)ifsave.size()) && corr<(int)ifsave[error]->size() && (*ifsave[error])[corr])
    {
        
        return (*count_prob[error])[corr];
    }

    int push_size=(error-(int)count_prob.size()+1);
    
    if (push_size)
    {
        
        for (int i=0;i<push_size;i++)
        {
            
            ifsave.push_back(new vector<bool>);
            count_prob.push_back(new vector<float>);
        }
    }
    
    push_size=(int)(corr-(int)(*count_prob[error]).size()+1);
    if (push_size)
    {
        
        for (int i=0;i<push_size;i++)
        {
            
            (*ifsave[error]).push_back(0);
            (*count_prob[error]).push_back(0.0);
        }
    }
    
    float confi=0.0;
    if (error_coefs[0])
    {
        
        float confi_from_error=error_coefs[0]*error*error+error_coefs[1]*error+error_coefs[2];
        
        float slope_corr=corr_coefs[0]*error*error+corr_coefs[1]*error+corr_coefs[2];
        
        confi_from_error-=permutation_log(1,error);
        
        float confi_from_corr=(median_coverage-corr)*slope_corr;
        
        confi=-confi_from_error-confi_from_corr;
    }
    
    (*ifsave[error])[corr]=confi;
    (*ifsave[error])[corr]=1;
    
    
    return confi;
    
}
