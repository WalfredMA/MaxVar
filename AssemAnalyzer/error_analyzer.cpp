//
//  stochastic_analyzer.cpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/27/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "error_analyzer.hpp"


using namespace std;

float Error_anal::estimate_error(vector<double> &frequencies, int input_coverage){
    
    double total_error=1;
    for (int i=0;i<input_coverage/4;i++)
    {
        
        total_error+=frequencies[i]*i;
    }
    
    total_error/=accumulate(frequencies.begin(), frequencies.end(), 1);
    total_error/=input_coverage;
    
    return total_error;
    
}


bool Error_anal::modeling(vector<double> &frequencies, int input_coverage, error& type_error)
{
    
    int lowconfi_cutoff=(int)frequencies.size();
    
    int highconfi_cutoff=(int)frequencies.size();
    
    int first_decrea=0, first_increa=0;
    
    if (find_Inflection(frequencies, &first_decrea, &first_increa)<highconfi_cutoff+1)
    {
        int fitting_start=first_decrea+1,fitting_end=first_decrea+5;
        
        fitting_end=(fitting_end>highconfi_cutoff)?highconfi_cutoff:fitting_end;
        
        vector<double> fitting_bins (frequencies.begin()+fitting_start, frequencies.begin()+fitting_end);
        
        if ((fitting_end-fitting_start)>2 && TailFitting(fitting_bins,fitting_start, type_error))
        {
            
            float error_slope=type_error.curve_quadratic[0];
            
            int iter=fitting_end-1;
            
            lowconfi_cutoff=iter;
            
            float proj_error=log(frequencies[iter]+1);
        
            for (;iter<highconfi_cutoff;iter++)
            {
                
                if (log(frequencies[iter]+1)-proj_error>log(10))
                {
                    break;
                }
                
                proj_error+=error_slope;
            }
            
            highconfi_cutoff=iter;
        }
    }
    
    double reverse_cdf=accumulate(frequencies.begin(), frequencies.end(), 0);
    double cutnumber=max(10,(int)reverse_cdf/100);
    
    int max_cutoff=0;
    int min_cutoff=0;
    int lowconfi_max_cutoff=0;
    
    for (int iter2=0;iter2<highconfi_cutoff;iter2++)
    {
        
        reverse_cdf-=frequencies[iter2];
        
        if (!min_cutoff && reverse_cdf<=2*cutnumber)
        { //the cutoff must filter 99% of the sites.
            min_cutoff=iter2;
        }
        
        if (!lowconfi_max_cutoff && reverse_cdf<=cutnumber)
        {//the cutoff reserve 0.5% of the sites.
            lowconfi_max_cutoff=iter2;
        }
        
        if (reverse_cdf<=cutnumber/10)
        { //it must reserve 0.05% of the sites
            max_cutoff=iter2;
            break;
        }
    }
    
    type_error.lowconfi_cut=max(min(lowconfi_cutoff,lowconfi_max_cutoff),min_cutoff);
    
    type_error.errorcut=max(min(highconfi_cutoff,max_cutoff),lowconfi_max_cutoff);
        
    return 1;
    
}

bool Error_anal::TailFitting(vector<double> &frequencies, int binstart, error& type_error)
{
    
    float coef_a, coef_b, coef_c;
    if (this->log_Gaussian(frequencies, binstart, &coef_a, &coef_b, &coef_c))
    {
        
        type_error.curve_quadratic=vector<float>{coef_a,coef_b,coef_c};

        return 1;
        
    }
    
    
    else
    {
        
        return 0;
        
    }
    
}


int Error_anal::find_Inflection(vector<double> &frequencies,int *first_decrea, int* first_increa){
    
    int binsize=(int)frequencies.size();
    if (binsize<2) return binsize;
    
    double lastpoint=frequencies[0];
    
    int i=1;
    for (; i<binsize; i++){
        
        if (frequencies[i]<lastpoint){
            *first_decrea=i-1;
            lastpoint=frequencies[i];
            break;
        }
        
        lastpoint=frequencies[i];
        
    }
    
    for (; i<binsize; i++){
        
        if (frequencies[i]>lastpoint){
            *first_increa=i-1;
            
            lastpoint=frequencies[i];
            break;
        }
        
        lastpoint=frequencies[i];
    }
    
    
    return i;
    
}
