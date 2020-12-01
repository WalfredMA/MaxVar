//
//  coverage_analyzer.cpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/27/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "coverage_analyzer.hpp"

using namespace std;

bool Cover_anal::modeling(vector<double> frequencies, int binstart, cover_curve& allele_cover){
    
    size_t binsize=frequencies.size();
    size_t left_end=binstart, right_end=binstart+binsize-1;
    
    long long highest_count=0;
    int highest_index=binstart;
    
    for (int i=binstart; i<binstart+binsize; i++){
        
        if (frequencies[i-binstart]>highest_count){
            
            highest_count=frequencies[i-binstart];
            highest_index=i;
        }
    }
    
    
    for (int i=highest_index; i>0; i--){
        
        if (frequencies[i-binstart-1]<10 ||  frequencies[i-binstart-1]>=frequencies[i-binstart]) {
            
            left_end=i;
            break;
        }
    }
    
    for (int i=highest_index; i<binstart+binsize-1; i++){
        
        if (frequencies[i-binstart+1]<10 ||
            frequencies[i-binstart+1]>= frequencies[i-binstart]) {
            
            right_end=i;
            break;
        }
    }
    
    size_t fitsize=highest_index/4;
    
    left_end=(left_end>highest_index-fitsize)?left_end:(highest_index-fitsize);
    
    right_end=(right_end<highest_index+fitsize)?right_end:(highest_index+fitsize);
    
    vector<double> fitting_bins (frequencies.begin()+left_end-binstart, frequencies.begin()+right_end+1-binstart);
    
    if (right_end-left_end<3)return 0;
    
    return PeakFitting(fitting_bins, (int)left_end, allele_cover);
    
}


bool Cover_anal::PeakFitting(vector<double> &frequencies, int binstart, cover_curve& current_model){
    
    float tempresult1, tempresult2;
    if (this->Poisson(frequencies, binstart,&tempresult1)){
        
        current_model.type="Poisson";
        
        current_model.mean=tempresult1;
        
        current_model.std_log=log(tempresult1)/2;
        
        return 1;
        
    }
    
    else if (this->Gaussian(frequencies, binstart, &tempresult1, &tempresult2)){
        
        current_model.type="Gaussian";
        
        current_model.mean=tempresult1;
        
        current_model.std_log=log(tempresult2);
        
        return 1;
        
    }
    
    else{
        
        return 0;
        
    }
    
}
