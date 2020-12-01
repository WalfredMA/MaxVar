//
//  qual_analyzer.cpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/27/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "qual_analyzer.hpp"


bool Qual_anal::modeling(vector<int> frequencies,int& qualcut, int& qualcut_low)
{
    int cutoff=0;
    
    double total=accumulate(frequencies.begin(), frequencies.end(), 0LL);
    double thesum=0;
    
    
    qualcut_low=0;
    for (;cutoff<frequencies.size();++cutoff)
    {
        
        thesum+=frequencies[cutoff];
        
        if (!qualcut_low && thesum>=total/500)
        {
            qualcut_low=cutoff-1;
            
        }
        
        if (thesum>=total/100) break;
        
    }
    
    qualcut=cutoff;
    
    return 1;
}

