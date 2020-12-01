//
//  coverage_analyzer.hpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/27/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef coverage_analyzer_hpp
#define coverage_analyzer_hpp

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

typedef struct
{
    std::string type;
    float mean;
    float std_log;
    
}cover_curve;


class Cover_anal: MLE{
    
public:
    
    Cover_anal(){};
    
    bool modeling(std::vector<double> frequencies, int binstart, cover_curve& allele_cover);
            
    
private:
    
    bool PeakFitting(std::vector<double> &frequencies, int binstart , cover_curve&);
        
};


#endif /* coverage_analyzer_hpp */
