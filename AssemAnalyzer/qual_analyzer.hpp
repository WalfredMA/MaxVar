//
//  qual_analyzer.hpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/27/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef qual_analyzer_hpp
#define qual_analyzer_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <map>

#define  Vprint( n )  for ( auto i: n ) { cout<<i<<" ";} cout<<endl;


#include "MLE.hpp"

class Qual_anal: MLE{
    
public:
    
    Qual_anal(){};
    
    bool modeling(std::vector<int> frequencies, int& qualcut, int& qualcut_low);
    
};

#endif /* qual_analyzer_hpp */
