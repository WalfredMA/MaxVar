//
//  MDTag.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/28/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef MDTag_hpp
#define MDTag_hpp

#include "Config.h"
#include <stdio.h>
#include <vector>

class MDTag{
    
public:
    
    MDTag(){};
    ~MDTag(){};
    
    bool GrepMD(std::string& des, const char *md_text);
    
};


#endif /* MDTag_hpp */
