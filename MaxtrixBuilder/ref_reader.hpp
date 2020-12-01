//
//  ref_reader.hpp
//  testcpp
//
//  Created by Wangfei MA on 8/13/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef ref_reader_hpp
#define ref_reader_hpp

#include "Config.h"

#include <stdio.h>


class ref_Reader{
    
public:
    
    ref_Reader(const char *file):inputfile(file){};
    ~ref_Reader(){};
    
    std::string Load_Chr(const char* chr ,int start,int end);
    
    int Masker(std::string& seq,int start,int end);
        
    const char *inputfile;
        
};

#endif /* ref_reader_hpp */
