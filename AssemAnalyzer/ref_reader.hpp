//
//  ref_reader.hpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/30/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef ref_reader_hpp
#define ref_reader_hpp

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>

class ref_Reader{
    
public:
    
    ref_Reader(const char *file):inputfile(file){};
    ~ref_Reader(){};
    
    std::string Load_Chr(const char* inputfile, const char* chr);
        
    const char *inputfile;
    
};

#endif /* ref_reader_hpp */
