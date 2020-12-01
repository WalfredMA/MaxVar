//
//  ref_reader.hpp
//  testcpp
//
//  Created by Wangfei MA on 8/13/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef ref_reader_hpp
#define ref_reader_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>

class refreader{
    
public:
    
    refreader(const char *file):inputfile(file){};
    
    std::string load_chr(const char* chr);
    
    const char *inputfile;
    
};

#endif /* ref_reader_hpp */
