//
//  Cigar.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef cigar_hpp
#define cigar_hpp

#include <stdio.h>
#include <vector>

#include "api/BamReader.h"

class Cigar_Reader{
    
public:
    
    Cigar_Reader(): cigar_counter(0), cigar_mallocsize(0) {}
    ~Cigar_Reader(){};
    
    int Load(std::vector<BamTools::CigarOp>* cigar_data);
    
    int Load(const char *cigar_text);
    
    int cigar_counter;
    
    std::vector<int> cigar_segments;
    
    std::vector<char> cigar_types;
    
protected:
    
    std::vector<BamTools::CigarOp>* cigarbam;
    
    const char* cigarText;
    
    char type;
    
    int cigar_mallocsize;
    
};

#endif /* cigar_hpp */
