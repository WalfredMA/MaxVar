//
//  Samloader.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/26/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef Samloader_hpp
#define Samloader_hpp

#include <stdio.h>
#include <fstream>

#include "Config.h"
#include "Fileloader.hpp"


class Samloader: public Fileloader 

{
    
public:
    
    Samloader(const char* input, const char* output, const char* ref):Fileloader(input, output, ref){};
    
    ~Samloader(){};
private:
    
    void Locate_region();
    
    void Start_Read();
    
    int GetNextAlignment(BamTools::BamAlignment &alignment, Cigar_Reader &cigar, int *barindex);
    
    void Close();
    
    long long Grep_Barcode(const char*);
    
    std::fstream Saminfile;
    
};

#endif /* Samloader_hpp */
