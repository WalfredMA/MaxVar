//
//  Bamloader.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/26/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef Bamloader_hpp
#define Bamloader_hpp

#include <iostream>
#include <stdio.h>

#include "api/BamReader.h"
#include "Config.h"
#include "Fileloader.hpp"


class Bamloader: public Fileloader
{
    
public:
    
    Bamloader(const char* input, const char* output, const char* ref):Fileloader(input, output, ref){};
    
    ~Bamloader(){};
private:
    
    void Locate_region();
    
    void Start_Read();
    
    int GetNextAlignment(BamTools::BamAlignment &alignment, Cigar_Reader &cigar, int *barindex);
    
    void Close();
    
    BamTools::BamReader Baminfile;
    
};

#endif /* Samloader_hpp */

