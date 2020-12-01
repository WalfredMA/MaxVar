//
//  read_reader.cpp
//  sample_analyzer
//
//  Created by Wangfei MA on 8/30/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "ref_reader.hpp"

using namespace std;
//reading reference
string ref_Reader::Load_Chr(const char *inputfile, const char* chr){ 
    
    string StrLine, reference;
    
    if (strlen(inputfile)<2) return reference;
    
    fstream fafile;
    fafile.open(inputfile, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << inputfile << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    while (getline(fafile,StrLine)){
        
        if (StrLine[0]!='>' || StrLine.compare(1, min(StrLine.find(' '),StrLine.find('\t')), string(chr))) continue;
        break;
    }
    
    printf("Start reading chromosome: %s \n",chr);
    
    while (getline(fafile,StrLine)){
        
        if (StrLine[0]=='>' || StrLine[0]==' ') break;
        
        reference+=StrLine;
        
    }
    
    fafile.close();
    
        
    return reference;
    
}


