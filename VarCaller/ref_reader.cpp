//
//  ref_reader.cpp
//  testcpp
//
//  Created by Wangfei MA on 8/13/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "ref_reader.hpp"

using namespace std;

string refreader::load_chr(const char* chr){
    
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
    
    while (getline(fafile,StrLine)){
        
        if (StrLine[0]=='>' || StrLine[0]==' ') break;
        
        reference+=StrLine;
        
    }
    
    return reference;
    
}

