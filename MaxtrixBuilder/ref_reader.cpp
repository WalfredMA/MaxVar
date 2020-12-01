//
//  ref_reader.cpp
//  testcpp
//
//  Created by Wangfei MA on 8/13/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "ref_reader.hpp"
#include <regex>


#define mask_exd 1

using namespace std;

int ref_Reader::Masker(string& seq,int region_start,int region_end){
    
    int len=(int)seq.length();
    
    region_start=max(region_start-10,1);
    region_end=min(region_end+10,len);
    
    bool* ifmask=(bool*)malloc((region_end-region_start+2*mask_exd)*sizeof(bool));
    memset(ifmask, 1, (region_end-region_start+2*mask_exd)*sizeof(bool));
    
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);// unmask
    
    regex rx("([^N]{1,4})\\1(\\w{0,1}\\1){2,}", std::regex_constants::icase);
    
    smatch forward_search;
    
    int allowgap;
    int left_pos,right_pos;
    string find_pattern;
    size_t forward_find;
    smatch m;
    
    int total_mask=0,num_mask=0;
    
    for (sregex_iterator i = sregex_iterator(seq.begin()+region_start+mask_exd, seq.begin()+region_end-mask_exd,rx); i != sregex_iterator(); ++i)
    {
        m = *i;
        find_pattern=m[1];
        allowgap=1+(int)m[1].length();
        
        right_pos=min(len,(int)(m.position()+m.length()));
        left_pos=(int)m.position();
        forward_find=0;
        
        while (left_pos>allowgap && (forward_find=seq.substr(left_pos-allowgap, allowgap).find(find_pattern))!=string::npos){
            
            left_pos+=(int)forward_find-allowgap;
        }
        
        total_mask+=right_pos-left_pos+2*mask_exd;
        ++num_mask;
        
        memset(ifmask+left_pos-mask_exd,false,(right_pos-left_pos+2*mask_exd)*sizeof(bool));
    }
        
    for (int i=0;i< region_end-region_start; ++i){
        
        if(!ifmask[i]) {seq.at(i+region_start)=tolower(seq[i+region_start]);}
    }
    
    free(ifmask);
    
    return total_mask;
}

string ref_Reader::Load_Chr(const char* chr,int start,int end){
    
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
    
    printf("Start masking high error-prone sequences and motifs on %s\n", chr);
    
    //printf("Finished masking: %d sites\n", Masker(reference, start, end));
        
    return reference;
    
}

