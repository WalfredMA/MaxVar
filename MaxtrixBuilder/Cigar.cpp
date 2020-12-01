//
//  cigar_Reader.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "Cigar.hpp"

using namespace std;

int Cigar_Reader::Load(const char *cigar_text){ 
    
    int cigarlen=(int )strlen(cigar_text);
    
    char current_type;
    
    int current_segment=0;
    
    int alignment_size=0;
    
    cigar_counter=0;
    
    for (int i = 0; i < cigarlen; ++i){
        
        if (cigar_mallocsize<=cigar_counter){
            
            cigar_segments.push_back(0);
            
            cigar_types.push_back('A');
            
            ++cigar_mallocsize;
            
        }
        
        current_type=cigar_text[i];
        
        if ((current_type>='a' && current_type<='z') || (current_type>='A' && current_type<='Z')){
            
            if (current_type=='I'){
                
                if (!cigar_counter) {
                    
                    cigar_segments[cigar_counter]=current_segment;
                    cigar_types[cigar_counter]='S';
                    ++cigar_counter;
                    
                }
                
                else {
                    
                    cigar_segments[cigar_counter]=current_segment+1;
                    cigar_types[cigar_counter]='I';
                    cigar_segments[cigar_counter-1]-=1;
                    ++cigar_counter;
                    ++alignment_size;
                    
                }
                current_segment=0;
                
            }
            
            else{
                
                cigar_segments[cigar_counter]=current_segment;
                cigar_types[cigar_counter]=current_type;
                ++cigar_counter;
            
                if (current_type!='S' &&current_type!='N'){
                    alignment_size+=current_segment;
                }
                
                current_segment=0;
                
            }
        }
        
        else if (current_type>='0'&&current_type<='9'){
            
            current_segment*=10;
            current_segment+=(int)(current_type-'0');
            
        }
        
        else {
            
            break;
        }
        
    }
    
    return alignment_size;
    
}


int Cigar_Reader::Load(std::vector<BamTools::CigarOp>* cigar_data) {
        
    char current_type;
    
    int current_segment;
    
    int alignment_size=0;
    
    cigar_counter=0;
        
    for (size_t i=0; i<cigar_data->size(); ++i){
        
        if (cigar_mallocsize<=cigar_counter){
            
            cigar_segments.push_back(0);
            
            cigar_types.push_back('A');
            
            ++cigar_mallocsize;
            
        }
                
        current_type=(*cigar_data)[i].Type;
        
        current_segment= (int) (*cigar_data)[i].Length;
        
        if ((current_type>='a' && current_type<='z') || (current_type>='A' && current_type<='Z')){
            
            if (current_type=='I'){
                
                if (!cigar_counter || cigar_types[cigar_counter-1]!='M') {
                    
                    cigar_segments[cigar_counter]=current_segment;
                    cigar_types[cigar_counter]='S';
                    
                }
                
                else {
                    
                    cigar_segments[cigar_counter]=current_segment+1;
                    cigar_types[cigar_counter]='I';
                    cigar_segments[cigar_counter-1]-=1;
                    ++cigar_counter;
                    
                }
                
                current_segment=0;
                
            }
            
        else{
                
            cigar_segments[cigar_counter]=current_segment;
            cigar_types[cigar_counter]=current_type;
            ++cigar_counter;
            
            if (current_type!='H' && current_type!='S'){
                alignment_size+=current_segment;
            }
            current_segment=0;
                
            }
        }
        
    }
    
    return alignment_size;
    
}
