//
//  snp_records.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef snp_record_hpp
#define snp_record_hpp

#include "Config.h"
#include "stats.hpp"
#include "Cigar.hpp"
#include "MDTag.hpp"
#include "realign.hpp"
#include "Base_buffer.hpp"

#include <time.h>

class Snp_Records
{
    
public:
    
    Snp_Records(const char* file);
    ~Snp_Records(){free(writeLine);}
    
    void set_Region(int phase_start);
    
    bool WriteSnp(int write_end);
    
    int write_start, write_end;
    
    const char* outfile;
    std::string samplefile;
            
private:
    
    void Count_Allele(int* moles, char *allquals, int* all_mole, int* goold_mole, int* qual, int len);
    
    void Count_SmallDel(Base* current_base, int* all_mole, int* goold_mole, int* qual, int len);
    
    bool Alelle_Stat(int current_index,FILE *samplefile,int notsamplex);
    
    int Load_Line(int current_index);
        
    void Correct_SNPs(int write_end);
    
    Realigner aligner;
    
    stats calculator;
    
    char *writeLine=(char*)malloc(1000000*sizeof(char)); 
    
    int Linesize;

    
        
};



#endif /* snp_record_hpp */
