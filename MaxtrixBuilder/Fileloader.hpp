//
//  Fileloader.hpp
//  testcpp
//
//  Created by Wangfei MA on 8/19/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef Fileloader_hpp
#define Fileloader_hpp

#include <iostream>
#include <stdio.h>

#include "api/BamReader.h"
#include "Config.h"
#include "Reads_Processor.hpp"
#include "barcodes_record.hpp"
#include "snp_records.hpp"
#include "ref_reader.hpp"
#include "MDTag.hpp"
#include "Cigar.hpp"


class Fileloader{
    
public:
    
    Fileloader(const char* input, const char* output, const char* ref): Barcodes_Records(output), snps_records(output), reader(ref), line_index(0), inputfile(input), outputfile(output), reffile(ref), baifilename( (std::string(inputfile)+".bai")) ,last_write_start(0),current_start(0) {};
    
    ~Fileloader(){}; 
    
    bool Process_region(const char *, int , int);
    
    int LoadAlignment(BamTools::BamAlignment &BamLine, Cigar_Reader &cigar, int* barindex, int pid);
    
    int phase_end, phase_start; 
    
    std::string phase_chrID;
    
    int last_write_start,current_start;
    
protected:
    
    virtual void Locate_region() {};
    
    virtual void Start_Read() {};
    
    virtual void Close() {};
        
    virtual int GetNextAlignment(BamTools::BamAlignment &alignment, Cigar_Reader &cigar, int* barindex) {return 0;};
        
    int Filtration(std::string, int Position, int MapQuality, int AlignmentFlag);
        
    bool Wait_Write(int);
        
    bool Write(int writeend);
    
    void checkstatus();
    
    int issueflag(int flag);
    
    const char *inputfile, *outputfile, *reffile;
    
    long long line_index;
    
    std::string baifilename, barnamefile, barcodefile, samplefile;
    
    Snp_Records snps_records;
    
    Barcode_Records Barcodes_Records;
    
    ref_Reader reader;
    
    MDTag MDtag;
    
};


#endif /* Fileloader_hpp */
