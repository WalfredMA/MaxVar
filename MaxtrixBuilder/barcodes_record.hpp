//
//  barcodes_record.hpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef barcodes_record_hpp
#define barcodes_record_hpp

#include <iostream>
#include <vector>
#include <unordered_map>

#include "Config.h"


class Barcode_Records
{
    
public:
    
    Barcode_Records(const char* file);
    ~Barcode_Records(){};
    
    unsigned int push_back(long long barcode);
    
    unsigned int push_back(const char* barname);
 
    void Update_Infor(int molesize, int start, int end);
    
    void Write_Barcodes(long long current_posi);
    
	void Write_Barnames();
    
    void Outputall();

    unsigned int barcode_counter, singlemole_counter;
    
    char *current_barname;
    
    unsigned int current_barindex;
    
    long long current_barcode;
    
    std::string barnamefile;
    
    std::string barcodefile;
    
private:
    
    void Addbarcode();
    
    void Addbarname();
    
    void Check_malloc();
    
    unsigned int Record_Search(long long barcode);
    
    bool Clear_Barcode(long long barcode);
   
    char* Allbarnames[1000];
    long long Allbarnames_barcodes[1000];
    
    unsigned int barname_counter, current_hash0, current_hash1;
    
    std::unordered_map<long long, unsigned int> barcode_hashcells;
    
    long long hash(const char*, int);
    
    unsigned int deleted_barcode;
 
    std::vector<Molecule*>Allbarcodes;
    std::vector<bool>ifvalide_barcodes;
};




#endif /* barcodes_record_hpp */
