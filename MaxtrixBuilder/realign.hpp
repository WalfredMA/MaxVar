//
//  realign.hpp
//  testcpp
//
//  Created by Wangfei MA on 8/10/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef realign_hpp
#define realign_hpp

#include <vector>
#include <set>
#include <map>

#include "Config.h"

class Realigner{

public:
    
    Realigner(){};
    ~Realigner(){};
    
    void Realign(int newtarget);
    
    void Buffer_Region(int start, int end);
        
private:
    
    void Realign_Allele();
    
    bool Find_SameSNP(std::vector<int>&,std::vector<std::vector<int>>&);
    
    bool Generate_newRef();
    
    void Realign(std::vector<int>& othersnps, std::vector<std::vector<int>>& othersnps_moles);
    
    int Distract_Mole(int pos, int mole, int allele, Seed*, int);
    
    bool Summary_Moles(std::vector<std::vector<int>>& ref_mole,std::vector<int>& ref);
    
    bool Realign_newRef(Seed, int);
    
    Base** buffered_region;
    
    std::atomic_uchar* buffered_region_ifinit;
    
    void Add_History(int target_pos, int othersnp, int iallele, std::vector<int>& moles);
    int Search_History(int findpos, int imole,std::vector<int>& des, int sign);
    int Search_Buffer(int index, int imole, std::vector<int>& des, int sign);
        
    Base *current_snp;
    
    int buffer_start, search_end, write_end, target_pos, target_allele, minor_allele, major_allele, num_ref_Molecules, deletion_start, deletion_end;
    
    std::set<int> ref_Molecules;
        
    Seed ref;
};

#endif /* realign_hpp */


