//
//  matrix_reader.hpp
//  snp_caller
//
//  Created by Wangfei MA on 8/28/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef matrix_reader_hpp
#define matrix_reader_hpp

#include <iostream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <string>
#include <set>
#include <algorithm>
#include <numeric>
#include <cstring>


typedef struct
{
    
    std::vector<int> f_moles;
    std::vector<int> r_moles;
    
    std::vector<char> f_map_scores;
    std::vector<char> r_map_scores;
    
    std::vector<int> f_base_scores;
    std::vector<int> r_base_scores;
    
    float f_base_score;
    float r_base_score;
    
    int f_unmask;
    int r_unmask;
    
}Allele_infor;

typedef struct
{
    int coordinate;
    int allele_counts[10];
    int alleles_blocked[10];
    Allele_infor Allele[8];
    int max_allele,max_count;//index(base type), molecule count, unmasked molecule count, mean base quality
    int second_allele,second_count;
    int third_allele,third_count;
    int total,total_f,total_r;
    int total_raw;
    std::vector<int> sv_del, sv_ins;
    int sv_del_f, sv_del_r,sv_ins_f, sv_ins_r;
    std::vector<char> del_map, ins_map;
    std::set<int> alleles_found;
    std::set<int> alleles_lowconfi;
    std::vector<char> alleles_label;
    
}Base_infor;


class matrix_reader
{
    
public:
    
    matrix_reader(){};
    
    void ReadLine(std::string& StrLine, Base_infor& snp_infor);
    
};

#endif /* matrix_reader_hpp */
