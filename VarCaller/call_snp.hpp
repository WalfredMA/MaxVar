//
//  call_snp.hpp
//  snp_caller
//
//  Created by Wangfei MA on 7/1/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef call_snp_hpp
#define call_snp_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <limits>


#include "probability_calculator.hpp"
#include "matrix_reader.hpp"
#include "ref_reader.hpp"

using namespace std;


class sv_buffer
{
    
public:
    
    sv_buffer(int coordi, std::vector<int> set1, std::vector<int> set2, std::string seq):line_coordi(coordi),allele_1(set1),allele_2(set2),sequences(seq) {};
    ~sv_buffer()
    {
        free(ifoutput);
        free(last_coordi);
    }
    
    std::string sequences;
    std::vector<int> allele_1;
    std::vector<int> allele_2;
    int line_coordi;
    bool *ifoutput;
    int *last_coordi;
    
};




class snp_caller
{
    
public:
    
    snp_caller(prob_calculator &cal);
    
    int snp_calling(const char* inputfile, const char* outputfile, const char* reffile, const char* chr);
    
    int phase_cutoff=40, baseconfi_cutoff=15, postphase_cutoff=35, indel_cutoff=40;
    
private:
    
    int Var_Call(std::string &inLine);

    int filtration(Base_infor& current_base);
    
    int Allele_Call(Base_infor& current_base, int cutoff);
    
    int SV_Call(Base_infor& current_base);
    
    int Indel_Call(Base_infor& current_base, int cutoff);
    
    float Base_evaluate(Base_infor& current_base);
    
    int Allele_evaluate(Allele_infor& current_allele, int coordinate, int allele_index, int f_consecutive, int r_consecutive, float &score);
    
    float Strd_evaluate(Allele_infor& current_allele);
    
    float Haplo_evaluate(int allele_index);
    
    float Prop_evaluate(int allele_index);
    
    float Mask_evaluate(Allele_infor& current_allele, int current_coordi, int allele_index);
    
    float Qual_evaluate(Allele_infor& current_allele, int allele_index, int f_cutoff, int r_cutoff, int f_hcutoff, int r_hcutoff); 
    
    float Freq_evaluate(Allele_infor& current_allele, int allele_index, int f_consecutive, int r_consecutive, int &f_cutoff, int &r_cutoff, int &f_hcutoff, int &r_hcutoff);
    
    int Base_phasegroup(Base_infor& current_base);
    
    float get_consecutive(int current_coordi, char allele);
    
    int motif_index(int current_coordi, int& f_consecutive_index, int& r_consecutive_index);
    
    void output_buffer(FILE *fwrite);
        
    void find_top_alleles(int *top_values, int *top_indexs, int new_value, int new_index);
    
    bool *add_buffer(sv_buffer& new_buffer);
    
    int MapFilter(Allele_infor& allele);
    
    prob_calculator calculator;
    
    matrix_reader reader;
    
    Base_infor current_base;
    
    std::vector<sv_buffer*> sv_buffers;
    
    sv_buffer* new_buffer;
    
    std::string refgenome;
    
    bool ifconnect=1;
    
    int valide_start, valide_end, valide_size;
    
    int max_label;
    float posi_errors[10], nega_errors[10], allele_scores[10];
    int allele_phases[10];
    float base_score;
    
    vector<float> output_scores;
    
};


#endif /* call_snp_hpp */
