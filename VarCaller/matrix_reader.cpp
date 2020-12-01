//
//  matrix_reader.cpp
//  snp_caller
//
//  Created by Wangfei MA on 8/28/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "matrix_reader.hpp"

#define  Vprint( n )  for ( auto i: n ) { cout<<(int)i<<" ";} cout<<endl;

using namespace std;

extern float mean_coverage;

float higher_mean(vector<int>& nums)
{
    
    if (!nums.size()) return 0;
    
    sort(nums.begin(), nums.end());
    
    return 1.0*(accumulate(nums.begin()+nums.size()/2, nums.end(), 0))/((int)(nums.size()+1)/2);
}


class top_alleles
{
    
public:
    top_alleles(int* highest, int* highest_value, int* second, int* second_value, int* third, int* third_value): highest_index(highest), second_index(second), third_index(third),
    highest_count(highest_value), second_count(second_value), third_count(third_value)
    {
        
        *highest_count=-1;
        *second_count=-1;
        *third_count=-1;
    };
    
    void find_top_alleles(int new_value, int new_index){
        
        if (new_value>*highest_count)
        {
            *third_count=*second_count;
            *third_index=*highest_index;
            
            *second_count=*highest_count;
            *second_index=*highest_index;
            
            *highest_count=new_value;
            *highest_index=new_index;
        }
        
        else if (new_value>*second_count)
        {
            *third_count=*second_count;
            *third_index=*highest_index;
            
            *second_count=new_value;
            *second_index=new_index;
            
        }
        
        else if (new_value>*third_count)
        {
            
            *third_count=new_value;
            *third_index=new_index;
        }
        
    }
    
    int *highest_index, *second_index,*third_index;
    
    int *highest_count, *second_count,*third_count;
};



void grep_molecules(string column, int &strindex ,int& total_allele, Allele_infor &allele_infor, std::vector<int> &sv_del, std::vector<char> &del_map, int &sv_del_f, int &sv_del_r , std::vector<int> &sv_ins, std::vector<char> &ins_map, int &sv_ins_f, int &sv_ins_r)
{
    
    std::vector<int> &f_moles = allele_infor.f_moles;
    std::vector<int> &r_moles = allele_infor.r_moles;
    f_moles.clear();r_moles.clear();
    
    std::vector<char> &f_map_scores = allele_infor.f_map_scores;
    std::vector<char> &r_map_scores = allele_infor.r_map_scores;
    f_map_scores.clear();r_map_scores.clear();
    
    std::vector<int> &f_base_scores = allele_infor.f_base_scores;
    std::vector<int> &r_base_scores = allele_infor.r_base_scores;
    f_base_scores.clear();r_base_scores.clear();
    
    float &f_base_score = allele_infor.f_base_score;
    float &r_base_score = allele_infor.r_base_score;
    
    int &f_unmask = allele_infor.f_unmask;
    int &r_unmask = allele_infor.r_unmask;
    
    
    f_unmask=0;r_unmask=0;

    char current_chr=NULL, last_chr=',', map_score='A';
    int curren_num=0,current_sign=1,current_qual=0;
    
    int sv_size=0,one_side_size=0;
    
    srand(time(0));
    
    for (;strindex< column.length();++strindex){
        
        current_chr=column[strindex];
        
        if (current_chr=='\t')
        {
            ++strindex;
            break;
        }
        
        switch (last_chr)
        {
            case ',':
                current_qual=(int)(current_chr-32);
                last_chr='\0';
                continue;
                break;
            
            case '\0':
                map_score=current_chr;
                last_chr='a';
                continue;
                break;
                
            case ':':
                switch (current_chr)
                {
                    case ',':
                        break;
                    case '-':
                        sv_size=one_side_size;
                        one_side_size=0;
                        continue;
                        break;
                    case '0'...'9':
                        one_side_size*=10;
                        one_side_size-=(int)(current_chr-'0');
                        continue;
                        break;
                    default:
                        sv_size+=2;
                        continue;
                        break;
                }
                break;
                
            default:
                break;
        }
        
        switch (current_chr)
        {
            
            case '0'...'9':
                curren_num*=10;
                curren_num+=(int)(current_chr-'0');
                break;
            
            case ':':
                last_chr=':';
                one_side_size=0;
                sv_size=0;
                break;
                
            case '-':
                current_sign=-1;
                break;
            
            case ',':
                
                sv_size+=one_side_size;
                last_chr=',';
                
                if (rand()%10>5){
                
                if (map_score<'a')
                {
                    
                    if (sv_size>10)
                    {
                        sv_ins.push_back(current_sign*curren_num);
                        ins_map.push_back(map_score);
                        sv_ins_f++;
                    }
                    else if (sv_size<-10)
                    {
                        sv_del.push_back(current_sign*curren_num);
                        del_map.push_back(map_score);
                        sv_del_f++;
                    }
                    
                    else
                    {
                        f_moles.push_back(current_sign*curren_num);
                        
                        f_map_scores.push_back(map_score);
                        
                        f_base_scores.push_back(current_qual);
                        
                        if (current_sign>0) f_unmask++;
                    }
                    
                }
                else
                {
                    
                    if (sv_size>10)
                    {
                        sv_ins.push_back(current_sign*curren_num);
                        ins_map.push_back(map_score);
                        sv_ins_r++;
                    }
                    else if (sv_size<-10)
                    {
                        sv_del.push_back(current_sign*curren_num);
                        del_map.push_back(map_score);
                        sv_del_r++;
                    }
                    
                    else
                    {
                        r_moles.push_back(current_sign*curren_num);
                        
                        r_map_scores.push_back(map_score);
                        
                        r_base_scores.push_back(current_qual);
                        
                        if (current_sign>0) r_unmask++;
                    }
                }
                }
                current_sign=1;
                curren_num=0;
                current_qual=0;
                map_score='A';
                break;
                
            default:
                continue;
                break;
        }
        
    }
    
    total_allele=(int)(f_moles.size()+r_moles.size());
    
    f_base_score=higher_mean(f_base_scores);
    
    r_base_score=higher_mean(r_base_scores);
    
        
}


void matrix_reader::ReadLine(string& StrLine, Base_infor& snp_infor)
{
    
    top_alleles tops(&snp_infor.max_allele, &snp_infor.max_count, &snp_infor.second_allele, &snp_infor.second_count, &snp_infor.third_allele, &snp_infor.third_count);
    
    snp_infor.sv_del.clear();
    snp_infor.sv_ins.clear();
    snp_infor.del_map.clear();
    snp_infor.ins_map.clear();
    snp_infor.alleles_found.clear();
    snp_infor.alleles_lowconfi.clear();
    snp_infor.alleles_label.clear();
    snp_infor.total=0;snp_infor.total_raw=0;snp_infor.total_f=0;snp_infor.total_r=0;
    snp_infor.sv_del_f=0; snp_infor.sv_del_r=0;
    snp_infor.sv_ins_f=0; snp_infor.sv_ins_r=0;
    snp_infor.max_allele=0;
    snp_infor.max_count=0;
    snp_infor.second_allele=0;
    snp_infor.second_count=0;
    snp_infor.third_allele=0;
    snp_infor.third_count=0;
    snp_infor.coordinate=0;
    memset(snp_infor.allele_counts,0,10*sizeof(int));
    memset(snp_infor.alleles_blocked,0,10*sizeof(int));
    
    
    vector<int> allele_molecules[8];
    vector<char> allele_maps[8];
    
    int seekpos=0; char current_chr;
    while(seekpos<StrLine.length())
    {
        current_chr=StrLine[seekpos++];
        
        if (current_chr=='\t') break;
        
        if (current_chr>='0' && current_chr<='9')
        {
            snp_infor.coordinate*=10;
            snp_infor.coordinate+=current_chr-'0';
        }
    }
    
    if (StrLine[0]=='-') snp_infor.coordinate*=-1;
    
    for (int column=0;column<8;column++)
    {
        
        grep_molecules(StrLine, seekpos, snp_infor.allele_counts[column], snp_infor.Allele[column], snp_infor.sv_del, snp_infor.del_map, snp_infor.sv_del_f, snp_infor.sv_del_r, snp_infor.sv_ins,  snp_infor.ins_map, snp_infor.sv_ins_f, snp_infor.sv_ins_r);
        
        if ( column<4 )
        {
           tops.find_top_alleles(snp_infor.allele_counts[column], column);
        }
                
        snp_infor.total+=snp_infor.Allele[column].f_unmask+snp_infor.Allele[column].r_unmask;
        snp_infor.total_raw+=snp_infor.allele_counts[column];
        
        snp_infor.total_f+=snp_infor.Allele[column].f_moles.size();
        snp_infor.total_r+=snp_infor.Allele[column].r_moles.size();
    }
    
    snp_infor.total_f+=snp_infor.sv_del_f+snp_infor.sv_ins_f;
    snp_infor.total_r+=snp_infor.sv_del_r+snp_infor.sv_ins_r;
    snp_infor.allele_counts[8]=(int)snp_infor.sv_del.size();
    snp_infor.allele_counts[9]=(int)snp_infor.sv_ins.size();

} 

