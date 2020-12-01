//
//  stats.cpp
//  snp_caller
//
//  Created by Wangfei MA on 6/30/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "sampling.hpp"

#define  NucleoIndex( base , index) for ( index=0;index<8;++index) {if ("ATCGDNMI"[index]==base) break;}

#define  NucleoBase(theindex) "ATCGDNMI"[theindex]

using namespace std;

extern float coverage;

extern float pre_errors[8];

extern long long total_alleles[8];

extern int sign;

extern bool ifMD;

extern float pre_qualcuts[8];

extern int turn;


int find_compli_base(int base)
{
    
    if (base>3) return base;
    
    return 2*(base/2)+(1-base%2);
    
}

void find_average_allele(int alleles[], int correct_base)
{
    
    alleles[6]+=(alleles[0]+alleles[1]+alleles[2]+alleles[3]+2);//all mutated single indels;
    if (correct_base<4 || correct_base==6)
    {
        
        alleles[6]-=alleles[correct_base];
        alleles[6]/=3;
    }
    else
    {
        alleles[6]/=4;
    }
}


int highest_intereface(int alleles[], int first_base, int mutated_base)
{
    
    int intereface=0;
    if (mutated_base==6 || first_base==6)
    {
        
        if (first_base==4 || mutated_base==4) return alleles[7];
        else if (first_base==7 || mutated_base==7) return alleles[4];
        else return (alleles[4]>alleles[7])?alleles[4]:alleles[7];
    }
    
    else
    {
        
        for (int i=0;i<8;i++)
        {
            
            if (i==mutated_base || i==first_base || i==6 || i==5) continue;
            if (intereface<alleles[i])
            {
                intereface=alleles[i];
            }
        }
    }
    
    return intereface;
}


sampling::sampling()
{
    
    valide_start=3*coverage/4;
    
    valide_end=2*coverage-valide_start+1;
    
    valide_size=valide_end-valide_start;
    
    for (int j=0;j< 16*64;j++)
    {
        error_counts.push_back((long long*) calloc(2*coverage, 2*sizeof(long long)));
    }
    
    
    for (int j=0;j<9;j++)
    {
        quals_errors.push_back((long long*) calloc(2000, 2*sizeof(long long)));
        
        quals_counts.push_back((long long*) calloc(2000, 2*sizeof(long long)));
        
        allele_counts.push_back((long long*) calloc(valide_end+1, 2*sizeof(long long)));
        
    }
}

int sampling::get_consecutive(int current_index, char allele)
{
        
    if (current_index<10 || current_index>refgenome.length()-10) return 9;
    
    int motif_index=0, left_base=0, right_base=0,lconsecutive=0,rconsecutive=0;
    
    for (int coordi=2;coordi<5;coordi++)
    {
        
        if (toupper(refgenome[current_index-coordi]) != allele)
        {
            NucleoIndex(toupper(refgenome[current_index-coordi]), left_base)
            break;
        }
        
        else
        {
            lconsecutive++;
        }
    }
    
    for (int coordi=-1;coordi<2;coordi++)
    {
        if (toupper(refgenome[current_index+coordi+1]) != allele)
        {
            NucleoIndex(toupper(refgenome[current_index-coordi]), right_base)
            break;
        }
        else
        {
            rconsecutive++;
        }
    }
    
    if (lconsecutive>2 || rconsecutive>2) return 9;
    
    motif_index=lconsecutive*3+rconsecutive;
    
    return motif_index;
    
}

void sampling::loadref(const char *reffile, const char *chr)
{
    
    ref_Reader Reader(reffile);
    
    refgenome=Reader.Load_Chr(reffile,chr);
    
}


int sampling::Read_alleles(string line, int line_data[40], int &allele_sum, int &allele_max, bool ifmask)
{
    
    int line_size=(int)line.length();
    
    if (!line_size) return -1;
    
    int max_index=0, new_allele=0, new_index=0;
    
    int char_index=0;
    char current_char=',';
    
    int alleles[8]={0};
    
    for (;char_index<line_size;++char_index)
    {
        
        current_char=line[char_index];
        
        if (current_char==',')
        {
            
            line_data[new_index]=new_allele;
            
            if (new_index<16) alleles[new_index/2]=new_allele;

            new_index++;
            new_allele=0;
        }
        
        else if (current_char>='0' && current_char<='9')
        {
            
            new_allele*=10;
            new_allele+=current_char-'0';
        }
    }
    
    max_index=(int)(max_element(alleles, alleles + 8) -alleles);
    
    allele_max=line_data[2*max_index]+line_data[2*max_index+1];
    
    allele_sum=(int) accumulate(alleles, alleles+8, 0);
    
    return max_index;
    
}

size_t sampling::readsamples(string samplefile)
{
    
    int alleles[40]={0};
    
    totalnum=0;
    
    
    cout<<"start loading"<<endl;
    
    string line;
    int correct_base, compli_base, error_count_f, error_count_r, sum_count, max_count, rest_count, mean_qual_f, mean_qual_r, lastcoordi=1, nextcoordi, lastbaseindex=-1, consecutive_index_f,consecutive_index_r, rest_cutoff=0;
    
    char ref_base, nextbase;
    
    fstream file;
    file.open(samplefile, ios_base::in);
    if (!file) {
        std::cerr << "ERROR: Could not open " << samplefile << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    
    while (getline(file,line,'\n'))
    {
        
        if (!(++totalnum%1000000))
        {
            printf("Processing line: %lld\n",totalnum);
        }
        
        max_count=0;sum_count=0;
        correct_base=Read_alleles(line, alleles, sum_count, max_count, 0);
        
        ref_base=toupper(refgenome[alleles[0]-1]);
        
        nextcoordi=alleles[33];
        nextbase= NucleoBase(alleles[34]);
        
        consecutive_index_f=4*lastbaseindex+alleles[34];
        
        consecutive_index_r=find_compli_base(lastbaseindex)+4*find_compli_base(alleles[34]);
        
        if (correct_base<0 ||correct_base==5 ||!max_count)
        {
            continue;
        }
        
        compli_base=find_compli_base(correct_base);
        
        rest_cutoff=min(max_count/4,(int)coverage/2);
        
        rest_count=sum_count-max_count;
        
        if (alleles[2*correct_base+1])
        {
            mean_qual_f=(10*alleles[2*correct_base+17]/((alleles[2*correct_base+1]+1)/2)-320);
        }
        else
        {
            mean_qual_f=0;
        }
        
        if (alleles[2*correct_base+2])
        {
            mean_qual_r=(10*alleles[2*correct_base+18]/((alleles[2*correct_base+2]+1)/2)-320);
        }
        else
        {
            mean_qual_r=0;
        }
        
        if (max_count>0.8*sum_count && max_count<=valide_end && rest_count<rest_cutoff)
        {
            
            allele_counts[correct_base][alleles[2*correct_base+1]]++;
            
            allele_counts[compli_base][alleles[2*correct_base+2]]++;
            
            if (min(mean_qual_f,mean_qual_r)>=0 && max(mean_qual_f,mean_qual_r)<2000)
            {
                quals_counts[correct_base][mean_qual_f]++;
                quals_counts[compli_base][mean_qual_r]++;
                
                if (correct_base<4)
                {
                    quals_counts[5][mean_qual_f]++;
                    quals_counts[5][mean_qual_r]++;
                }
                
            }
            
        }
        
        if (lastcoordi!=alleles[0]-1 || nextcoordi!=alleles[0]+1 || toupper(refgenome[nextcoordi-1])!= nextbase || toupper(refgenome[lastcoordi-1])!= NucleoBase(lastbaseindex) || ref_base!=NucleoBase(correct_base) || max_count<sum_count/2 || max_count>=2*coverage  /*mean_qual_f<pre_qualcuts[correct_base] || mean_qual_r<pre_qualcuts[compli_base] */)
            
        {
            lastcoordi=alleles[0];
            lastbaseindex=correct_base;
            continue;
        }
        
        lastcoordi=alleles[0];
        lastbaseindex=correct_base;
        
        //find_average_allele(alleles, correct_base);
        
        for (int mutated_base=0;mutated_base<5;mutated_base++){
            
            if (mutated_base==5 || (mutated_base==correct_base) || (mutated_base==4&&line[0]=='-')){continue;}
            
            error_count_f=alleles[2*mutated_base+1];
            error_count_r=alleles[2*mutated_base+2];
            
            if (error_count_f || error_count_r)
            {
                
                if (error_count_f)
                {
                    mean_qual_f=(10*alleles[mutated_base+17]/((error_count_f+1)/2)-320);
                }
                else
                {
                    mean_qual_f=0;
                }
                
                if (error_count_r)
                {
                    mean_qual_r=(10*alleles[mutated_base+18]/((error_count_r+1)/2)-320);
                }
                else
                {
                    mean_qual_r=0;
                }

                
                if (min(mean_qual_f,mean_qual_r)>=0 && max(mean_qual_f,mean_qual_r)<2000)
                {
                    quals_errors[correct_base][mean_qual_f]++;
                    quals_errors[compli_base][mean_qual_r]++;
                    
                    if (correct_base<4)
                    {
                        quals_errors[5][mean_qual_f]++;
                        quals_errors[5][mean_qual_r]++;
                    }
                    
                }
                
            }
            
            //if (highest_intereface(alleles, correct_base, mutated_base)>rest_cutoff){continue;}
            
            error_counts[consecutive_index_f*64+8*correct_base+mutated_base][error_count_f]++;
            
            error_counts[consecutive_index_f*64+8*6+mutated_base][error_count_f]++;
            
            error_counts [consecutive_index_r*64+8*compli_base+find_compli_base(mutated_base)][error_count_r]++;
            
            error_counts[consecutive_index_r*64+8*6+find_compli_base(mutated_base)][error_count_r]++;
            
        }
    }
    
    file.close();
    
    for (int i=valide_start; i<valide_end; i++)
    {
        
        //get all alleles coverage
        for (int j=0;j<8;j++)
        {
            //sum all alleles to N
            allele_counts[5][i-valide_start]+=allele_counts[j][i-valide_start];
        }
        
        //sum all ATCG coverage to M
        for (int j=0;j<4;j++)
        {
            allele_counts[6][i-valide_start]+=allele_counts[j][i-valide_start];
        }
        
    }
    
    for (int i=0; i<8; i++)
    {
        
        for (int j=0; j<valide_size; j++)
        {
            
            total_correct_alleles[i]+=allele_counts[i][j];
            total_correct_moles[i]+=allele_counts[i][j]*(j+valide_start);
            
        }
    }
    
    return totalnum;
    
}




bool sampling::statistics(){
    
    
    char alleles[8]={'A','T','C','G','D','N','M','I'};
    
    int consecutive, correct_allele, mutated_allele;
    
    cout<<"start modeling"<<endl;
    
    for (int i=0; i<8; i++)
    {
        
        int qualcut, qualcut_low;
        
        //Aprint(quals_counts[i], long long, 2000)
        
        //Aprint(quals_errors[i], long long, 2000)
        
        qual_anal.modeling(vector <int>(quals_counts[i], quals_counts[i]+2000), qualcut,qualcut_low);
        
        qual_cuts[alleles[i]]=qualcut;
        qual_cuts_low[alleles[i]]=qualcut_low;
    }
    
    for (int i=0; i<8; i++)
    {
        cover_curve allele_cover;
        
        if (!cover_anal.modeling(vector <double> (allele_counts[i],allele_counts[i]+valide_size), valide_start, allele_cover) && total_correct_alleles[i]>coverage/4)
        {
            
            allele_cover.type="estimate";
            allele_cover.mean=total_correct_moles[i]/total_correct_alleles[i];
            allele_cover.std_log=log(allele_cover.mean)/2;
            
        }
        
        if (fabs(log(allele_cover.mean)-log(coverage))>0.5)
        {
            
            allele_cover.type="default";
            allele_cover.mean=coverage;
            allele_cover.std_log=log(allele_cover.mean)/2;
            
        }
        
        allele_coverages[alleles[i]]=allele_cover;
        
    }
    
    for (int i=0; i<16*64; i++)
    {
        error  type_error;
        
        consecutive=i/64;
        correct_allele=(i%64)/8;
        mutated_allele=i%8;
        
        //if (correct_allele==5 || mutated_allele==5 || (correct_allele!=6 && correct_allele==mutated_allele))
        if (correct_allele>3)
        {
            continue;
        }
        
        vector <double>  frequencies(error_counts[i], error_counts[i]+2*(int)coverage);
        
        string error_type= string(1,alleles[consecutive/4])+string(1,alleles[correct_allele])+string(1,alleles[consecutive%4]) +'-'+string(1,alleles[mutated_allele]);
        
        //if (error_type.compare("CGG-A")) continue;
                
        type_error.errorrate=error_anal.estimate_error(frequencies, (int) coverage);
        
        error_anal.modeling(frequencies,(int)coverage, type_error);
        
        type_errors[error_type]=type_error;
        
    }
    cout<<"end modeling"<<endl;
    
    return 1;
}


void sampling::output(string outputfile){
    
    const char *sign_chr="";
    if (sign<0)
    {
        sign_chr="-";
    }
    
    char alleles[8]={'A','T','C','G','D','N','M','I'};
    
    ofstream file;
    
    file.open (outputfile,std::ofstream::app);
    
    file<<endl;
    
    for (int i=0; i<8; i++)
    {
        cover_curve *current_model=&allele_coverages[alleles[i]];
        
        if (!current_model->type.compare("Poisson"))
        {
            
            file<<sign_chr <<alleles[i]<<"\tqualcut:"<<0.1*qual_cuts[alleles[i]]<<"\tlowcut:"<<0.1*qual_cuts_low[alleles[i]]<< "\tcoverage:"<<current_model->mean<<endl;
            
        }
        else
        {
            file<<sign_chr<<alleles[i]<<"\tqualcut:"<<0.1*qual_cuts[alleles[i]]<<"\tlowcut:"<<  0.1*qual_cuts_low[alleles[i]] <<"\tcoverage:"<<current_model->mean<< "\tstdev:"<<current_model->std_log<<endl;
        }
    }
    
    int correct_allele, mutated_allele,consecutive;
    for (int i=0; i<16*64; i++)
    {
        
        consecutive=i/64;
        
        correct_allele=(i%64)/8;
        mutated_allele=(i%64)%8;
        
        if (correct_allele==5 || mutated_allele==5 || (correct_allele!=6 && correct_allele==mutated_allele))
        {
            
            continue;
        }
        
        
        string error_type{alleles[consecutive/4] ,alleles[correct_allele],alleles[consecutive%4],'-',alleles[mutated_allele]};
        
        error* current_model=&type_errors[error_type];
        
        file <<sign_chr<<error_type<<"\terror_cutoff:"<< current_model->errorcut<< "\tlowconfi_cutoff:" << current_model->lowconfi_cut<< "\terror_rate:"<<current_model->errorrate<<"\terror_curve:"<<current_model->curve_quadratic[0]<<","<<current_model->curve_quadratic[1]<<endl;
        
        
    }
    
    file.close();

}
