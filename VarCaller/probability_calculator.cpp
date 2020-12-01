//
//  probability_calculator.cpp
//  snp_caller
//
//  Created by Wangfei MA on 7/8/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "probability_calculator.hpp"

#define  NucleoIndex( base , index) for ( index=0;index<8;++index) {if ("ATCGDNMI"[index]==base) break;}

#define  Vprint( n )  for ( auto i: n ) { cout<<i<<" ";} cout<<endl;


extern float median_coverage;
float median_coverage=0.0;

extern float mean_coverage;
float mean_coverage=0.0;

extern int base_qualcuts[8];

extern int base_hqualcuts[8];

extern float confi_cutoff;
float confi_cutoff=7;

extern int confidence;

extern int icycle;

using namespace std;

float log_add(float x, float y)
{
    
    
    float higher=(x>y)?x:y;
    
    float diff=fabs(x-y);
    
    if (diff>3){
        
        return higher;
        
    }
    
    return higher+log(1+exp(-diff));
    
}

void prob_calculator::loadStat(string statfile){
    
    for (int i=0;i<8;++i)
    {
        
        alleles[i]=new allele;
    }
    
    for (int i=0;i<1024;++i)
    {
        mutations[i]=new mutation;
    }
    
    allele* current_allele;
    mutation* current_error;
    
    fstream file;
    file.open(statfile, ios_base::in);
    
    if (!file) cout<<"Could not find file "<<statfile<<endl;
    
    vector<long long> counters;
    
    string line;
    
    getline(file,line);
    
    
    int i=0;
    while(getline(file,line)&&i<10)
    {
        
        if (line.size()) counters.push_back((long long) atoll(line.c_str()));
        i++;
    }
    
    long long allelecounts=accumulate(counters.begin()+2, counters.end(), 0L);
    
    mean_coverage=(float)(allelecounts)/counters[0];
    
    int default_cutoff=min((int)mean_coverage/8,(int)floor(21/confidence));
    
    int default_hqcutoff=min((int)mean_coverage/4,(int)floor(31/confidence));
    
    size_t pos1, allele_index, mutation_index, last_index, next_index;
    
    int length;
    while (getline(file,line))
    {
        
        length=(int)line.length();
        if (length<2)continue;
        
        if (line[3]!='-')
        {
            
            NucleoIndex(line[0],allele_index);
            
            current_allele=alleles[allele_index];
            
            int &current_qualcut=base_qualcuts[allele_index];
            int &current_hqualcut=base_hqualcuts[allele_index];
            
            if ((pos1 = (int)line.find("coverage:",2)) != std::string::npos )
            {
                
                current_allele->coverage_exp=atof(line.c_str()+pos1+9);
                
                
            }
            
            if ((pos1 = (int)line.find("qualcut:",2)) != std::string::npos)
            {
                
                current_allele->qualcut= atof(line.c_str()+pos1+8);
                current_hqualcut=current_allele->qualcut;
                
            }
            
            if ((pos1 = (int)line.find("lowcut:",2)) != std::string::npos)
            {
                
                current_allele->lowqual= atof(line.c_str()+pos1+7);
                current_qualcut=current_allele->lowqual;
                
            }
            
            if ((pos1 = (int)line.find("stdev:",0)) != std::string::npos)
            {
                
                current_allele->coverage_errorlog= atof(line.c_str()+pos1+6);
                
            }
            
            else
            {
                current_allele->coverage_errorlog=log(current_allele->coverage_exp)/2;
                current_allele->ifpoisson=0;
            }
            
        }
        
        else
        {
            
            NucleoIndex(line[0],last_index);
            NucleoIndex(line[1],allele_index);
            NucleoIndex(line[2],next_index);
            NucleoIndex(line[4],mutation_index);
            
            current_error= mutations[256*last_index+64*next_index+8*allele_index+mutation_index];
            
            const char *grep_string;
            
            grep_string="lowconfi_cutoff:";
            if ((pos1 = (int)line.find(grep_string,2)) != std::string::npos)
            {
                
                current_error->errorcut= atoi(line.c_str()+pos1+(int)strlen(grep_string));
                
                if (current_error->errorcut<=0 || current_error->errorcut<default_cutoff)
                {
                    
                    current_error->errorcut=default_cutoff;
                    
                }
                
            }
            
            grep_string="error_cutoff:";
            if ((pos1 = (int)line.find(grep_string,2)) != std::string::npos)
            {
                
                current_error->highcut= atoi(line.c_str()+pos1+(int)strlen(grep_string));
                
                if (current_error->highcut<=0 || current_error->highcut<default_hqcutoff)
                {
                    
                    current_error->highcut=default_hqcutoff;
                    
                }
                
            }
        }
    }
    
    file.close();
    
    for (allele *current_allele: alleles)
    {
        if (current_allele->qualcut<alleles[5]->qualcut) current_allele->qualcut=alleles[5]->qualcut;
        if (current_allele->lowqual<alleles[5]->lowqual) current_allele->lowqual=alleles[5]->lowqual;
    }
    
    
}


void prob_calculator::loadSNP(string snpfile)
{
    
    map<int, vector<int>> snp_sites;
    
    fstream file1;
    file1.open(snpfile, ios_base::in);
    
    if (!file1)
    {
        cout<<"Could not find file "<<snpfile<<endl;
    }
    
    string line;
    
    while(getline(file1,line))
    {
        
        char chr='-';
        int coordinate=0,i=0;
        for (i=0;i<line.size();++i)
        {
            chr=line[i];
            
            if (chr>='0' &&chr<='9')
            {
                coordinate*=10;
                coordinate+=(chr-'0');
            }
            else
            {
                i++;
                break;
            }
        }
        
        int current_num=0;
        for (;i<line.size();++i)
        {
            chr=line[i];
            
            if (chr>='0' &&chr<='9')
            {
                current_num*=10;
                current_num+=(chr-'0');
            }
            else if (chr==',')
            {
                snp_sites[current_num].push_back(coordinate);
                current_num=0;
            }
            else if (chr=='\t'||chr=='\n')
            {
                i++;
                break;
            }
        }
        
        if (current_num) snp_sites[current_num].push_back(coordinate);
        
        for (;i<line.size();++i)
        {
            
            chr=line[i];
            
            if (chr=='\t'|| chr=='\n') break;
            
            if ((chr>='A' && chr<='Z') || (chr>='a' && chr<='z'))
            {
                labels[coordinate].push_back(chr);
            }
                
        }
        
    }
    
    for (auto &group: snp_sites)
    {
        
        if (group.second.size()>3 && (*max_element(group.second.begin(), group.second.end())-*min_element(group.second.begin(), group.second.end())-1)>100)
        {
            used_group.insert(group.first);
            std::copy(group.second.begin(), group.second.end(), std::inserter(refsnp, refsnp.end()));
        }
        else
        {
            std::copy(group.second.begin(), group.second.end(), std::inserter(lowsnp, lowsnp.end()));
        }
        
    }
    
}


void prob_calculator::loadMole(string molefile)
{
    
    fstream file2;
    file2.open(molefile, ios_base::in);
    
    if (!file2)
    {
        cout<<"Could not find file "<<molefile<<endl;
    }
    
    string line;
    while(getline(file2,line))
    {
        
        char chr;
        int molecule=0,i=0;
        for (i=0;i<line.size();i++)
        {
            chr=line[i];
            
            if (chr>='0' &&chr<='9')
            {
                molecule*=10;
                molecule+=(chr-'0');
            }
            else
            {
                i++;
                break;
            }
        }
        
        int current_num=0,sign=1;
        for (;i<line.size();i++)
        {
            chr=line[i];
            
            if (chr>='0' &&chr<='9')
            {
                current_num*=10;
                current_num+=(chr-'0');
            }
            else if (chr=='-')
            {
                sign=-1;
            }
            else if (chr==',' || chr=='\t')
            {
                
                if (current_num && used_group.count(current_num))
                {
                    mole_phase[molecule].push_back(current_num*sign);
                }
                
                current_num=0;
                sign=1;
            }
        }
        if (current_num && used_group.count(current_num))
        {
            mole_phase[molecule].push_back(current_num*sign);
        }
    }
    
    file2.close();
}

int prob_calculator::BaseHaplo(vector<pair<vector<int>&,vector<int>&>> &all_moles, float *posi_label, float *nega_label) 
{
    float total[(int)all_moles.size()];
    map<int, int> total_counts;
    vector<map<int, pair<int, int>>> all_label_counts;
    
    for (int i=0;i<all_moles.size();++i)
    {
        
        auto &allele_labels = all_moles[i];
        total[i]=0;
        map<int, pair<int, int>> label_counts;
        
        for (int mole: allele_labels.first)
        {
            total[i]+=allele_labels.first.size();
            
            for (int label: mole_phase[mole])
            {
                if (!label) continue;
                
                int label_abs = abs(label);
                total_counts[label_abs]++;
                
                if (label>0)
                {
                    label_counts[label_abs].first++;
                }
                else
                {
                    label_counts[label_abs].second++;
                }
            }
            
        }
        
        for (int mole: allele_labels.second)
        {
            total[i]+=allele_labels.second.size();
            
            for (int label: mole_phase[mole])
            {
                if (!label) continue;
                
                int label_abs = abs(label);
                total_counts[label_abs]++;
                
                if (label>0)
                {
                    label_counts[label_abs].first++;
                }
                else
                {
                    label_counts[label_abs].second++;
                }
            }
        }
        
        all_label_counts.push_back(label_counts);
    }
    
    int max_label=0, max_count=0;
    for (auto label: total_counts)
    {
        if (label.second>max_count)
        {
            max_count=label.second;
            max_label=label.first;
        }
    }
    
    if (!max_label || max_count<5 )
    {
        return 0;
    }
    
    for (int i=0; i<all_moles.size(); ++i)
    {
        posi_label[i]=all_label_counts[i][max_label].first;
        nega_label[i]=all_label_counts[i][max_label].second;
    }

    
    return max_label;
}

float prob_calculator::AlleleHaplo(int allele_index, int *allele_count, float *posi_label, float *nega_label) 
{

    float first_posi = posi_label[allele_index];
    float first_nega = nega_label[allele_index];
    
    float *false_negas, phase_confi;
    if (first_posi>first_nega)
    {
        false_negas=posi_label;
        phase_confi=first_posi;
    }
    else
    {
        false_negas=nega_label;
        phase_confi=first_nega;
    }
    
    float max_false=0.0;int max_allele=0;
    for (int i=0;i<10;++i)
    {
        
        if (i!=allele_index && max_false<false_negas[i])
        {
            max_false=false_negas[i];
            max_allele=i;
        }
        
    }
    
    float score=0;
    if (phase_confi>2*max_false || phase_confi<max_false)
    {
        score = phase_confi-2*max_false;
    }
    
    max_false*=(allele_count[max_allele]+1)/(posi_label[max_allele]+nega_label[max_allele]+1);
    max_false/=(allele_count[allele_index]+1)/(posi_label[allele_index]+nega_label[allele_index]+1);
    
    float score_corr=0;
    if (phase_confi>2*max_false || phase_confi<max_false)
    {
        score_corr = phase_confi-2*max_false;
    }
                
    return confidence*fmax(score_corr, score);
    
}



int prob_calculator::BaseQualfilter(int f_index, int r_index, float f_qual, float r_qual)
{
    
    float f_cutoff=alleles[f_index]->lowqual;
    float r_cutoff=alleles[r_index]->lowqual;
    
    
    return (f_qual>=f_cutoff) + (r_qual>=r_cutoff);
    
}


float prob_calculator::systemfilter(int first_alleleindex, int second_alleleindex, int first_obs, int second_obs, int &f_freq_cut, int &r_freq_cut, int &f_high_cut, int &r_high_cut)
{
    
    float freq_score=0.0;
    if (first_alleleindex<1024)
    {
        f_freq_cut=mutations[first_alleleindex]->errorcut;
        f_high_cut=mutations[first_alleleindex]->highcut;
    }
    else
    {
        f_freq_cut=mean_coverage/8;
        f_high_cut=mean_coverage/4;
    }
    
    if (second_alleleindex<1024)
    {
        r_freq_cut=mutations[second_alleleindex]->errorcut;
        r_high_cut=mutations[second_alleleindex]->highcut;
    }
    else
    {
        r_freq_cut=mean_coverage/8;
        r_high_cut=mean_coverage/4;
    }
    
    if (first_obs>=f_high_cut)
    {
        freq_score+=10;
    }
    else if (first_obs>=f_freq_cut)
    {
        freq_score+=7;
    }
    
    if (second_obs>=r_high_cut)
    {
        freq_score+=10;
    }
    else if (second_obs>=r_freq_cut)
    {
        freq_score+=7;
    }
    
    return  freq_score;
}


float  prob_calculator::factorial_log(int end)
{
    
    static float results_storage[1000];
    static bool results_index[1000]={0};
    
    if (end<1000 && results_index[end])
    {
        return results_storage[end];
    }
    
    float factorial=0.0;
    
    for (int i=2;i<end+1;i++)
    {
        factorial+=log2f(i);
    }
    
    
    if (end<1000)
    {
        results_index[end]=1;
        results_storage[end]= factorial;
    }
    
    return  factorial;
}

float  prob_calculator::permutation_log(int start,int end)
{
    return factorial_log(end)-factorial_log(start);
}

float  prob_calculator::log_add(float x, float y)
{
    
    static int log2addtable[]={100, 100, 99, 99, 98, 98, 97, 97, 96, 96, 95, 95, 94, 94, 93, 93, 92, 92, 91, 91, 90, 90, 89, 89, 88, 88, 88, 87, 87, 86, 86, 85, 85, 84, 84, 84, 83, 83, 82, 82, 81, 81, 81, 80, 80, 79, 79, 78, 78, 78, 77, 77, 76, 76, 76, 75, 75, 74, 74, 73, 73, 73, 72, 72, 72, 71, 71, 70, 70, 70, 69, 69, 68, 68, 68, 67, 67, 67, 66, 66, 65, 65, 65, 64, 64, 64, 63, 63, 63, 62, 62, 62, 61, 61, 61, 60, 60, 60, 59, 59, 58, 58, 58, 58, 57, 57, 57, 56, 56, 56, 55, 55, 55, 54, 54, 54, 53, 53, 53, 52, 52, 52, 52, 51, 51, 51, 50, 50, 50, 49, 49, 49, 49, 48, 48, 48, 47, 47, 47, 47, 46, 46, 46, 46, 45, 45, 45, 44, 44, 44, 44, 43, 43, 43, 43, 42, 42, 42, 42, 41, 41, 41, 41, 40, 40, 40, 40, 39, 39, 39, 39, 38, 38, 38, 38, 38, 37, 37, 37, 37, 36, 36, 36, 36, 36, 35, 35, 35, 35, 34, 34, 34, 34, 34, 33, 33, 33, 33, 33, 32, 32, 32, 32, 32, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29, 29, 29, 29, 29, 29, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 26, 26, 26, 26, 26, 26, 25, 25, 25, 25, 25, 25, 24, 24, 24, 24, 24, 24, 23, 23, 23, 23, 23, 23, 23, 22, 22, 22, 22, 22, 22, 22, 21, 21, 21, 21, 21, 21, 21, 21, 20, 20, 20, 20, 20, 20, 20, 19, 19, 19, 19, 19, 19, 19, 19, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    int x_100 = (int) (x*100);
    int y_100 = (int) (y*100);
    
    int higher=max(x_100, y_100);
    int lower = min(x_100, y_100);
    
    int diff=higher-lower;
    
    if (diff>800)
    {
        return lower/100.0;
        
    }
    
    return 1.0*(lower-log2addtable[diff])/100.0;
    
}

double prob_calculator::bionomial(int total, int obs)
{
    
    int notobs=total-obs;
    static double pro_log=log2f(0.5);
    static double unpro_log=log2f(0.5);
    
    double result= obs* pro_log+notobs* unpro_log + factorial_log(total)-factorial_log(notobs)- factorial_log(obs);
    return -result;
}


float prob_calculator::cul_bionom(int total, int obs /*minor obs*/)
{
    
    static float results_storage[51][51];
    static bool results_history[51][51];
    
    if (total>50)
    {
        obs=obs*50/total;
        total=50;
    }
    
    if (results_history[total][obs])
    {
        return results_storage[total][obs];
    }
    
    float result=1000.0;
    for (int i=0;i<=obs;i++)
    {
        if (2*i==total) {result=log_add(result, bionomial(total, i));}
        else {result=log_add(result, bionomial(total, i)-1);}
        
    }
    
    result=(result<0.0)?0.0:result;
    
    if (total<51 && obs<51)
    {
        results_storage[total][obs]=result;
        results_history[total][obs]=1;
    }
    return result;
}


double prob_calculator::bionomial_error(int total, int obs)
{
    
    int notobs=total-obs;
    static double pro_log=log2f(1.00-(float)pow(2,-confidence));
    static double unpro_log=-confidence;
    
    double result= obs* pro_log+notobs* unpro_log +factorial_log(total)-factorial_log(notobs)- factorial_log(obs);
    
    return -result;
}


float prob_calculator::cul_bionom_error(int total, int obs /*minor obs*/)
{
    
    static float results_storage[51][51];
    static bool results_history[51][51];
    
    if (total>50)
    {
        obs=obs*50/total;
        total=50;
    }
    
    if (results_history[total][obs]) return results_storage[total][obs];
    
    obs=total-obs;
    
    float result=1000.0;
    for (int i=0;i<=obs;i++)
    {
        result=log_add(result, bionomial_error(total, i));
    }
    
    result=(result<0.0)?0.0:result;
    
    if (total<51 && obs<51)
    {
        results_storage[total][obs]=result;
        results_history[total][obs]=1;
    }
    return result;
}

float prob_calculator::RandomError(int total, int obs)
{
    return cul_bionom_error(total,obs)-cul_bionom(total,obs);
}
