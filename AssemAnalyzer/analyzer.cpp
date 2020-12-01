//
//  main.cpp
//  snp_caller
//
//  Created by Wangfei MA on 6/30/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "sampling.hpp"

#define  NucleoIndex( base , index) for ( index=0;index<8;++index) {if ("ATCGDNMI"[index]==base) break;}


using namespace std;

extern float coverage;
float coverage=0;

extern float coverage_raw;
float coverage_raw=0;

extern float pre_errors[8];
float pre_errors[8]={0};

extern float pre_qualcuts[8];
float pre_qualcuts[8]={0.0};

extern int turn;
int turn=0;

extern long long total_alleles[8];
long long total_alleles[8]={0};

extern bool ifMD;
bool ifMD=0;

extern int sign;
int sign=1;

using namespace std;

void readinfor(string statsfile){
    
    cout<<statsfile<<endl;
    
    fstream file;
    file.open(statsfile);
    string line;
    
    vector<long long> counters;
    
    getline(file,line);
    
    int i=0;
    int pos1, allele_index;
    while(getline(file,line))
    {
        
        i++;
        if (i<10 && line.size())
        {
           counters.push_back((long long) atoll(line.c_str()));
        }
        else if ((pos1 = (int)line.find("lowcut:",2)) != std::string::npos)
        {
            NucleoIndex(line[0],allele_index);
            pre_qualcuts[allele_index] = atof(line.c_str()+pos1+8);
            turn=1;
        }
        
    }
    
    long long allelecounts=accumulate(counters.begin()+2, counters.end(), 0L);
    
    for (int i=0; i<counters.size()-2 && i<8;i++)
    {
        
        total_alleles[i]=counters[i+2]*counters[0]/counters[1];
        
        pre_errors[i]=log(((float)counters[i+2]+1)/counters[1]);
                
    }
    
    coverage=(float)(allelecounts)/counters[0];
    
    coverage_raw=(float)(counters[1])/counters[0];
    
    file.close();
    
    
    
}



int main(int argc, const char * argv[]) {
    
    cout << "Starting Program" <<endl;
    const char* inputfile=argv[1];
    const char* reffile=argv[2];
    const char* chr=argv[3];
    
    sign=atoi(argv[4]);
    
    string samplefile=(string(inputfile)+"_sampling");
    string statsfile=(string(inputfile)+"_stat");
        
    readinfor(statsfile);
    
    sampling SamplePool;
    
    SamplePool.loadref(reffile,chr);
    
    SamplePool.readsamples(samplefile);
            
    SamplePool.statistics();
    
    SamplePool.output(statsfile);

    return 0;
}
