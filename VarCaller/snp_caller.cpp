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
#include <algorithm>
#include "call_snp.hpp"
#include "probability_calculator.hpp"

using namespace std;

extern float coverage;
float coverage=0;

extern float coverage_raw;
float coverage_raw=0;

extern float pre_errors[8];
float pre_errors[8]={0};

extern int base_qualcuts[8];
int base_qualcuts[8]={0};

extern int base_hqualcuts[8];
int base_hqualcuts[8]={0};

extern int confidence;
int confidence=4;

extern int icycle;
int icycle=0;

int main(int argc, const char * argv[]) {
    
    cout << "Starting Program" <<endl;
    const char* inputfile=argv[1];
    
    const char* outputfile=argv[2];
    
    const char* reffile=argv[3];
    
    const char* chr=argv[4];
    
    if (argc>5) confidence=atoi(argv[5]);
    
    prob_calculator calculator;
    
    const char* snpfile="";
    
    if (argc>6)
    {
        snpfile=argv[6];
        
        calculator.loadSNP(snpfile);
    }
    
    const char* molefile="";
    if (argc>7)
    {
        molefile=argv[7];
        
        calculator.loadMole(molefile);
        
        icycle=1;
    }
    
    calculator.loadStat(string(inputfile)+"_stat");
    
    snp_caller snpcaller(calculator);
    
    snpcaller.snp_calling(inputfile, outputfile, reffile, chr);
        
    cout << "Program Finished" <<endl;
    return 0;
}
