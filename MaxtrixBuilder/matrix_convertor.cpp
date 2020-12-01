//
//  matrix_converter.cpp
//  bam_matrix_builder
//
//  Created by Wangfei MA on 6/17/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>

#include "Config.h"
#include "Bamloader.hpp"
#include "Samloader.hpp"
#include "Base_buffer.hpp"

using namespace std;


//global variables for options
extern std::atomic_uchar coverage;
std::atomic_uchar coverage(50);

extern int writing_mode;
int writing_mode=1;

extern int qualcut;
int qualcut=20;

extern bool read_nobarcode;
bool read_nobarcode=1;

extern bool ifcheckflag;
bool ifcheckflag=0;

extern volatile bool useMDtag;
volatile bool useMDtag=1;


//global variables for multi-threads
extern volatile int nthreads;
int volatile nthreads=1;

extern std::atomic_int Threads_Posi[256];
std::atomic_int Threads_Posi[256];

extern std::atomic_int Threads_Next_Posi[256];
std::atomic_int Threads_Next_Posi[256];

std::mutex file_read_lock;
std::shared_timed_mutex buffer_move_lock;


std::mutex Threads_Posi_lock;

std::mutex Minilock_lock;
std::mutex cell_init_lock;

std::mutex realign_record_lock;
std::mutex barcode_record_lock;

extern std::condition_variable write_cd;
std::condition_variable write_cd;

extern std::mutex writing_lock;
std::mutex writing_lock;

extern std::unique_lock<std::mutex> writing_lk;
std::unique_lock<std::mutex> writing_lk(writing_lock);

//global variables for buffer data
extern std::string refgenome;
std::string refgenome;

extern std::unordered_map<int, std::unordered_map<int, std::intpair>> Corr_history;
std::unordered_map<int, std::unordered_map<int, std::intpair>> Corr_history;

extern Buff_ctrl Base_Buffer_ctrl;
Buff_ctrl Base_Buffer_ctrl;

extern volatile int buff_start;
volatile int buff_start=1;

extern volatile int buff_end;
volatile int buff_end=1;

extern volatile int last_read_start;
volatile int last_read_start=1;

extern volatile int last_realign_end;
volatile int last_realign_end=1;
//global variables for statistics
extern long long total_counts[10];
long long total_counts[10]={0,0,0,0,0,0,0,0,0,0};


int main(int argc, const char * argv[]) {
    
    cout << "Starting Program" <<endl;
    const char *inputfile=argv[1];
    const char *outputfile=argv[2];
    const char *reffile=argv[3];

    const char *input_chr="";
    if (argc>4){
        input_chr=argv[4];
    }
    
    int input_start=-1;
    if (argc>5){
        input_start=(int) atoi(argv[5])-1;
    }
    
    int input_end=INT_MAX;
    if (argc>6){
        input_end=(int) atoi(argv[6])+1;
    }
    
    if (argc>7){
        nthreads=(int) atoi(argv[7]);
    }

    if (argc>8){
        writing_mode=(int) atoi(argv[8]);
    }
    
    if (argc>9){
        qualcut=(int) atoi(argv[9]);
    }
    
    printf("Converting %s to %s at region %s:%d-%d\n", inputfile,outputfile,input_chr,input_start,input_end);
    
	long long barcode_counter;
    if ( string(inputfile).find( ".bam\0" ) == string(inputfile).size()-4 )
    {
        cout<<"Reading Bam file"<<endl;
        
        
        Bamloader Bamloader(inputfile, outputfile, reffile);
        
        barcode_counter=Bamloader.Process_region(input_chr, input_start, input_end);
           
    }
    else {
        
        cout<<"Reading Sam file\n"<<endl;
        
        Samloader Samloader(inputfile, outputfile, reffile);
        
        barcode_counter=Samloader.Process_region(input_chr, input_start, input_end); 
        
    }
    
    ofstream finfor;
    
    finfor.open (string(outputfile)+"_stat",std::ofstream::out);
    
    finfor<<"global_stats:\n";
    for (int i=0;i<10;++i){
        
        finfor<<total_counts[i]<<"\n";
    }
    finfor<<barcode_counter;
    
    finfor.close();
    
    cout << "Program Finished" <<endl;

    return 0;
}
