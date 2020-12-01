//
//  Config.h
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef Config_h
#define Config_h

#include <iostream>
#include <stdio.h>
#include <limits.h>
#include <sstream>   
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <numeric>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <string>
#include <chrono>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include <atomic>


#define hash0 300
#define hash1 299

#define max_Base_Buffer 1001000
#define blocksize 10
#define Reserve_size (10+2*coverage)

#define  Reserve_alloc( prt , type , rsize , size ) prt=(type*) malloc(rsize*sizeof(type));  prt=(type*) realloc(prt, size*sizeof(type));
#define  Nucleotides( base ) "ATCGDNMI"[base]
#define  NucleoIndex( base , index) for ( index=0;index<8;++index) {if ("ATCGDNMI"[index]==base) break;}
#define  UNIQ( n, len )  ( (int)(unordered_set<unsigned int> ((n),(n)+(len))).size() )
#define  Vprint( n )  for ( auto Vprinti: n ) { cout<<Vprinti<<" ";} cout<<endl;
#define  Aprint( n , type ,len )  for ( type Aprinti: vector<type> (n, n+ len ) ) { cout<<Aprinti<<" ";} cout<<endl;
#define  intpair pair<int, int>
#define  maxvalue ( dict , max, func) for (auto pair: dict) {if (func(dict.second)>max) max=func(dict.second);}
#define  Wait(n) std::this_thread::sleep_for(std::chrono::nanoseconds(n));
#define  tid std::this_thread::get_id()


typedef struct
{
    int lens[8];
    int *snps[8];
    char *quals[8];
    short *inserts;
    std::pair<int, int>* deletions;
}Base;

typedef struct
{
    std::vector<int> up;
    std::vector<int> back;
    std::vector<int> var;
}Seed;

typedef struct
{
    long long barcode;
    char* barname;
    int firstposi;
    int lastposi;
    int molesize;
    size_t molecounter;
    
}Molecule;

extern std::atomic_uchar coverage;

extern int writing_mode;

extern int qualcut;

extern bool read_nobarcode;

extern bool ifcheckflag;

extern volatile bool useMDtag;

extern volatile int nthreads;

extern std::atomic_int Threads_Posi[256];
extern std::atomic_int Threads_Next_Posi[256];

extern long long total_counts[10];

extern Base * Base_Buffer[max_Base_Buffer];

extern std::atomic_uchar Base_IfInit[max_Base_Buffer];

extern std::unordered_map<int, std::unordered_map<int, std::intpair>> Corr_history;

extern std::string refgenome;

extern volatile int buff_start;

extern volatile int buff_end;

extern volatile int last_read_start;

extern volatile int last_realign_end;

extern std::condition_variable write_cd;
extern std::mutex writing_lock;
extern std::unique_lock<std::mutex> writing_lk;

extern std::mutex Threads_Posi_lock;
extern std::mutex Minilock_lock;;
extern std::mutex cell_init_lock;
extern std::mutex file_read_lock;
extern std::mutex barcode_record_lock;
extern std::mutex realign_record_lock;

extern std::shared_timed_mutex buffer_move_lock;

#endif /* Config_h */
