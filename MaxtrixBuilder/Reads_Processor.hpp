//
//  Reads_Processor.hpp
//  testcpp
//
//  Created by Wangfei MA on 8/22/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef Reads_Processor_hpp
#define Reads_Processor_hpp

#include "Config.h"
#include "realign.hpp"
#include "Cigar.hpp"
#include "Base_buffer.hpp"


class Processor{
    
public:
    
    Processor(){};
    ~Processor(){free(error_flags);}
    
    template <class Loader>  
    static void Initiate(Processor*, Loader* Readloader, int id);
    
    template <class Loader>
    int Process_get(Loader* Readloader);
    
    int Process_read(int );
    
    Cigar_Reader cigar;
    
    BamTools::BamAlignment alignment;
    
    int barindex;
    
    int pid;
    
protected:
    
    int Read_Alignment(const char* seq, const char* qual, const Cigar_Reader&, int bar_index, int alignment_start);
    
    int Read_MDtag(const char* seq, const char* qual, const Cigar_Reader&, int bar_index, int alignment_start);
    
    int Realign(int alignment_start);
        
    void Pause();
    
    void checkstatus(int);

    char loader_type;
            
    Realigner aligner;
        
    signed char* error_flags=(signed char*)malloc(max_Base_Buffer*sizeof(signed char));
    
    int flagsize=max_Base_Buffer;
    
};

template <class Loader>
void Processor::Initiate(Processor* myself, Loader* Readloader, int i) {
    
    Threads_Posi[i+2]=1;
    
    myself->pid=i+2;
    
    int if_reachend=0;
    while(!if_reachend){
        
        if_reachend=myself->Process_read(myself->Process_get(Readloader));
    }
    
}


template <class Loader>
int Processor::Process_get(Loader* Readloader){
    
    int ifget=Readloader->LoadAlignment(alignment, cigar,  &barindex, pid);
    
    return ifget; //Get next cycle Quit the thread
}


#endif /* Reads_Processor_hpp */
