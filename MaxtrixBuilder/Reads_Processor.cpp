//
//  Reads_Processor.cpp
//  testcpp
//
//  Created by Wangfei MA on 8/22/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "Reads_Processor.hpp"

extern Buff_ctrl Base_Buffer_ctrl;

using namespace std;

static void Flag_Errors(const char *seq, signed char *des, int len){
    
    char current_char='X';
    int current_index;
    
    int ATCG_start[4]={1,1,1,1};
    int ATCG_last[4]={1,1,1,1};
    int ATCG_counter[4]={0,0,0,0};
    
    for (int i=0;i<len;++i){
        
        current_char=seq[i];
        
        NucleoIndex(toupper(current_char), current_index)
        
        if (current_index>3) continue;
        
        if (i-ATCG_last[current_index]>2)
        {
            
            if (ATCG_counter[current_index]>3)
            {
                memset(des+ATCG_start[current_index]-1,(signed char)-1, (ATCG_last[current_index]-ATCG_start[current_index]+3)*sizeof(signed char));
                
            }
            
            ATCG_start[current_index]=i;
            ATCG_counter[current_index]=0;
        }
        
        ATCG_last[current_index]=i;
        ++ATCG_counter[current_index];
        
    }
    
    
    for ( int current_index=0;current_index<4;++current_index){
        
        if (ATCG_counter[current_index]>3){
            
            memset(des+ATCG_start[current_index]-1,(signed char)-1,(ATCG_last[current_index]-ATCG_start[current_index]+2)*sizeof(signed char));
            
        }
    }
}



static char Mean_Qual(const char qual[], int charsize, int start, int end){
    
    end=(end<charsize-1)?end:(charsize-1);
    if (end<=start) return ' ';
    
    int qualsum=0;
    for (int i=start;i<end;++i){
        
        qualsum+=qual[i]-'!';
    }
    qualsum/=(end-start);
    
    return '!'+qualsum;
    
}

static short Insert_Hash(const char *seq, int size){
    
    
    bool use_size=0;
    char base;
    short hash=1;
    if (size<5){
        
        for (int i=size-1; i>=0;--i){
            
            hash<<=2;
            base=seq[i];
            
            if (base=='A') {
                hash+=0;
            }
            
            else if (base=='T'){
                hash+=1;
            }
            
            else if (base=='C'){
                hash+=2;
            }
            
            else if (base=='G'){
                hash+=3;
            }
            
            else{
                use_size=1;
            }
        }
    }
    else {
        
        for (int i=size-1; i>size-3;--i){
            
            hash<<=2;
            base=seq[i];
            
            if (base=='A') {
                hash+=0;
            }
            
            else if (base=='T'){
                hash+=1;
            }
            
            else if (base=='C'){
                hash+=2;
            }
            
            else if (base=='G'){
                hash+=3;
            }
            
            else{
                use_size=1;
            }
        }
        
        for (int i=2; i>=0;--i){
            
            hash<<=2;
            base=seq[i];
            if (base=='A') {
                hash+=0;
            }
            
            else if (base=='T'){
                hash+=1;
            }
            
            else if (base=='C'){
                hash+=2;
            }
            
            else if (base=='G'){
                hash+=3;
            }
            
            else{
                use_size=1;
            }
        }
        
        
    }
    
    
    if (use_size) {return min(9,size);}
    else {return hash+10;}
    
}



void Processor::Pause(){
    
    int alignment_end=alignment.Position+ alignment.Length;
    
    write_cd.wait(writing_lk, [&] {return(buff_start > alignment_end-max_Base_Buffer);} );
    
}


int Processor::Process_read(int ifget){
    
    
    if (alignment.Position+ alignment.Length-buff_start>max_Base_Buffer) Pause();
    
    if (ifget<=0) {
        
        return 1;
    }
    
    int alignment_start=alignment.Position;
    
    int alignment_size=alignment.Length;
    
    if (alignment.QueryBases.length()>2){
        
        alignment_size=Read_Alignment( alignment.QueryBases.c_str(), alignment.Qualities.c_str(), cigar, barindex, alignment_start);
    }
    
    else if (useMDtag){
        
        alignment_size=Read_MDtag( alignment.QueryBases.c_str(), alignment.Qualities.c_str(), cigar,barindex, alignment_start);
        
    }
    
    Realign(last_read_start);
        
    
    return 0; //Get next cycle Quit the thread
    
}



int Processor::Realign( int allalign_start){
    
    std::shared_lock<std::shared_timed_mutex> Read_Lock(buffer_move_lock);
    
    Threads_Posi_lock.lock();
    
    int l_buff_start=buff_start;

    int get_last_realign=last_realign_end;
    
    int realign_start=max(l_buff_start,get_last_realign);
    
    int realign_end=max(realign_start, allalign_start-100);
    
    last_realign_end=realign_end;
    
    Threads_Posi_lock.unlock();
    
    int valide_realign_end=min(int(buff_end),allalign_start-100);
    
    if (realign_end-realign_start>0) {
        
        aligner.Buffer_Region(l_buff_start,allalign_start);
        
        for (int current_index=realign_start; current_index<valide_realign_end; ++current_index){
            
            aligner.Realign(current_index);
        }
        
    }
    
    return 0;

}



int Processor::Read_Alignment(const char *seq, const char* qual,const Cigar_Reader& cigar_data, int bar_index, int alignment_start){
    
    int qualsize=(int)strlen(qual);
    
    if (qualsize<2) {qualsize=0;}
    
    int num_cigar=cigar_data.cigar_counter;
    
    std::vector<int> cigar_segments=cigar_data.cigar_segments;
    
    std::vector<char> cigar_types=cigar_data.cigar_types;
    
    int current_segment=0;
    char current_type, current_qual;
    
    int current_qposi=0;
    
    int current_rposi=alignment_start;
    
    int Loopend=0;
    
    int basetype = 0;
    
    int asize=0,error_flag;
    
    if (flagsize<=qualsize+4) {
        error_flags=(signed char*)realloc(error_flags,(qualsize+4)*sizeof(signed char));
        flagsize=qualsize+4;
    };
    
    memset(error_flags,0,qualsize*sizeof(signed char));
    
    Flag_Errors(seq, error_flags, qualsize);
    
    for (int i=0; i< num_cigar; ++i){
        
        
        current_segment=cigar_segments[i];
        current_type=cigar_types[i];
        
        if (current_type=='H'){continue;}
        
        else if (current_type=='S'){
            
            current_qposi+=current_segment;
            
            continue;
            
        }
        
        else if (current_type=='N'){
            
            current_qposi+=current_segment;
            
            current_rposi+=current_segment;
            
            continue;
        }
        
        if (current_type=='I'){
            
            error_flag=(int)(*max_element(&error_flags[current_qposi], &error_flags[current_qposi]+current_segment))*2+1;
            
            if (current_qposi+current_segment<qualsize)   current_qual=Mean_Qual(qual, qualsize, current_qposi, current_qposi+2);
            else {current_qual='!';}

            Base_Buffer_ctrl.update(current_rposi, 7, bar_index, Mean_Qual(qual, qualsize, current_qposi, current_qposi+current_segment),Insert_Hash(seq+current_qposi, current_segment));

            current_qposi+=current_segment;
            ++current_rposi;
            
        }
        
        else if (current_type=='D'){
            
            if (current_qposi<qualsize)   current_qual=Mean_Qual(qual, qualsize, current_qposi, current_qposi+2);
            else {current_qual='!';}
            
            error_flag=error_flags[current_qposi]*2+1;
            
            intpair deletion=make_pair(current_rposi,current_rposi+current_segment);
            
            Loopend=current_rposi+current_segment;
            for (;current_rposi<Loopend;++current_rposi){
                
                Base_Buffer_ctrl.update(current_rposi, 4, error_flag*bar_index, current_qual,deletion);
            }
            
        }
        
        else if (current_type=='M'){
            
            asize+=current_segment;
            
            Loopend=current_rposi+current_segment;
            for (;current_rposi<Loopend ;++current_rposi) {
                
                error_flag=error_flags[current_qposi]*2+1;
                
                NucleoIndex(seq[current_qposi],basetype);
                
                Base_Buffer_ctrl.update(current_rposi, basetype, error_flag*bar_index, qual[current_qposi],0);
                
                ++current_qposi;
                
            }
            
        }
        
    }
    
    return asize;
}

int Processor::Read_MDtag(const char* seq, const char* qual, const Cigar_Reader& cigar_data, int bar_index, int alignment_start){
    
    int qualsize=(int)strlen(qual);
    if (qualsize<2) {qualsize=0;}
    
    int num_cigar=cigar_data.cigar_counter;
    
    std::vector<int> cigar_segments=cigar_data.cigar_segments;
    
    std::vector<char> cigar_types=cigar_data.cigar_types;
    
    int current_segment=0;
    char current_type, current_qual;
    
    int current_rposi=alignment_start;
    
    int current_qposi=0;
    
    int Loopend=0;
    
    int basetype = 0;
    
    int insertsize=0,asize=0;
    
    short insert;
    
    for (int i=0; i< num_cigar; ++i){
        
        current_segment=cigar_segments[i];
        current_type=cigar_types[i];
        
        if (current_type=='H'){continue;}
        
        else if (current_type=='S'){
            
            current_qposi+=current_segment;
            continue;
            
        }
        
        else if (current_type=='N'){
            
            current_qposi+=current_segment;
            
            current_rposi+=current_segment;
            
            continue;
        }
        
        
        if (current_type=='I'){
            
            insertsize=(current_segment>9)?9:current_segment;
            
            insert=insertsize;
            
            Base_Buffer_ctrl.update(current_rposi, 7, bar_index, Mean_Qual(qual, qualsize, current_qposi, current_qposi+current_segment),insert);
            
            ++current_rposi;
            current_qposi+=current_segment;
            
        }
        
        else if (current_type=='D'){
            
            current_qual=Mean_Qual(qual, qualsize, current_qposi, current_qposi+2);
            
            intpair deletion(current_rposi,current_rposi+current_segment);
            
            Loopend=current_rposi+current_segment;
            for (;current_rposi<Loopend;++current_rposi){
                
                Base_Buffer_ctrl.update(current_rposi, 4, bar_index,current_qual,deletion);
                
            }
                        
        }
        
        else if (current_type=='M'){
            
            asize+=current_segment;
            
            Loopend=current_rposi+current_segment;
            
            for (;current_rposi<Loopend ;++current_rposi) {
                
                if (current_rposi<refgenome.size())NucleoIndex(seq[current_rposi],basetype);
                
                Base_Buffer_ctrl.update(current_rposi,  basetype, bar_index, qual[current_qposi],0);
                
                ++current_qposi;
                
            }
            
        }
    
    }
    
    return asize;
}



