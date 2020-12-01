//
//  Base_buffer.cpp
//  testcpp
//
//  Created by Wangfei MA on 8/26/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "Base_buffer.hpp"

using namespace std;



void Buff_ctrl::Mini_unLock(int index){
    
    lock_guard<mutex>guard(Minilock_lock);
    
    Base_IfInit[index]=1;
}


void Buff_ctrl::Mini_Lock(int index){
    
    atomic_uchar *lock=&Base_IfInit[index];
    
    do {
        
        if (Minilock_lock.try_lock()){
            
            if (*lock>(unsigned char)1) {
                
                Minilock_lock.unlock();
            }
            
            else if (!*lock){
                
                *lock=2;
                Minilock_lock.unlock();
                //InitBase(index);
                return;
            }
            
            else{
                *lock=2;
                Minilock_lock.unlock();
                return;
            }
            
        }
        
        this_thread::sleep_for(chrono::microseconds(1));
        
    }while (1);
    
}


Base* Buff_ctrl::Lock_Base(int index){
    
    std::shared_lock<std::shared_timed_mutex> Read_Lock(buffer_move_lock);
    
    Mini_Lock(index);
    
    return Base_Buffer[index-buff_start];
    
}


void Buff_ctrl::ReleaseAllBases(int release_size){
    
    int release_index=0;
    for(release_index=0; release_index<release_size; ++release_index)
    {
        Release(release_index);
        
    }
    
}


void Buff_ctrl::move(int movesize){
    
    std::unique_lock<std::shared_timed_mutex> Move_Lock(buffer_move_lock);
    
    ReleaseAllBases(movesize);
    
    Base* temp;
    for (int index=movesize; index<buff_end-buff_start; ++index){
        if (!(unsigned char)Base_IfInit[index]){continue;}
        
        
        
        temp=Base_Buffer[index-movesize];
        
        Base_Buffer[index-movesize]=Base_Buffer[index];
        
        Base_Buffer[index]=temp;
        
        Base_IfInit[index-movesize]=1;
        
        Base_IfInit[index]=0;
        
    }
    
}

void Buff_ctrl::update(int init_index, int column, int barindex, char qual, int extra_data){
    
    if (column==6 && init_index<refgenome.length()) {NucleoIndex(toupper(refgenome[init_index]), column)}
    
    if (column>7) column=5;
    
    Base * Base_Data=Lock_Base(init_index);
    
    int current_columnnum=Base_Data->lens[column]++;
    
    if (!current_columnnum){
        
        if (column<5) {
            
            Reserve_alloc(Base_Data->snps[column], int, Reserve_size, blocksize)
            
            Reserve_alloc(Base_Data->quals[column], char, Reserve_size, blocksize)
            
        }
        
        else {
            
            Base_Data->snps[column]=(int*)malloc(blocksize*sizeof(int));
            
            Base_Data->quals[column]=(char*)malloc(blocksize*sizeof(char));
            
        }
        
    }
    
    else if(!((current_columnnum-blocksize)%coverage)){
        
        Base_Data->snps[column]=(int*)realloc(Base_Data->snps[column],(coverage+current_columnnum+1)*sizeof(int));
        
        Base_Data->quals[column]=(char*)realloc(Base_Data->quals[column],(coverage+current_columnnum+1)*sizeof(char));
    }
    
    Base_Data->snps[column][current_columnnum]=barindex;
    
    Base_Data->quals[column][current_columnnum]=qual;
    
    Mini_unLock(init_index);
    
}

void Buff_ctrl::update(int init_index, int column, int barindex, char qual, intpair extra_data){
    
    Base * Base_Data=Lock_Base(init_index);
    
    int current_columnnum=Base_Data->lens[4]++;
    
    if (!current_columnnum){
        
        Reserve_alloc(Base_Data->snps[4], int, Reserve_size, blocksize)
        
        Reserve_alloc(Base_Data->quals[4], char, Reserve_size, blocksize)
        
        Reserve_alloc( Base_Data->deletions, intpair, Reserve_size, blocksize)
        
    }
    
    else if(!((current_columnnum-blocksize)%coverage)){
        
        Base_Data->snps[4]=(int*)realloc(Base_Data->snps[4],(coverage+current_columnnum+1)*sizeof(int));
        
        Base_Data->quals[4]=(char*)realloc(Base_Data->quals[4],(coverage+current_columnnum+1)*sizeof(char));
        
        Base_Data->deletions=(intpair*) realloc(Base_Data->deletions,(coverage+current_columnnum+1)*sizeof(intpair));
        
    }
    
    Base_Data->snps[4][current_columnnum]=barindex;
    
    Base_Data->quals[4][current_columnnum]=qual;
    
    Base_Data->deletions[current_columnnum]=extra_data;
    
    Mini_unLock(init_index);
    

}

void Buff_ctrl::update(int init_index,int column, int barindex, char qual, short extra_data){
    
    Base * Base_Data=Lock_Base(init_index);
    
    int current_columnnum=Base_Data->lens[7]++;
    
    if (!current_columnnum){
        
        Base_Data->snps[7]=(int*)malloc(blocksize*sizeof(int));
        
        Base_Data->quals[7]=(char*)malloc(blocksize*sizeof(char));
        
        Base_Data->inserts=(short*)malloc(blocksize*sizeof(short));
        
    }
    
    else if(!((current_columnnum-blocksize)%coverage)){
        
        Base_Data->snps[7]=(int*)realloc(Base_Data->snps[7],(coverage+current_columnnum+1)*sizeof(int));
        
        Base_Data->quals[7]=(char*)realloc(Base_Data->quals[7],(coverage+current_columnnum+1)*sizeof(char));
        
        Base_Data->inserts=(short*)realloc(Base_Data->inserts,(coverage+current_columnnum+1)*sizeof(short));
        
    }
    
    Base_Data->snps[7][current_columnnum]=barindex;
    
    Base_Data->quals[7][current_columnnum]=qual;
    
    Base_Data->inserts[current_columnnum]=extra_data;
    
    Mini_unLock(init_index);
    
}



void Buff_ctrl::Release(int release_index){
    
    if (!(unsigned char)Base_IfInit[release_index]) {
        
        return ;
    }
    
    Base * current_base=Base_Buffer[release_index];
    
    if (current_base->lens[7]) {
        
        free(current_base->inserts);
    }
    
    if (current_base->lens[4]) {
        
        free(current_base->deletions);
    }
    
    for (int j=0;j<8; ++j){
        
        if (current_base->lens[j]){
            
            free(current_base->snps[j]);
            free(current_base->quals[j]);
        }
        
    }
    
    memset(current_base->lens,(int)0,8*sizeof(int));
    
    Base_IfInit[release_index]=0;
    
    
}


void Buff_ctrl::resize(Base* move_to, int target_allele, int addsize){ 
    
    int current_size=move_to->lens[target_allele], current_block_num=(current_size-blocksize)/coverage+1;
    
    int new_size=current_size+addsize, new_block_num=(new_size-blocksize)/coverage+1;
    
    if(new_block_num!=current_block_num){
        
        int realloc_size=(coverage*new_block_num+1+blocksize);
        
        move_to->snps[target_allele]= (int*)realloc(move_to->snps[target_allele],realloc_size*sizeof(int));
        
        move_to->quals[target_allele]= (char*)realloc(move_to->quals[target_allele],realloc_size*sizeof(char));
        
        if (target_allele==7) {move_to->inserts=(short*)realloc(move_to->inserts,realloc_size*sizeof(short));}
        
        if (target_allele==4) {move_to->deletions=(intpair*)realloc(move_to->deletions,realloc_size*sizeof(intpair));}
    }
    
    move_to->lens[target_allele]=new_size;
    
}


//correction and move data
void Buff_ctrl::transfer(int old_pos, int new_pos, int target_allele, std::vector<int>& moles){
    
    Base* move_from=Base_Buffer[old_pos-buff_start], *move_to=Base_Buffer[new_pos-buff_start];
    
    //realloc the destinations
    resize(move_to,target_allele,(int) moles.size());
    
    int current_size=move_to->lens[target_allele];
    
    //move the data
    int num_alladdmoles=move_from->lens[target_allele];
    int* alladdmoles=move_from->snps[target_allele];
    char* allquals=move_from->quals[target_allele];
    short* allinserts=move_from->inserts;
    intpair* alldeletions=move_from->deletions;
    
    int addmole, addindex;
    for (int index=0;index<moles.size();++index){
        
        addmole=moles[index];
        
        addindex=(int)(find(alladdmoles,alladdmoles+num_alladdmoles, addmole)-alladdmoles);
        if (addindex==num_alladdmoles) addindex=(int)(find(alladdmoles,alladdmoles+num_alladdmoles, -addmole)-alladdmoles);
        
        move_to->snps[target_allele][current_size+index]=addmole;
        
        move_to->quals[target_allele][current_size+index]=allquals[addindex];
        allquals[addindex]=0;
        
        if (target_allele==7) move_to->inserts[current_size+index]=allinserts[addindex];
        
        if (target_allele==4) move_to->deletions[current_size+index]=alldeletions[addindex];
    }
    
    if (moles.size()>current_size){
        
        move_from->lens[target_allele]=0;
        
        free(allquals);
        free(alladdmoles);
        
        if (target_allele==7) {
            free(allinserts);
        }
        
        if (target_allele==4){
            
            free(alldeletions);
        }
    }
    
    else {
        
        
        int num_reserve=0, imole;
        for (int index=0;index<num_alladdmoles;++index){
            
            imole=allquals[index];
            
            if (!imole) continue;
            
            alladdmoles[num_reserve]=alladdmoles[index];
            allquals[num_reserve]=allquals[index];
            if (target_allele==7) allinserts[num_reserve]=allinserts[index];
            if (target_allele==4) alldeletions[num_reserve]=alldeletions[index];
            ++num_reserve;
        }
        
        move_from->lens[target_allele]=num_reserve;
        
    }
}
