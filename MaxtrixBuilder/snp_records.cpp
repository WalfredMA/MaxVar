//
//  snp_records.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "snp_records.hpp"

using namespace std;

extern Base * Base_Buffer[max_Base_Buffer];
Base * Base_Buffer[max_Base_Buffer];

extern std::atomic_uchar Base_IfInit[max_Base_Buffer];
std::atomic_uchar Base_IfInit[max_Base_Buffer];

extern Buff_ctrl Base_Buffer_ctrl;

static bool Ifsnp(int total, int first, int second){
    
    if (4*second>max(8,total)) {
        
        return 1;
    }
    
    
    //if (calculator.cul_bionom(total, first)>5){
    
    //return 1;
    //}
    
    return 0;
}


static short Insert_Unhash(short inserts_hash, char insert[]){
    
    int insert_len=0;
    if (inserts_hash>=10){
        
        inserts_hash-=10;
        while(inserts_hash){
            
            insert[insert_len++]=Nucleotides(inserts_hash%4);
            inserts_hash>>=2;
        }
        --insert_len;
        insert[insert_len]='\0';
    }
    
    else{
        
        insert[0]='0'+inserts_hash;
        insert[1]='\0';
    }
    
    return insert_len;
}


Snp_Records::Snp_Records(const char* file):outfile(file), write_start(0), Linesize(10000000){
    
    for (int i=0;i<max_Base_Buffer;i++){
        
        Base_Buffer[i]=new Base;
        memset(Base_Buffer[i]->lens,0,8*sizeof(int));
        
        Base_IfInit[i]=0;
        
    }
        
    std::unique_lock<std::shared_timed_mutex> Move_Lock(buffer_move_lock,defer_lock);
    
    samplefile=std::string(file)+string("_sampling");

}

void Snp_Records::Count_SmallDel(Base* current_base, int* all_mole, int* good_mole, int* qual, int len){
    
    int* Molecules=current_base->snps[4];
    
    char* quals=current_base->quals[4];
    
    int Good_moles[len];
    int Bad_moles[len];
    
    int Good_moles_counter=0;
    int Bad_moles_counter=0;

    for (int j=0;j<len;++j){
        
        if (current_base->deletions[j].second-current_base->deletions[j].first<10){
            
            if (Molecules[j]>0){
                
                Good_moles[Good_moles_counter++]=Molecules[j];
                *qual+=(int)quals[j];
            }
            
            else {
                
                Bad_moles[Bad_moles_counter++]=Molecules[j];
            }
        }
    }
    
    if (Good_moles_counter) {*good_mole=(int)unordered_set<int>(Good_moles, Good_moles+Good_moles_counter).size();}
    else {*good_mole=0;}
    
    *all_mole=*good_mole;
    if (Bad_moles_counter) *all_mole+=(int)unordered_set<int>(Bad_moles, Bad_moles+Bad_moles_counter).size();
    
    if (Good_moles_counter) *qual=(*qual)*(*good_mole)/Good_moles_counter;
    
}

void Snp_Records::Count_Allele(int* Molecules, char* quals, int* all_mole, int* good_mole, int *qual, int len){
    
    int Good_moles[len];
    int Bad_moles[len];
    
    int Good_moles_counter=0;
    int Bad_moles_counter=0;
    
    for (int j=0; j<len;++j){
        
        if (Molecules[j]>0){
            
            Good_moles[Good_moles_counter++]=Molecules[j];
            *qual+=(int)quals[j];
        }
        
        else {
            
            Bad_moles[Bad_moles_counter++]=Molecules[j];
        }
    }
    
    if (Good_moles_counter) {*good_mole=(int)unordered_set<int>(Good_moles, Good_moles+Good_moles_counter).size();}
    else {*good_mole=0;}
    
    *all_mole=*good_mole;
    if (Bad_moles_counter) *all_mole+=(int)unordered_set<int>(Bad_moles, Bad_moles+Bad_moles_counter).size();
    
    if (Good_moles_counter) *qual=(*qual)*(*good_mole)/Good_moles_counter;
    
}

bool Snp_Records::Alelle_Stat(int current_index,FILE *sampling_file,int randnum)
{
    Base * current_base=Base_Buffer[current_index]; 
    
    int themax=0, thesecond=0, maxindex=0, total=0,errorcount=0;
    int allele_counts[8]={0}, sumquals[8]={0};
	int tempint, tempposi;

    for (int i=0; i<8; ++i){
       
        if (i==5){
            continue;
        }
        
        tempint=current_base->lens[i];
        
        tempposi=0;
        
        if (tempint){
            
            
            if (i!=4){
                
                Count_Allele(current_base->snps[i], current_base->quals[i], &tempint, &tempposi, &sumquals[i], tempint);
                
            }
            
            else{
                
                Count_SmallDel(current_base, &tempint, &tempposi, &sumquals[4], tempint);
                
            }
        }
        
        allele_counts[i]=tempposi;
        
        errorcount+=tempint-tempposi;
        total+=tempint;
        
        if (tempint>themax){
            maxindex=i;
            thesecond=themax;
            themax=tempint;
        }
        
        else if (tempint>thesecond){
            thesecond=tempint;
        }
    }
    
    if (total){
      
        total_counts[maxindex+2]+=themax;
        
        ++total_counts[0];
        
        total_counts[1]+=total;
                
        if(!randnum && errorcount<total/5){
            
            if (current_index+buff_start<refgenome.length() && refgenome[current_index+buff_start-1]<'a'){
                
                fprintf(sampling_file,"%d,",current_index+buff_start);
            }
            else {
                
                fprintf(sampling_file,"-%d,",current_index+buff_start);
            }
            
            
            for (int i=0;i<8;++i){
                
                fprintf(sampling_file,"%d,",allele_counts[i]);
            }
            
            for (int i=0;i<8;++i){
                
                fprintf(sampling_file,"%d,",sumquals[i]);
            }
            
            fprintf(sampling_file,"\n");
            
        }
        
    }
    
    if (writing_mode>1) return 1;
    
    bool ifvar=Ifsnp(total, themax, thesecond);
    
    if (writing_mode){
        
        
        ifvar=(ifvar || (themax>5 && maxindex!=5 && refgenome.length()>current_index+buff_start && Nucleotides(maxindex)!=toupper(refgenome[current_index+buff_start-1]) && toupper(refgenome[current_index+buff_start-1])!='N'));
        
    }
    
    return ifvar;
}

int Snp_Records::Load_Line(int current_index){
    
    Base * current_base=Base_Buffer[current_index];
    
    int *current_snp; char *current_quals;
    char insert[20]="";
    int current_num;
    
    int new_Linesize=(int)(12*accumulate(current_base->lens, current_base->lens+8, 100)) + 10*(current_base->lens[7])+22*(current_base->lens[4]);
    
    if (new_Linesize>Linesize){
        
        writeLine=(char*)realloc(writeLine,new_Linesize*sizeof(char));
        
        Linesize=new_Linesize;
    }
    
    int current_Linesize=0;
    
    if (current_index+buff_start<refgenome.length() && refgenome[current_index+buff_start-1]>='a'){
        current_Linesize+=snprintf(writeLine,66, "-%d\t",current_index+buff_start);
    }
    else{
        current_Linesize+=snprintf(writeLine,66, "%d\t",current_index+buff_start);
    }
    
    short *inserts_prt=current_base->inserts;
    
    intpair* deletions_prt=current_base->deletions;
    
    for (int j=0;j<8; ++j){
        
        current_num=current_base->lens[j];
        
        current_quals=current_base->quals[j];
        
        if (current_num)
        {
            current_snp=current_base->snps[j];
            
            for (int k=0; k<current_num;++k)
            {
                if (j==4) {
                    
                    current_Linesize+=snprintf(writeLine+current_Linesize,66, "%c%d:%d-%d,",current_quals[k],current_snp[k],current_index+buff_start-deletions_prt[k].first,deletions_prt[k].second-current_index-buff_start);
                }
                
                else if (j==7) {
                    
                    Insert_Unhash(inserts_prt[k], insert);
                    
                    current_Linesize+=snprintf(writeLine+current_Linesize,66, "%c%d:%s,",current_quals[k],current_snp[k],insert);
                }
                
                else{
                    
                    current_Linesize+=snprintf(writeLine+current_Linesize,66, "%c%d,",current_quals[k],current_snp[k]);
                    
                }
                
            }
        }
        
        current_Linesize+=snprintf(writeLine+current_Linesize,66, "\t");
        
    }
    
    return current_Linesize;
    
}




bool Snp_Records::WriteSnp(int write_end) {
        
    write_start=max((int)buff_start,write_start); 
    write_end=write_end-100;
    
    int valide_write_end=min((int)buff_end, write_end);
    
    int write_size=write_end-write_start;
    int valide_write_size=min((int)buff_end, write_end)-write_start;
    
    int release_end=write_end-100;
    int valide_release_size=min((int)buff_end, release_end)-buff_start;
    
    if (write_size<=0) return 1;
    
    int realign_start=max(buff_start,last_realign_end);
        
    if (valide_write_end>realign_start) {
        
        aligner.Buffer_Region(realign_start, valide_write_end);
        
        for (int current_index=realign_start;current_index<min((int)buff_end, release_end);current_index++){
            
            aligner.Realign(current_index);
        
        }
        
        last_realign_end=release_end;
    }
    
    if (valide_write_size>=0) {
        
        Correct_SNPs(valide_write_end);
        
        srand((unsigned)time(NULL));
        
        FILE *f = fopen(outfile, "a");
        
        FILE *sampling_file = fopen(samplefile.c_str(), "a");
        
        
        if(f==NULL)
        {
            
            fclose(f);
            
            std::cerr << "ERROR: Could not open file " << outfile << " for writing.\n" <<std::endl;
            
            std::_Exit(EXIT_FAILURE);
            
            return false;
        }
        
        //write snp data in buff memory within write_size region;
        for (int current_index=write_start-buff_start; current_index<valide_write_end-buff_start; ++current_index){
            
            if (!(unsigned char)Base_IfInit[current_index]){continue;}
            
            if (!Alelle_Stat(current_index, sampling_file, (rand() % 10))) continue;
            
            Load_Line(current_index);
            
            fprintf(f,"%s\n", writeLine);
        }
        
        fclose(f);
        
        fclose(sampling_file);
                
    }
    
    if (valide_release_size>=0) Base_Buffer_ctrl.move(valide_release_size);
    
    buff_start=release_end;
    
    write_start=write_end;
    
    return 1;
    
}


void Snp_Records::set_Region(int phase_start){
    
    buff_start=max(1,phase_start-100);
    write_start=max(1,phase_start-1);
    buff_end=write_start+1;
    
}




void Snp_Records::Correct_SNPs(int write_end){
    
    if (!Corr_history.size()) return;
    
    unordered_map<int,unordered_map<int,vector<int>>> Move_records;
    
    vector<int> corrected;
    //std::unordered_map<int, std::unordered_map<int, std::intpair>>
    for (auto record:Corr_history){
        
        int move_from=record.first;

        if (move_from>=write_end) continue;
        
        auto move_record=record.second;
        
        for (auto mole_record: move_record) {
            
            int imole=mole_record.first;
            int iallele=mole_record.second.first;
            int move_to=mole_record.second.second;
            
            Move_records[move_to][iallele].push_back(imole);
        }
        
        for (auto record1:Move_records){
            
            for (auto record2: record1.second){
                
                Base_Buffer_ctrl.transfer(move_from, record1.first, record2.first, record2.second);
            }
        }
        
        corrected.push_back(move_from);
    }
    
    for (int move_from:corrected) {
        
        Corr_history.erase(move_from);
        
    }
    
}



