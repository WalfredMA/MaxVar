//
//  realign.cpp
//  testcpp
//
//  Created by Wangfei MA on 8/10/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "realign.hpp"
#include <vector>
#include <set>
#include <map>
#include <algorithm>

using namespace std;



//find the minor allele index, only allows biallellic
int FindSNP(int* alleles, int *max_allele, int *second_allele, int cutoff){
    
    int alleles_cut=0;
    
    for (int i=0;i<8;++i){
        
        alleles_cut+=alleles[i];
    }
    alleles_cut=max(3,alleles_cut/cutoff);
    
    
    *max_allele=-1;
    *second_allele=-1;
    int max_count=0,snp_count=0;
    
    for (int i=0;i<8;++i){
        
        if (alleles[i]>alleles_cut){
            ++snp_count;
            
            
            if (alleles[i]>max_count){
                
                *second_allele=*max_allele;
                *max_allele=i;
                max_count=alleles[i];
            }
            else{
                *second_allele=i;
            }
        }
    }
    
    return snp_count;
}

//find intersection of two sets
vector<int> set_match(vector<int> set1, set<int> set2) {
    
    vector<int> intersect;
    
    set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(),
                     std::inserter(intersect,intersect.begin()));
    
    return intersect;
}

//find compliment of two sets
vector<int> set_unmatch(vector<int> set1, set<int> set2) {
    
    vector<int> intersect;
    
    set_difference(set1.begin(),set1.end(),set2.begin(),set2.end(),
                     std::inserter(intersect,intersect.begin()));
    
    return intersect;
}

//continues alignments in order
bool consecutive_comp(vector<int> &x, int x_start, vector<int> &y, int y_start, int size){
    
    
    for (int index=0; index<size; ++index){
        
        if (x[x_start+index]!=y[y_start+index]) return 0;
    }
    
    return 1;
}

//unzip the hashed insert, output to des according to oritation (orit)
int unzip_insert(unsigned short inserts_hash, vector<int>& des, int orit){
    
    int insert[10];
    int insert_len=0;
    
    if (inserts_hash>=10){
        
        inserts_hash-=10;
        while(inserts_hash){
            
            insert[insert_len++]=inserts_hash%4;
            inserts_hash>>=2;
        }
        --insert_len;
    }
    
    else{
        
        return 0;
        
    }
    
    if ( !orit) {
        
        for (int index=0; index<insert_len; ++index){
            
            des.push_back(insert[index]);
        }
        
        return insert_len;
    }
    
    else if (orit<0){
 
        insert_len=min(5,insert_len);
		int push_size=insert_len;
		push_size=(insert_len+1)/2;
       
        //negative oritation to storage, only right side two bases are storaged
        for (int index=insert_len-1; index>=push_size; --index){
            
            des.push_back(insert[index]);
        }
        
        return 5-push_size;
    }
    
    else {
 
        insert_len=min(5,insert_len);
        int push_size=insert_len;
        push_size=(insert_len+1)/2;

        //positive oritation to storage, only left side two molecules are storaged
        for (int index=0; index<push_size; ++index){
            
            des.push_back(insert[index]);
        }
        
        return push_size;
    }
}

//pairwise alignment of 4 bases
int realign_pairwise(vector<int>& query, vector<int>& ref){
    
    int align_size=(int)min(ref.size(),query.size());
    
    int align_size_l=(int)min(ref.size(),query.size()-1);
    
    int align_size_r=(int)min(ref.size()-1,query.size());
    
	align_size=min(4,align_size);

	align_size_l=min(5,align_size);

	align_size_r=min(5,align_size);

    if (align_size<3) return 4-align_size;
    
    int error=0;
    
    for (int index=0; index<align_size; ++index){
        
        if (query[index]==ref[index]) continue;
        
        else {
            
            if (++error>1) return 2;
            
            if (index==3) break;
            
            if (consecutive_comp(query, index+1, ref,index+1, align_size-1-index)) break; //substitution
            
            if (consecutive_comp(query, index+2, ref,index+1, align_size_l-1-index)) break; //insertion
            
            if (consecutive_comp(query, index, ref,index+1, align_size_r-index)) break; //deletion
            
            return 2;
            
        }
    }
    
    return error; //return error count
}


void Realigner::Realign_Allele(){
    
    num_ref_Molecules=current_snp->lens[target_allele];

    int temp_abs[num_ref_Molecules];
    for (int i=0;i<num_ref_Molecules;++i){
        
        temp_abs[i]=abs(current_snp->snps[target_allele][i]);
    }

    ref_Molecules=set<int> (temp_abs, temp_abs+num_ref_Molecules);
    num_ref_Molecules=(int)ref_Molecules.size();

    vector<int> othersnps;
    vector<vector<int>> othersnps_moles;
    if (!Find_SameSNP(othersnps,othersnps_moles)) return;

    if (!Generate_newRef())return ;

    Realign(othersnps,othersnps_moles);
}


//find all similar snp and realign them to current snp
void Realigner::Realign(int newtarget){
    
    target_pos=newtarget-buffer_start;
    
    if (!(unsigned char)buffered_region_ifinit[target_pos]) return;
    
    
    current_snp=buffered_region[target_pos];
    
    ref.up.clear();
    ref.back.clear();
    ref.var.clear();
    
    if (refgenome[newtarget]>='a') return;
    
    if(FindSNP(current_snp->lens, &major_allele, &minor_allele,5)!=2) return;
    
    if (major_allele<0 || minor_allele<0 || major_allele==5 || major_allele==6 || minor_allele==5 || minor_allele==6) return;
    
    //realign on major allele
    target_allele=major_allele;
    
    Realign_Allele();
    
    //realign on minor allele
    target_allele=minor_allele;
    
    Realign_Allele();
    
    
}

//find the same type snp nearby.
bool Realigner::Find_SameSNP(vector<int>& othersnps, vector<vector<int>>& othersnps_moles){
    
    Base* flanking_base;
    int allele_count,flanking_major,flanking_minor;
    for (int index=max(0,target_pos-50);index<min(search_end,target_pos+50);++index){
        
        if (index==target_pos || !(unsigned char)buffered_region_ifinit[index] || refgenome[index+buffer_start]>='a') continue;
        
        flanking_base=buffered_region[index];
                
        if(FindSNP(flanking_base->lens, &flanking_major, &flanking_minor, 10)<2 || (target_allele!=flanking_major && target_allele!=flanking_minor)) continue;
        
        allele_count=flanking_base->lens[target_allele];
        
        //must be more than 3 Molecules or 20% of coverage to be consider realignment
        //if (allele_count<=3 || /* (*max_element(flanking_base->lens,flanking_base->lens+8))==target_allele ||  */ allele_count<=(int)accumulate(flanking_base->lens, flanking_base->lens+8, 0)/5) continue;
        
        vector<int> current_moles (flanking_base->snps[target_allele],flanking_base->snps[target_allele]+allele_count);
        
        if (target_allele==4){
            
            intpair* alldeletions=flanking_base->deletions;
            
            current_moles.clear();
            
            for (int j=0; j<allele_count;++j){
                
                if (alldeletions[j].first>target_pos || alldeletions[j].second<target_pos){
                    
                    current_moles.push_back(abs(flanking_base->snps[target_allele][j]));
                }
            }
            
            allele_count=(int)current_moles.size();
            
        }
        
        else{
            for (int j=0; j<allele_count;++j){
                
                current_moles[j]=abs(current_moles[j]);
                
            }
        }
        
        //must not overlap Molecules of current snp
        vector<int> unmatch=set_unmatch (current_moles, ref_Molecules);

        if ((int)(allele_count-unmatch.size())>= 0.2*min(allele_count,num_ref_Molecules)) continue;

 
        othersnps.push_back(index);
        othersnps_moles.push_back(unmatch);
        
    }
    
    return othersnps.size();
}

bool Realigner::Generate_newRef(){
    
    vector<vector<int>> ref_ups, ref_backs, ref_vars;
    Seed ref_mole;
    
    for (int mole: ref_Molecules){
        
        ref_mole.up.clear();
        ref_mole.back.clear();
        ref_mole.var.clear();
        
        //distract data of current mole from buffer
        if (!Distract_Mole(target_pos,abs(mole), target_allele, &ref_mole, 4)) continue;
        
        ref_ups.push_back(ref_mole.up);
        ref_backs.push_back(ref_mole.back);
        
        if (target_allele==7) ref_vars.push_back(ref_mole.var);
    }
    
    //generate reference from all distracted molecules
    if (!Summary_Moles(ref_ups,ref.up)) return 0;
    if (!Summary_Moles(ref_backs,ref.back)) return 0;
    //if (target_allele==7 && !Summary_Moles(ref_vars,ref.var)) return 0;
    return 1;
    
}


//realign candidate molecules to the reference
void Realigner::Realign(vector<int>& othersnps, vector<vector<int>>& othersnps_moles){
    
    int othersnp, total_mole;
    vector<int> othersnp_moles,allmove_moles;
    
    //iter each candidate snp sites
    for (int index=0;index<othersnps.size();++index){
        
        allmove_moles.clear();
       
        othersnp=othersnps[index];
        
        total_mole=buffered_region[othersnp]->lens[target_allele];
        
        othersnp_moles=othersnps_moles[index];

        //distract candidate molecules for realignment
        for (int mole: othersnp_moles) {

            mole=abs(mole);
            
            Seed mole_data;
            
            if (!Distract_Mole(othersnp, mole, target_allele, &mole_data,5) && Realign_newRef(mole_data, abs(othersnp-target_pos))) {
				    
                //if align, correct
                allmove_moles.push_back(mole);
                
                cout<<buffer_start+target_pos<<"<-"<<buffer_start+othersnp<<":"<<target_allele<<":"<<mole<<endl;
                Vprint(mole_data.up)
                Vprint(ref.up)
                
                Vprint(mole_data.back)
                Vprint(ref.back)
            }
        }
        
        //correct all molecules, if most molecules need to be corrected, then clear whole snp
        if (allmove_moles.size()<=min(2,total_mole/10)) continue;
        
        printf("Correct %d molecule's Allele %c on POS:%d to POS: %d\n",(int)allmove_moles.size(),Nucleotides(target_allele), othersnp+buffer_start, target_pos+buffer_start);
        
        Add_History(target_pos+buffer_start, othersnp+buffer_start, target_allele, allmove_moles);
        
    }
    
}


bool Realigner::Realign_newRef(Seed query, int distance){
    
    int error=0;
	int errorcut=1;

    error+=realign_pairwise(query.up, ref.up);
    if (error>errorcut) return 0;
    //error+=realign_pairwise(query.var, ref.var)+max(0,(int)(query.var.size()-ref.var.size()-1));
    //if (error>1) return 0;
    error+=realign_pairwise(query.back, ref.back);
    if (error>errorcut) return 0;
    
    return 1;
}

//set buffer region
void Realigner::Buffer_Region(int start, int end){
    
    buffer_start=start;
    search_end=end-start;
    
    buffered_region=&Base_Buffer[start-buff_start];
    buffered_region_ifinit=&Base_IfInit[start-buff_start];
    
}

//generate reference
bool Realigner::Summary_Moles(vector<vector<int>>& ref_moles,vector<int>& des){
    
    int counts[4][4]={0};
    
    for (vector<int> ref_mole: ref_moles){
        
        for (int i=0;i<min(4,(int)ref_mole.size());++i){
            if (ref_mole[i]>=0) ++counts[i][ref_mole[i]];
        }
    }
    
    int max_allele, max_count, total_count;
    for (int i=0;i<4;++i){
        
        max_allele=-1;
        total_count=0;
        max_count=0;
        for (int j=0;j<4;++j){
            
            total_count+=counts[i][j];
            if (counts[i][j]>3 && counts[i][j]>max_count){
                
                max_count=counts[i][j];
                max_allele=j;
            }
        }
        
        if (max_count<=total_count/2) return 0;
        
        des.push_back(max_allele);
    }
    
    
    
    return 1; 
}

//distract molecule data from buffer
int Realigner::Distract_Mole(int pos,int mole, int allele, Seed *results, int size){
    
    Base* flanking_base;
    
    int allele_count,ifind;
    int *allele_moles, index;
    int iffind=-1;
    
    int left_size=0,right_size=0,allele_size=1;
    
    //if it is insertion, then storage its hash
    if (allele==7){
        
        flanking_base=buffered_region[pos];
        allele_count=flanking_base->lens[allele];
        allele_moles=flanking_base->snps[allele];
        
        ifind=(int)(find(allele_moles,allele_moles+allele_count,mole)-allele_moles);
        if (ifind==allele_count) ifind=(int)(find(allele_moles,allele_moles+allele_count,-mole)-allele_moles);
        
        vector<int> insert;
        if (ifind!=allele_count) {
            allele_size=unzip_insert(flanking_base->inserts[ifind], insert,0);
        }
        
        for (int index=(allele_size+1)/2-1; index>=0; --index){
                
            results->up.push_back(insert[index]);
        }

		left_size+=(allele_size+1)/2;    
    
        for (int index=(allele_size+1)/2; index<allele_size; ++index){
                
            results->back.push_back(insert[index]);
        }

		right_size+=allele_size/2;
        
    }
    
    //front 4 bases of currect snp
    for (index=pos+1;index<search_end;++index){
        
        iffind=Search_Buffer(index, mole,results->back, 1);
            
        if (iffind>=0) right_size+=iffind;
        
        if (iffind<0 || iffind>1 || right_size>=size) break;
        
    }
    //distract back 4 bases of current snp
    for (index=pos-1;index>=0;--index){
        
        iffind=Search_Buffer(index, mole,results->up, -1);
        
        if (iffind>=0) left_size+=iffind;
        if (iffind<0 || iffind>1 || left_size>=size) break;
        
    }
    if ((left_size+right_size)<6) return 0;

    return allele_size;
    
}


void Realigner::Add_History(int target_pos, int othersnp, int iallele, vector<int>& moles){
    
    lock_guard<mutex> His(realign_record_lock);
    
    for (int imole: moles){
        
        Corr_history[othersnp][imole]=make_pair(iallele,target_pos); 
    }
    
}


int Realigner::Search_Buffer(int index, int imole, vector<int>& des, int sign){
    
    if (!(unsigned char)buffered_region_ifinit[index]) {
        
        return -1;
    }

    int *allele_moles;
    Base * flanking_base=buffered_region[index];
    int allele_count, findindex, iffind=-1;
    
    for (int iallele=0; iallele<7; ++iallele) {
        
        allele_count=flanking_base->lens[iallele];
        
        if (allele_count) {
            
            allele_moles=flanking_base->snps[iallele];
            
            findindex=(int)(find(allele_moles,allele_moles+allele_count,imole)-allele_moles);
            if (findindex==allele_count) findindex=(int)(find(allele_moles,allele_moles+allele_count,-imole)-allele_moles);
            
            if (findindex < allele_count){
                
                if (iallele!=4) {
                    des.push_back(iallele);
                    iffind=1;
                }
                else {
                    iffind=0;
                }
                break;
            }
        }
    }
    
    
    if (iffind>=0) return iffind;
    
    allele_count=flanking_base->lens[7];
    
    if (allele_count) {
        
        allele_moles=flanking_base->snps[7];
        
        findindex=(int)(find(allele_moles,allele_moles+allele_count,imole)-allele_moles);
        if (findindex==allele_count) findindex=(int)(find(allele_moles,allele_moles+allele_count,-imole)-allele_moles);
        
        if (findindex < allele_count){
        
            iffind=unzip_insert(flanking_base->inserts[findindex], des,sign);
        }
    }
    
    return iffind;
}

