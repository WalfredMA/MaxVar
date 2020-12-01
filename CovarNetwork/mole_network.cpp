//
//  mole_network.cpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "mole_network.hpp"
#include "snp_network.hpp"

using namespace std;

Mole_network::Mole_network()
{
    
    vector<int>bar_p, bar_n;
    
    
    mole_matrix_p.push_back(bar_p);
    mole_matrix_n.push_back(bar_n);
    
    mole_barcode.push_back(0);
    badindex.push_back(0);
    
    hash_coordi[0][0].push_back(0);
    hash_index[0][0].push_back(0);
    
    mole_counter=1;
     
    badbar=NULL;
    
    
}


void Mole_network::set_hash(vector<int> &vecx_p, vector<int> &vecx_n, vector<bool> &badbar)
{
    
    for (int current_index=0; current_index<40; current_index++)
    {
        
        hash_dic[current_index].clear();
        hash_dic_sign[current_index].clear();
    }
    
    
    int current_num=0;
    for (int current_index=0; current_index<vecx_p.size(); current_index++)
    {
        
        current_num = vecx_p[current_index];
        if (badbar[current_num]){continue;}
        hash_dic[current_num%40].push_back(current_num);
        hash_dic_sign[current_num%40].push_back(1);
        
    }
    
    for (int current_index=0; current_index<vecx_n.size(); current_index++)
    {
        
        current_num = vecx_n[current_index];
        if (badbar[current_num]){continue;}
        hash_dic[current_num%40].push_back(current_num);
        hash_dic_sign[current_num%40].push_back(0);
        
    }
    
    
}

void Mole_network::itsect(vector<int> &vecy_p, vector<int> &vecy_n, vector<bool> &badbar, int *posi_counter, int *nega_counter)
{
    
    vector<int>* hash_cell;
    vector<bool>* hash_cell_sign;
    
    
    int current_num=0;
    for (int current_index=0; current_index<vecy_p.size(); current_index++)
    {
        
        current_num = vecy_p[current_index];
        
        if (badbar[current_num]){continue;}
        
        hash_cell=&hash_dic[current_num%40];
        hash_cell_sign=&hash_dic_sign[current_num%40];
        
        vector<int>::iterator iter = find((*hash_cell).begin(),(*hash_cell).end(), current_num);
        
        if(iter!=(*hash_cell).end())
        {
            
            if ((*hash_cell_sign)[iter-(*hash_cell).begin()])
            {
                
                (*posi_counter)++;
                
            }
            
            else
            {
                
                (*nega_counter)++;
            }
        }
    }
    
    for (int current_index=0; current_index<vecy_n.size(); current_index++)
    {
        
        current_num = vecy_n[current_index];
        
        if (badbar[current_num]){continue;}
        
        hash_cell=&hash_dic[current_num%40];
        hash_cell_sign=&hash_dic_sign[current_num%40];
        
        vector<int>::iterator iter = find((*hash_cell).begin(),(*hash_cell).end(), current_num);
        
        if(iter!=(*hash_cell).end())
        {
            
            if (!(*hash_cell_sign)[iter-(*hash_cell).begin()])
            {
                
                (*posi_counter)++;
                
            }
            
            else
            {
                
                (*nega_counter)++;
            }
        }
    }
}


int Mole_network::findindex(int coordi)
{
    
    int hash0_=coordi%hash0;
    int hash1_=coordi%hash1;
    
    vector<int> hash_cell=hash_coordi[hash0_][hash1_];
    
    
    vector<int>::iterator iter = find(hash_cell.begin(),hash_cell.end(), coordi);
    
    if(iter==hash_cell.end())
    {
        
        hash_coordi[hash0_][hash1_].push_back(coordi);
        hash_index[hash0_][hash1_].push_back(mole_counter);
        mole_barcode.push_back(coordi);
        badindex.push_back(0);
        covar_num.push_back(0);
        covar_dom.push_back(0);
        
        vector<int> newvec0, newvec1;
        mole_matrix_p.push_back(newvec0);
        mole_matrix_n.push_back(newvec1);
        
        return mole_counter++;
    }
    
    return hash_index[hash0_][hash1_][iter-hash_cell.begin()];
    
}



void Mole_network::load_mole(int coordi, vector<int> &bar_p, vector<int> &bar_n){
    
    int size_p=(int)bar_p.size();
    int size_n=(int)bar_n.size();
    
    for (int i=0;i<size_p;i++)
    {
        
        mole_matrix_p[bar_p[i]].push_back(coordi);
        
    }
    
    for (int i=0;i<size_n;i++)
    {
        
        mole_matrix_n[bar_n[i]].push_back(coordi);
        
    }
    
}

void Mole_network::covar(float covar_cutoff){
    
    int posi_cov=0, nega_cov=0, current_abs=0,current_match=0, total_cov=0, match_cov=0;
    
    for (int i=1; i<mole_counter; i++)
    {
        
        if (!i%1000) cout<<i<<endl;
        
        total_cov=0;
        match_cov=0;
        
        set_hash(mole_matrix_p[i], mole_matrix_n[i], *badbar);
        
        for (int connect_index: mole_mole[i])
        {
            
            posi_cov=0;
            nega_cov=0;
            
            if (connect_index==i) continue;
            
            if (badindex[connect_index]){continue;}
            
            itsect(mole_matrix_p[connect_index], mole_matrix_n[connect_index],*badbar,  &posi_cov, &nega_cov);
            
            current_abs=posi_cov+nega_cov-1;
            
            if (current_abs<=0) continue;
            
            current_match= max(0,abs(posi_cov-nega_cov)-1);
            
            if (current_abs>maxvalue)
            {
                
                current_match=current_match*maxvalue/current_abs;
                current_abs=maxvalue;
                
            }
            
            total_cov+=current_abs;
            match_cov+=current_match;
            
        }
        
        //cout<<"mole:"<<","<<total_cov<<","<<match_cov<<endl;
        if (match_cov<covar_cutoff*total_cov)
        {
            
            covar_num[i]=0;
            covar_dom[i]=0;
            badindex[i]=1;
            
        }
        else
        {
            covar_num[i]=match_cov;
            covar_dom[i]=total_cov;
            badindex[i]=0;
        }
        
    }
    
    cout<<"filtered molecules "<<accumulate(badindex.begin(), badindex.end(), 0.0)<<" out of "<<badindex.size()<<endl;

    
}


