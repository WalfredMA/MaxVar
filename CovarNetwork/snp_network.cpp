//
//  snp_network.cpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "snp_network.hpp"



using namespace std;
Snp_network::Snp_network()
{
    
    vector<int>bar_p, bar_n;
    
    snp_matrix_p.push_back(bar_p);
    snp_matrix_n.push_back(bar_n);
    
    snp_coordi.push_back(0);
    badindex.push_back(0);
    
    hash_coordi[0][0].push_back(0);
    hash_index[0][0].push_back(0);
    
    snp_counter=1;
    
    badbar=NULL;
    
}

void Snp_network::find_uniq(vector<int>& list)
{
    
    sort(list.begin(), list.end());
    
    list=vector<int>(list.begin(),unique(list.begin(), list.end()));
    
    
}


void Snp_network::set_hash(vector<int> &vecx_p, vector<int> &vecx_n, vector<bool>& badbar){
    
    for (int current_index=0; current_index<40; current_index++){
        
        hash_dic[current_index].clear();
        hash_dic_sign[current_index].clear();
    }
    
    
    int current_num=0;
    for (int current_index=0; current_index<vecx_p.size(); current_index++){
        
        current_num = vecx_p[current_index];
        if (badbar[current_num]){continue;}
        hash_dic[current_num%40].push_back(current_num);
        hash_dic_sign[current_num%40].push_back(1);
        
    }
    
    for (int current_index=0; current_index<vecx_n.size(); current_index++){
        
        current_num = vecx_n[current_index];
        if (badbar[current_num]){continue;}
        hash_dic[current_num%40].push_back(current_num);
        hash_dic_sign[current_num%40].push_back(0);
        
    }
    
    
}

void Snp_network::itsect(vector<int> &vecy_p, vector<int> &vecy_n, vector<bool>& badbar, int *posi_counter, int *nega_counter)
{
    
    vector<int>* hash_cell;
    vector<bool>* hash_cell_sign;
    
    for (int pair_index : vecy_p)
    {
        
        if (badbar[pair_index]){continue;}
        
        hash_cell=&hash_dic[pair_index%40];
        hash_cell_sign=&hash_dic_sign[pair_index%40];
        
        vector<int>::iterator iter = find((*hash_cell).begin(),(*hash_cell).end(), pair_index);
        
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
    
    for (int pair_index : vecy_n)
    {
        
        if (badbar[pair_index]){continue;}
                
        hash_cell=&hash_dic[pair_index%40];
        hash_cell_sign=&hash_dic_sign[pair_index%40];
        
        vector<int>::iterator iter = find((*hash_cell).begin(),(*hash_cell).end(), pair_index);
        
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



int Snp_network::findindex(int coordi)
{
    
    int hash0_=coordi%hash0;
    int hash1_=coordi%hash1;
    
    vector<int> hash_cell=hash_coordi[hash0_][hash1_];
    
    
    vector<int>::iterator iter = find(hash_cell.begin(),hash_cell.end(), coordi);
    
    if(iter==hash_cell.end()) return 0;

    
    return hash_index[hash0_][hash1_][iter-hash_cell.begin()];
    
}




int Snp_network::load_snp(int coordi, vector<int>& bar_p, vector<int>& bar_n)
{
    int hash0_=coordi%hash0;
    int hash1_=coordi%hash1;
    
    snp_matrix_p.push_back(bar_p);
    snp_matrix_n.push_back(bar_n);
    
    snp_coordi.push_back(coordi);
    badindex.push_back(0);
    covar_num.push_back(0);
    covar_dom.push_back(0);
    
    hash_coordi[hash0_][hash1_].push_back(coordi);
    hash_index[hash0_][hash1_].push_back(snp_counter);
        
    return snp_counter++;
}


void Snp_network::covar(float covar_cutoff)
{
    
    int posi_cov=0, nega_cov=0, current_abs=0, current_match=0;
    
    int coordi=0, last_coordi=-100;
    
    int *total_cov=new int, *match_cov=new int, *sv_counter=new int;
        
    vector<tuple<int*, int*, int*>> cov_buffer;
    cov_buffer.push_back(make_tuple(total_cov, match_cov, sv_counter));
    
    for (int i=1; i<snp_counter; i++)
    {
        
        if (!i%1000) cout<<i<<endl;
        
        coordi=snp_coordi[i];
        
        if (coordi > last_coordi  /* + combine_dis */)
        {
            total_cov=new int;
            match_cov=new int;
            sv_counter=new int;
            
            *total_cov=0;
            *match_cov=0;
            *sv_counter=0;
        }
        
        last_coordi=coordi;
        (*sv_counter)++;
        
        set_hash(snp_matrix_p[i], snp_matrix_n[i], *badbar);
        
        for (int connect_index: snp_snp[i])
        {
            if (connect_index==i) continue;
            
            if (badindex[connect_index]){continue;}
            
            posi_cov=0; nega_cov=0;
            itsect(snp_matrix_p[connect_index], snp_matrix_n[connect_index],*badbar, &posi_cov, &nega_cov);
            
            current_abs=posi_cov+nega_cov-1;
            
            if (current_abs<=0) continue;
            
            current_match= max(0,abs(nega_cov-posi_cov)-1);
            
            if (current_abs>maxvalue)
            {
                
                current_match=current_match*maxvalue/current_abs;
                current_abs=maxvalue;
                
            }
            
            *total_cov+=current_abs;
            *match_cov+=current_match;
    
        }
        
        cov_buffer.push_back(make_tuple(total_cov, match_cov,sv_counter));
    }
    
    
    int *lasttotal_cov=NULL, *lastmatch_cov=NULL, *lastsv_counter=NULL;
    for (int i=1; i<snp_counter; i++)
    {
        total_cov=get<0>(cov_buffer[i]);
        match_cov=get<1>(cov_buffer[i]);
        sv_counter=get<2>(cov_buffer[i]);
        
        if (!(*total_cov && *sv_counter<5)||(*match_cov)<covar_cutoff*(*total_cov))
        {
            covar_num[i]=0;
            covar_dom[i]=0;
            badindex[i]=1;
        }
        
        else
        {
            covar_num[i]=*match_cov;
            covar_dom[i]=*total_cov;
            badindex[i]=0;
        }
        
        if (lasttotal_cov!=total_cov)
        {
            delete lasttotal_cov;
            delete lastmatch_cov;
            delete lastsv_counter;
        }
        
        lasttotal_cov=total_cov;
        lastmatch_cov=match_cov;
        lastsv_counter=sv_counter;
        
    }
    delete lasttotal_cov;
    delete lastmatch_cov;
    delete lastsv_counter;
    
    cout<<"filtered snps "<<accumulate(badindex.begin(), badindex.end(), 0.0)<<" out of "<<badindex.size()<<endl;
    
}







