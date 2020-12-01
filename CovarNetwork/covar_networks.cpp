//
//  networks.cpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "covar_networks.hpp"

int set_match(vector<int>& set1, set<int>& set2)
{
    
    vector<int> intersect;
    
    set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(),
                     std::inserter(intersect,intersect.begin()));
    
    return (int)intersect.size();
}



void covar_networks::sv_innerfilter(vector<tuple<int, vector<int>, vector<int>>>&sv_block)
{
    
    map<int,tuple<int,int,int>> moles;
    

    for (auto sv: sv_block)
    {
        
        for (int mole: get<1>(sv))
        {
            get<0>(moles[mole])++;
            
        }
        
        for (int mole: get<2>(sv))
        {
            get<1>(moles[mole])++;
        }
        
    }
    
    for (auto &mole: moles)
    {
        
        if (get<0>(mole.second) > max(5,3*get<1>(mole.second)))
        {
            get<2>(mole.second)=1;
        }
        else if (get<1>(mole.second) > max(5,3*get<0>(mole.second)))
        {
            get<2>(mole.second)=-1;
        }
        
    }
    
    for (auto &sv: sv_block)
    {
        
        int posi=0, nega=0, zero=0;
        
        for (int mole: get<1>(sv))
        {
            
            if (get<2>(moles[mole])==1)
            {
                posi++;
            }
            else if (get<2>(moles[mole])==-1)
            {
                nega++;
            }
            else
            {
                zero++;
            }
        }
        
        for (int mole: get<2>(sv))
        {
            if (get<2>(moles[mole])==1)
            {
                nega++;
            }
            else if (get<2>(moles[mole])==-1)
            {
                posi++;
            }
            else
            {
                zero++;
            }
        }

        if (!(abs(posi-nega)>5 && abs(posi-nega)>0.6*(posi+nega) && abs(posi-nega)>zero))
        {
            get<1>(sv).clear();
            get<2>(sv).clear();
        }
        
    }
    
}



int covar_networks::read_line(string StrLine, vector<int> &p_col, vector<int> &n_col)
{
    
    size_t str_size=StrLine.length();
    
    char current_chr;
    
    int coordinate=0, curren_num=0, strindex=0;
    
    for (strindex=0;strindex< str_size;strindex++)
    {
        
        current_chr=StrLine[strindex];
        
        if (current_chr>='0'&&current_chr<='9')
        {
            coordinate=coordinate*10+(int)(current_chr-'0');
        }
        
        else if (current_chr=='\t')
        {
            strindex++;
            break;
        }
    }
    
    for (;strindex< str_size;strindex++)
    {
        
        current_chr=StrLine[strindex];
        
        if (current_chr==',')
        {
            p_col.push_back(curren_num);
            
            curren_num=0;
        }
        
        else if (current_chr>='0'&&current_chr<='9')
        {
            
            curren_num=curren_num*10+(int)(current_chr-'0');
            
        }
        
        else if (current_chr=='\t')
        {
            
            if (curren_num) p_col.push_back(curren_num);
            
            curren_num=0;
            
            strindex++;
            
            break;
            
        }
        
    }
    
    for (;strindex< str_size;strindex++)
    {
                
        current_chr=StrLine[strindex];
        
        if (current_chr==',')
        {
            
            n_col.push_back(curren_num);
            
            curren_num=0;
            
        }
        
        else if (current_chr>='0'&&current_chr<='9')
        {
            
            curren_num=curren_num*10+(int)(current_chr-'0');
            
        }
        
        else if (current_chr=='\t'||current_chr=='\n')
        {
            
            if (curren_num) n_col.push_back(curren_num);
            
            curren_num=0;
            
            strindex++;
            
            break;
            
        }
        
    }
    
    if (curren_num) n_col.push_back(curren_num);
    
    for (;strindex< str_size;strindex++)
    {
        
        current_chr=StrLine[strindex];
        
        if (current_chr=='\t'||current_chr=='\n') break;
            
        labels[coordinate].push_back(current_chr);
        
    }
    

    return coordinate;
}



void covar_networks::load_data(string datafile)
{
    
    int line_size;
    
    int line_index=0;
    
    vector<tuple<int, vector<int>, vector<int>>> sv_block;
    
    int sv_start=-100, sv_end=-100, sv_sites=0;
    
    fstream file;
    file.open(datafile);
    
    string line;
    
    int coordinate=0, coordi_index=0;
    while (getline(file,line))
    {
        
        line_size=(int)line.length();
        
        if (!line_size){continue;}
        
        vector<int> p_column, n_column;
        
        coordinate=read_line(line, p_column, n_column);
        
        for (int i=0; i< p_column.size(); i++)
        {
            p_column[i]=mole_network.findindex(p_column[i]);
        }
        
        for (int i=0; i< n_column.size(); i++)
        {
            n_column[i]=mole_network.findindex(n_column[i]);
            
        }
        
        if (coordinate > sv_end + combine_dis)
        {
            
            if (sv_sites>4 && sv_sites>(sv_end-sv_start)/2)
            {
                sv_innerfilter(sv_block);
            }
            
            
            for (auto idata: sv_block)
            {
                
                if (!(++line_index%1000)) cout<<line_index<<endl;

                coordi_index=snp_network.load_snp(get<0>(idata), get<1>(idata), get<2>(idata));
                
                mole_network.load_mole(coordi_index, get<1>(idata), get<2>(idata));
                
            }
            
            sv_start = coordinate;
            sv_sites = 0;
            sv_block.clear();
            
        }
        
        sv_block.push_back(make_tuple(coordinate, p_column, n_column));
        sv_sites++;
        sv_end=coordinate;
        
    }
    
    
    if (sv_sites>10 && sv_sites>(coordinate-sv_start)/2)
    {
        sv_innerfilter(sv_block);
    }
    
    
    for (auto idata: sv_block)
    {
        
        coordi_index=snp_network.load_snp(get<0>(idata), get<1>(idata), get<2>(idata));
        
        mole_network.load_mole(coordi_index, get<1>(idata), get<2>(idata));
        
    }
    
    cout<<"end loading"<<endl;
    
}


void covar_networks::connect_snp()
{
    
    vector<int>* vec_a;
    vector<int>* vec_b;
    
    vector<int> sv_block_index;
    
    size_t size_a;
    int sv_start=-100, sv_end=-100;
    
    vector<int> connected;
    snp_network.snp_snp.push_back(connected);
    
    for (int i=1; i< snp_network.snp_counter; i++)
    {
        
        if (!(i%1000))
        {
            cout<<i<<endl;
        }
        
        int coordi=snp_network.snp_coordi[i];
        
        if (coordi > sv_start + combine_dis) sv_start = coordi;
        
        vector<int> connected;
        
        connected.reserve(sizeof(int)*100);
        
        vec_a=&snp_network.snp_matrix_p[i];
        size_a=(*vec_a).size();
        for (int j:(*vec_a))
        {
            
            vec_b=&mole_network.mole_matrix_p[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
            vec_b=&mole_network.mole_matrix_n[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
        }
        
        vec_a=&snp_network.snp_matrix_n[i];
        for (int j:(*vec_a))
        {
            
            vec_b=&mole_network.mole_matrix_p[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
            vec_b=&mole_network.mole_matrix_n[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
        }
        
        if (connected.size()>1)
        {
            sort( connected.begin(), connected.end() );
            connected.erase( unique( connected.begin(), connected.end() ), connected.end() );
        }
        
        vector<int> front_sites;
        //vector<int> after_sites;
        sv_end=coordi;
                
        int coordi_temp;
        for (int n:connected)
        {
            coordi_temp=snp_network.snp_coordi[n];
            
            if (coordi_temp<sv_start-combine_dis)
            {
                front_sites.push_back(n);
            }
            
            else if (coordi_temp>sv_end+combine_dis)
            {
                front_sites.push_back(n);
            }
            
            else
            {
                sv_end=coordi_temp;
            }
        }
        
        //front_sites.insert(front_sites.end(), after_sites.begin(), after_sites.end());
        
        snp_network.snp_snp.push_back(front_sites);
        
    }
}


void covar_networks::connect_mole()
{
    
    vector<int>* vec_a;
    vector<int>* vec_b;
    
    vector<int> zero(0);
    mole_network.mole_mole.push_back(zero);
    
    for (int i=1; i< mole_network.mole_counter; i++)
    {
        
        if (!(i%1000))
        {
            cout<<i<<endl;
        }
        
        vector<int> connected;
        
        connected.reserve(sizeof(int)*100);
        
        vec_a=&mole_network.mole_matrix_p[i];
        
        for (int j:(*vec_a))
        {
            
            vec_b=&snp_network.snp_matrix_p[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
            vec_b=&snp_network.snp_matrix_n[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
        }
        
        vec_a=&mole_network.mole_matrix_n[i];
        for (int j:(*vec_a))
        {
            
            vec_b=&snp_network.snp_matrix_p[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
            vec_b=&snp_network.snp_matrix_n[j];
            
            connected.insert( connected.end(), (*vec_b).begin(), (*vec_b).end() );
            
        }
        
        sort( connected.begin(), connected.end() );
        connected.erase( unique( connected.begin(), connected.end() ), connected.end() );
        
        mole_network.mole_mole.push_back(connected);
        
    }
    
}


void covar_networks::init(string datafile)
{
    
    cout<<"loading"<<endl;
    
    load_data(datafile);
    
    cout<<"initiating"<<endl;
    
    connect_snp();
    
    connect_mole();
    
    snp_network.badbar=&mole_network.badindex;
    
    mole_network.badbar=&snp_network.badindex;
    
    
}

void covar_networks::cycle(float cutoff)
{
    
    snp_network.covar(cutoff);
    
    mole_network.covar(0.4);
    
}

void covar_networks::pre_filtration()
{
    cout<<"start filter"<<endl;
    
    float global_covar=0.0;
    mole_network.covar(global_covar);
    global_covar=(double)accumulate(mole_network.covar_num.begin(), mole_network.covar_num.end(), 1)/(double)accumulate(mole_network.covar_dom.begin(), mole_network.covar_dom.end(), 1);
        
    
    int theturn=0;
    while (global_covar<0.6)
    {
        
        if (theturn++>5) break;
        
        cout<<"mole"<<endl;
        mole_network.covar(0.4);
        
        global_covar=(double)accumulate(mole_network.covar_num.begin(), mole_network.covar_num.end(), 1)/(double)accumulate(mole_network.covar_dom.begin(), mole_network.covar_dom.end(), 1);
        
        cout<<"snp"<<endl;
        cout<<(double)accumulate(mole_network.covar_num.begin(), mole_network.covar_num.end(),1)<<";"<<(double)accumulate(mole_network.covar_dom.begin(), mole_network.covar_dom.end(), 1)<<endl;
        
        global_covar=(global_covar>0.6)?0.6:global_covar;
        
        snp_network.covar(global_covar);
        
        cout<<theturn<<endl;
        
        
    }
    
}


void covar_networks::run(vector<float> cutoffs)
{
    
    pre_filtration();
    
    for (int icycle=0; icycle<cutoffs.size(); icycle++)
    {
        
        cout<<"filtering cutoff: "<<cutoffs[icycle]<<endl;
        
        cycle(cutoffs[icycle]);
        
        cout<<(double)accumulate(mole_network.covar_num.begin(), mole_network.covar_num.end(),1)<<";"<<(double)accumulate(mole_network.covar_dom.begin(), mole_network.covar_dom.end(), 1)<<endl;
    }
    
    mole_network.covar(cutoffs[cutoffs.size()-1]);
    
}

void covar_networks::output(string outfile)
{
    
    int posicount, negacount;
    
    ofstream file;
    
    string line;
    
    file.open (outfile,std::ofstream::app);
    
    int current_bar;
    
    int output_counter=0;
    for (int i=1; i<snp_network.snp_counter; i++)
    {
        
        if (snp_network.badindex[i] || !snp_network.snp_matrix_p[i].size() || !snp_network.snp_matrix_n[i].size()){continue;}
        
        posicount=0;
        negacount=0;
        
        line.clear();
        
        int coordinate=snp_network.snp_coordi[i];
        
        line+=to_string(coordinate);
        line+="\t";
        
        for (int j=0; j<snp_network.snp_matrix_p[i].size();j++)
        {
            
            current_bar=snp_network.snp_matrix_p[i][j];
            
            if ((*snp_network.badbar)[current_bar]) continue;
            
            posicount++;
            
            line+=to_string(mole_network.mole_barcode[current_bar]);
            line+=",";
        }
        
        line+="\t";
        
        for (int j=0; j<snp_network.snp_matrix_n[i].size();j++)
        {
            
            current_bar=snp_network.snp_matrix_n[i][j];
            
            if ((*snp_network.badbar)[current_bar]) continue;
            
            negacount++;
            
            line+=to_string(mole_network.mole_barcode[current_bar]);
            line+=",";
        }
        
        if ( abs(log2f(posicount+1)-log2f(negacount+1))>2 || min(negacount,posicount)<3) continue;
        
        string label = labels[coordinate];
        
        if (label.length())
        {
            line.push_back('\t');
            line+=label;
        }
        
        file<<line<<endl;
        
        output_counter++;
        
    }
    
    file.close();
    
    cout<<"output snp:"<<output_counter<<endl;
}


float log_add(float x, float y)
{
    
    static int log2addtable[]={100, 100, 99, 99, 98, 98, 97, 97, 96, 96, 95, 95, 94, 94, 93, 93, 92, 92, 91, 91, 90, 90, 89, 89, 88, 88, 88, 87, 87, 86, 86, 85, 85, 84, 84, 84, 83, 83, 82, 82, 81, 81, 81, 80, 80, 79, 79, 78, 78, 78, 77, 77, 76, 76, 76, 75, 75, 74, 74, 73, 73, 73, 72, 72, 72, 71, 71, 70, 70, 70, 69, 69, 68, 68, 68, 67, 67, 67, 66, 66, 65, 65, 65, 64, 64, 64, 63, 63, 63, 62, 62, 62, 61, 61, 61, 60, 60, 60, 59, 59, 58, 58, 58, 58, 57, 57, 57, 56, 56, 56, 55, 55, 55, 54, 54, 54, 53, 53, 53, 52, 52, 52, 52, 51, 51, 51, 50, 50, 50, 49, 49, 49, 49, 48, 48, 48, 47, 47, 47, 47, 46, 46, 46, 46, 45, 45, 45, 44, 44, 44, 44, 43, 43, 43, 43, 42, 42, 42, 42, 41, 41, 41, 41, 40, 40, 40, 40, 39, 39, 39, 39, 38, 38, 38, 38, 38, 37, 37, 37, 37, 36, 36, 36, 36, 36, 35, 35, 35, 35, 34, 34, 34, 34, 34, 33, 33, 33, 33, 33, 32, 32, 32, 32, 32, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29, 29, 29, 29, 29, 29, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 26, 26, 26, 26, 26, 26, 25, 25, 25, 25, 25, 25, 24, 24, 24, 24, 24, 24, 23, 23, 23, 23, 23, 23, 23, 22, 22, 22, 22, 22, 22, 22, 21, 21, 21, 21, 21, 21, 21, 21, 20, 20, 20, 20, 20, 20, 20, 19, 19, 19, 19, 19, 19, 19, 19, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    int x_100 = (int) (x*100);
    int y_100 = (int) (y*100);
    
    int higher=max(x_100, y_100);
    int lower = min(x_100, y_100);
    
    int diff=higher-lower;
    
    if (diff>800)
    {
        return lower/100.0;
        
    }
    
    return 1.0*(lower-log2addtable[diff])/100.0;
    
}


float covar_networks::factorial_log(int end)
{
    
    static float results_storage[1000];
    static bool results_index[1000]={0};
    
    if (end<1000 && results_index[end])
    {
        
        return results_storage[end];
        
    }
    
    float factorial=0.0;
    
    for (int i=2;i<end+1;i++)
    {
        factorial+=log2f(i);
    }
    
    
    if (end<1000)
    {
        
        results_index[end]=1;
        results_storage[end]= factorial;
        
    }
    
    return  factorial;
}


double covar_networks::bionomial(int total, int obs)
{
    
    int notobs=total-obs;
    double pro_log=log2f(0.5);
    double unpro_log=log2f(0.5);
    
    double result= obs* pro_log+notobs* unpro_log +factorial_log(total)-factorial_log(notobs)- factorial_log(obs);
    return -result;
}


float covar_networks::cul_bionom(int total, int obs /*minor obs*/)
{
    
    static float results_storage[51][51];
    static bool results_history[51][51];
    
    if (total>50)
    {
        obs=obs*50/total;
        total=50;
    }
    
    obs=(obs<total-obs)?obs:(total-obs);
    
    if (results_history[total][obs])
    {
        return results_storage[total][obs];
    }
    
    float result=1000.0;
    for (int i=0;i<=obs;i++)
    {
        if (2*i==total) {result=log_add(result, bionomial(total, i));}
        else {result=log_add(result, bionomial(total, i)-1);}
        
    }
    
    result=(result<0.0)?0.0:result;
    
    if (total<51 && obs<51)
    {
        results_storage[total][obs]=result;
        results_history[total][obs]=1;
    }
    return result;
}

