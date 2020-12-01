//
//  call_snp.cpp
//  snp_caller
//
//  Created by Wangfei MA on 7/1/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//


#include "call_snp.hpp"

#define  Vprint( n )  for ( auto i: n ) { cout<<i<<" ";} cout<<endl;
#define  Aprint( n, len )  for (int index=0; index<len; ++index) { cout<<n[index]<<" ";} cout<<endl;
#define  NucleoIndex( base , index) for ( index=0;index<8;++index) {if ("ATCGDNMI"[index]==base) break;}
#define  NucleoBase(theindex) "ATCGDNMI"[theindex]

extern float mean_coverage;

extern int base_qualcuts[8];

extern int base_hqualcuts[8];

extern int confidence;

extern int icycle;

static long long variants_scores[1000]={0};


int find_compli_base(int base)
{
    
    if (base>3) return base;
    
    return 2*(base/2)+(1-base%2);
    
}

bool set_match(vector<int> set1, vector<int> set2)
{
    
    vector<int> intersect;
    
    sort(set2.begin(),set2.end());
    
    set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(),
                     std::inserter(intersect,intersect.begin()));
    
    int size1=(int)set1.size();
    int size2=(int)set2.size();
    
    int size_min=(size1<size2)?size1:size2;
    
    if (intersect.size()>0.8*size_min) return 1;
    
    return 0;
}

bool goodmap(vector<char>& nums)
{
    
    if (!nums.size()) return 0.0;
    
    int map_scores[21]={0};
    
    for (char n:nums)
    {
        map_scores[toupper(n)-'A']++;
    }
    
    int mod_value=(int)(max_element(map_scores, map_scores+20)-map_scores);
    mod_value=1;
    
    int accep_num=(int)accumulate(map_scores, map_scores+mod_value+2, 0);
    
    return (accep_num>nums.size()/2);
    
}


snp_caller::snp_caller(prob_calculator &cal): calculator(cal)
{
    
    valide_start=3*mean_coverage/4;
    
    valide_end=8*mean_coverage/4;
 
    valide_size=valide_end-valide_start;
    
    ifconnect=0;
    
}


float snp_caller::get_consecutive(int current_index, char allele)
{
    
    if (current_index<10 || current_index>refgenome.length()-10) return 0;
    
    int lconsecutive=-1,rconsecutive=-1;
    
    char check;
    for (int coordi=2;coordi<4+21/confidence;++coordi)
    {
        if (toupper(refgenome[current_index-coordi]) != allele)
        {
            if (--lconsecutive<=0)
            {
                break;
            }
            
        }
        
        else
        {
            lconsecutive++;
        }
    }
    
    for (int coordi=0;coordi<2+21/confidence;++coordi)
    {
        check=toupper(refgenome[current_index-coordi]);
        if (toupper(refgenome[current_index+coordi]) != allele)
        {
            if (--rconsecutive<=0)
            {
                break;
            }
        }
        else
        {
            rconsecutive++;
        }
    }
    
    return (float)min(10,max(0,lconsecutive+rconsecutive));
    
}

int snp_caller::motif_index(int current_coordi, int& f_consecutive_index, int& r_consecutive_index)
{
    if (!current_coordi || current_coordi>=refgenome.size()) return 0;
    
    int f_last_base=0, f_next_base=0;
    int r_last_base=0, r_next_base=0;
    
    NucleoIndex(toupper(refgenome[current_coordi-2]),f_last_base);
    NucleoIndex(toupper(refgenome[current_coordi]),f_next_base);
    
    r_last_base=find_compli_base(f_last_base);
    r_next_base=find_compli_base(f_next_base);
    
    f_consecutive_index=f_last_base*4+f_next_base;
    r_consecutive_index=r_next_base*4+r_last_base;
    
    return 0;
}

int snp_caller::MapFilter(Allele_infor& allele)
{
    return goodmap(allele.f_map_scores)+goodmap(allele.r_map_scores);
}

float snp_caller::Strd_evaluate(Allele_infor& current_allele)
{
    
    int f_count = (int) current_allele.f_moles.size();
    int r_count = (int) current_allele.r_moles.size();
    
    float score = abs(log2f(f_count+1)-log2f(r_count+1)-(log2f(current_base.total_f+1)-log2f(current_base.total_r+1)));
    
    return fmax(5-5*score,-5);
}

int snp_caller::SV_Call(Base_infor& current_base)
{
    
    int max_count = current_base.total_raw;
    int out=0;
    
    if (current_base.sv_del.size()>max_count/4)
    {
        
        int f_sv=current_base.sv_del_f, r_sv=current_base.sv_del_r;
        
        
        if (min(f_sv,r_sv)>min(2,(int)current_base.sv_del.size()/4) || Haplo_evaluate(8)>0)
        {
            current_base.alleles_found.insert(4);
            out+=4;
        }
    }
    
    if (current_base.sv_ins.size()>max_count/4)
    {
        int f_sv=current_base.sv_ins_f, r_sv=current_base.sv_ins_r;
        
        if (min(f_sv,r_sv)>min(2,(int)current_base.sv_ins.size()/4) || Haplo_evaluate(9)>0)
        {
            current_base.alleles_found.insert(7);
            out+=7;
        }
    }
    
    return (out);
    
}


float snp_caller::Base_evaluate(Base_infor& current_base)
{
    
    if (current_base.coordinate>0)
    {
        base_score=0.0;
    }
    else
    {
        base_score=-5.0;
    }
    
    int indel_count = max(current_base.allele_counts[4],current_base.allele_counts[7]);
    
    float indel_score = (log2f(current_base.total_raw-indel_count+1) - log2f(indel_count+1) -2.0);
    
    indel_score = (indel_score<0.0)?indel_score:0.0; indel_score = (indel_score>=-5.0)?indel_score:-5.0;

    base_score += indel_score;
    
    return (current_base.max_count<=2*mean_coverage);
}

int snp_caller::Base_phasegroup(Base_infor& current_base)
{
    memset(allele_phases,0,10*sizeof(int));
    memset(posi_errors,0,10*sizeof(float));
    memset(nega_errors,0,10*sizeof(float));
    
    vector<int> sort_index{current_base.max_allele, current_base.second_allele};
    
    vector<pair<vector<int>&,vector<int>&>> all_moles;
    for (int i=0; i<8; ++i)
    {
        pair<vector<int>&,vector<int>&> each_alelle(current_base.Allele[i].f_moles,current_base.Allele[i].r_moles);
        
        all_moles.push_back(each_alelle);
    }
    
    vector<int> sv_empty1,sv_empty2;
    pair<vector<int>&,vector<int>&> del(current_base.sv_del,sv_empty1);
    pair<vector<int>&,vector<int>&> ins(current_base.sv_ins,sv_empty2);
    all_moles.push_back(del);
    all_moles.push_back(ins);
    
    max_label = calculator.BaseHaplo(all_moles, posi_errors, nega_errors);
    
    return max_label;
}

float snp_caller::Haplo_evaluate(int allele_index)
{
    
    if (current_base.allele_counts[allele_index]<2) return 0.0;
    
    
    return calculator.AlleleHaplo(allele_index,current_base.allele_counts ,posi_errors, nega_errors);
}

float snp_caller::Freq_evaluate(Allele_infor& current_allele, int allele_index, int f_consecutive, int r_consecutive, int &f_cutoff, int &r_cutoff, int &f_hcutoff, int &r_hcutoff)
{
    
    int f_motif = f_consecutive+allele_index;
    int r_motif = r_consecutive+find_compli_base(allele_index);
    
    int f_count = (int) current_allele.f_moles.size();
    int r_count = (int) current_allele.r_moles.size();
    
    float score = (calculator.systemfilter(f_motif, r_motif, f_count, r_count, f_cutoff, r_cutoff,f_hcutoff, r_hcutoff));
    
    return score;
    
}

float snp_caller::Qual_evaluate(Allele_infor& current_allele, int allele_index, int f_cutoff, int r_cutoff, int f_high_cut, int r_high_cut)
{
    
    int f_index=allele_index;
    int r_index=find_compli_base(allele_index);
    
    int f_qualcut=base_qualcuts[f_index], r_qualcut=base_qualcuts[r_index];
    int f_hqualcut=base_hqualcuts[f_index], r_hqualcut=base_hqualcuts[r_index];
    
    vector<int> f_moles=current_allele.f_moles,r_moles=current_allele.r_moles;
    vector<int> f_base_scores=current_allele.f_base_scores,r_base_scores=current_allele.r_base_scores;
    vector<char> f_map_scores=current_allele.f_map_scores,r_map_scores=current_allele.r_map_scores;
    
    int f_high_count=0, r_high_count=0;
    
    int f_count1=0,f_count2=0;
    for (int i=0; i<f_base_scores.size();++i)
    {
        if (f_base_scores[i]>f_qualcut)
        {
            f_count1++;
        }
        if (toupper(f_map_scores[i])-'A'<3)
        {
            f_count2++;
        }
        
        if (f_base_scores[i]>f_hqualcut && toupper(f_map_scores[i])-'A'==0 && f_moles[i]>0)
        {
            f_high_count++;
        }
    }
    
    int r_count1=0, r_count2=0;
    for (int i=0; i<r_base_scores.size();++i)
    {
        if (r_base_scores[i]>r_qualcut)
        {
            r_count1++;
        }
        
        if (toupper(r_map_scores[i])-'A'<3)
        {
            r_count2++;
        }
        
        if (r_base_scores[i]>r_hqualcut && toupper(r_map_scores[i])-'A'==0 && r_moles[i]>0)
        {
            r_high_count++;
        }
    }
    
    float score=0.0;
    
    score+= fmax(0,5*max(log2f(f_high_count+1)-log2f(f_high_cut+1),log2f(r_high_count+1)-log2f(r_high_cut+1)));
    
    if (f_count1< f_cutoff && r_count1<r_cutoff && !calculator.BaseQualfilter(f_index,r_index,current_allele.f_base_score , current_allele.r_base_score))
    {
        score-=5.0;
    }
    
    if (f_count2< f_cutoff && r_count2<r_cutoff && !MapFilter(current_allele))
    {
        score-=10.0;
    }
    
    return  score;
    
}



float snp_caller::Mask_evaluate(Allele_infor& current_allele, int current_coordi, int allele_index)
{
    
    float mask_score=0.0;
    
    if (allele_index<4)
    {
        
        if (current_coordi>0 || toupper(refgenome[-current_coordi-1])==NucleoBase(allele_index)) return 0.0;

        mask_score=fmax(get_consecutive(abs(current_coordi), NucleoBase(3-allele_index))/4,get_consecutive(abs(current_coordi), NucleoBase(allele_index)));
    }
    
    else
    {
        for (int i=0;i<4;++i)
        {
            mask_score=fmax(5*get_consecutive(abs(current_coordi), NucleoBase(i)),mask_score);
        }
    }
    
    return -mask_score;
    
}

float snp_caller::Prop_evaluate(int allele_index)
{
    int allele_f[10],allele_r[10];
    for (int i=0;i<8;++i)
    {
        allele_f[i]=(int)current_base.Allele[i].f_moles.size();
        allele_r[i]=(int)current_base.Allele[i].r_moles.size();
    }
    
    int max_count_f=*max_element(allele_f, allele_f+4);
    int max_count_r=*max_element(allele_r, allele_r+4);
    
    int allele_count_f=(int)current_base.Allele[allele_index].f_moles.size();
    int allele_count_r=(int)current_base.Allele[allele_index].r_moles.size();
    
    float prop_score_f = log2f(allele_count_f+1)-log2f(max_count_f+1);
    float prop_score_r = log2f(allele_count_r+1)-log2f(max_count_r+1);
    
    return fmin(10,5*(prop_score_f+prop_score_r)+10);
}

int snp_caller::Indel_Call(Base_infor& current_base, int cutoff=12)
{
    
    vector<int> index_found;
    
    int f_consecutive=0, r_consecutive=0;
    int f_consecutive_0=0, r_consecutive_0=0;
    
    int max_index=current_base.max_allele;
    int second_index=current_base.second_allele;
    
    int coordinate = current_base.coordinate;
    
    set<int>& snp_found = current_base.alleles_found;
    set<int>& snp_lowconfi = current_base.alleles_lowconfi;
    
    motif_index(abs(coordinate), f_consecutive, r_consecutive);
    
    vector<int> indels{7,4};
    for (int i: indels)
    {
        
        if (i!=max_index)
        {
            f_consecutive_0=f_consecutive*64+max_index*8;
            r_consecutive_0=r_consecutive*64+find_compli_base(max_index)*8;
        }
        else
        {
            f_consecutive_0=f_consecutive*64+second_index*8;
            r_consecutive_0=r_consecutive*64+find_compli_base(second_index)*8;
        }
        
        pair<vector<int>&,vector<int>&> molecules(current_base.Allele[i].f_moles,current_base.Allele[i].r_moles);
        
        int ifoutput = Allele_evaluate(current_base.Allele[i], current_base.coordinate, i, f_consecutive_0, r_consecutive_0, allele_scores[i]);
                
        if (ifoutput>0 && allele_scores[i]>indel_cutoff)
        {
            snp_found.insert(i);
        }
        else if (icycle && ifoutput)
        {
            snp_lowconfi.insert(i);
        }
        
    }
    
    return (int)snp_found.size()+(int)snp_lowconfi.size();
    
}

int snp_caller::Allele_evaluate(Allele_infor& current_allele, int coordinate, int allele_index, int f_consecutive, int r_consecutive, float &essential_score)
{
    
    float additional_score = 0.0, phase_score=0.0;
    int f_cutoff, r_cutoff,f_hcutoff, r_hcutoff;
    
    additional_score+=Freq_evaluate(current_allele, allele_index,f_consecutive, r_consecutive, f_cutoff, r_cutoff,f_hcutoff, r_hcutoff);
    
    if (icycle) phase_score=Haplo_evaluate(allele_index);
    
    if (phase_score!=0)
    {
        essential_score+=phase_score;
        
        if (allele_index<4 && phase_score+essential_score>=0 && max(posi_errors[allele_index],nega_errors[allele_index])>5) return 2;
        
        if (allele_index>4 && phase_cutoff<=0) return 0;
    }
    
    cout<<endl<<current_base.coordinate<<",";
    
    essential_score=fmax(essential_score,-10);
    
    cout<<additional_score<<",";
    
    additional_score+=Mask_evaluate(current_allele,coordinate,allele_index);
    
    cout<<additional_score<<",";
    
    additional_score += base_score;
    
    cout<<additional_score<<",";
    
    additional_score+=Qual_evaluate(current_allele, allele_index, f_cutoff, r_cutoff,f_hcutoff, r_hcutoff);
    
    cout<<additional_score<<",";
    
    additional_score+=Strd_evaluate(current_allele);
    
    cout<<additional_score<<",";
    
    additional_score+=Prop_evaluate(allele_index);
    
    cout<<additional_score<<",";
    
    essential_score+=additional_score;
    
    cout<<essential_score<<",";
    
    if (essential_score >= postphase_cutoff)
    {
        return 1;
    }
    else if (essential_score >= baseconfi_cutoff )
    {
        return -1;
    }
    
    return 0;
}


int snp_caller::Allele_Call(Base_infor& current_base, int cutoff=12)
{
    int f_consecutive=0, r_consecutive=0;
    int f_consecutive_0=0, r_consecutive_0=0;
    
    int max_index=current_base.max_allele;
    int second_index=current_base.second_allele;
    
    int coordinate = current_base.coordinate;
    
    set<int>& snp_found = current_base.alleles_found;
    set<int>& snp_lowconfi = current_base.alleles_lowconfi;
    
    motif_index(abs(coordinate), f_consecutive, r_consecutive);
    
    vector<int> sort_index{current_base.max_allele, current_base.second_allele};
    
    for (int i: sort_index)
    {
        if (i!=max_index)
        {
            f_consecutive_0=f_consecutive*64+max_index*8;
            r_consecutive_0=r_consecutive*64+find_compli_base(max_index)*8;
        }
        else
        {
            f_consecutive_0=f_consecutive*64+second_index*8;
            r_consecutive_0=r_consecutive*64+find_compli_base(second_index)*8;
        }
        
        int ifoutput= Allele_evaluate(current_base.Allele[i], current_base.coordinate, i, f_consecutive_0, r_consecutive_0, allele_scores[i]);
        
        if (ifoutput>0)
        {
            snp_found.insert(i);
        }
        else if (icycle && abs(ifoutput))
        {
            snp_lowconfi.insert(i);
        }
        
        if (i==current_base.second_allele && allele_scores[i]>0 && allele_scores[i]<100) variants_scores[(int)(allele_scores[i]*10)]++;
        
    }

    
    return (int)snp_found.size()+(int)snp_lowconfi.size();
}


int snp_caller::filtration(Base_infor& current_base)
{
    
    current_base.alleles_found.clear();
    
    if (icycle) max_label=Base_phasegroup(current_base);
    
    if (SV_Call(current_base))
    {
        for (int allele: current_base.alleles_found)
        {
            current_base.alleles_label.push_back(tolower(NucleoBase(allele)));
        }

        return 0;
        return -(int)2*(bool)current_base.alleles_found.size();
    }
    
    memset(allele_scores,0,(10)*sizeof(float));
    
    if (!Base_evaluate(current_base)) return 0;
    
    Allele_Call(current_base,baseconfi_cutoff);
    
    //if (icycle && current_base.alleles_found.size()<2 && current_base.coordinate>0 && current_base.allele_counts[4]>current_base.second_count) Indel_Call(current_base);
    
    if (calculator.refsnp.size() && (calculator.refsnp.count(abs(current_base.coordinate))))
    {
        
        if (calculator.labels.count(abs(current_base.coordinate))) current_base.alleles_label=calculator.labels[abs(current_base.coordinate)];
        
        int allele_index;
        for (char &allele: current_base.alleles_label)
        {
            allele=toupper(allele);
            
            NucleoIndex(toupper(allele),allele_index)
            
            output_scores.push_back(allele_scores[allele_index]);
        }
        
        return 2;
    }
    else if (calculator.lowsnp.size() && calculator.lowsnp.count(abs(current_base.coordinate)))
    {
        
        if (calculator.labels.count(abs(current_base.coordinate))) current_base.alleles_label=calculator.labels[abs(current_base.coordinate)];
        
        int allele_index;
        for (char &allele: current_base.alleles_label)
        {
            allele=toupper(allele);
            
            NucleoIndex(toupper(allele),allele_index)
            
            output_scores.push_back(allele_scores[allele_index]);
        }
        
        return 2;
    }
    
    
    for (int allele: current_base.alleles_found)
    {
        
        char allele_chr;
        if (allele==4)
        {
            allele_chr = 'E';
        }
        
        else if (allele==7)
        {
            allele_chr = 'N';
        }
        
        else
        {
            allele_chr = toupper(NucleoBase(allele));
        }
        
        current_base.alleles_label.push_back(allele_chr);
        
        output_scores.push_back(allele_scores[allele]);
    }
    
    for (int allele: current_base.alleles_lowconfi)
    {
        
        char allele_chr;
        if (allele==4)
        {
            allele_chr = 'e';
        }
        
        else if (allele==7)
        {
            allele_chr = 'n';
        }
        
        else
        {
            allele_chr = tolower(NucleoBase(allele));
        }
       
        current_base.alleles_label.push_back(allele_chr);
        
        output_scores.push_back(allele_scores[allele]);
    }
    
    return (int) current_base.alleles_label.size();
}


int snp_caller::Var_Call(string &StrLine)
{
    
    reader.ReadLine(StrLine, current_base);
    
    int num_allele = filtration(current_base);
    
    if (abs(num_allele)==2)
    {
        
        StrLine+='\t';
        for (char allele: current_base.alleles_label)
        {
            StrLine.push_back(allele);
            StrLine.push_back('-');
        }
        
        StrLine+='\t';
        for (float score: output_scores)
        {
            StrLine+=to_string(score);
            StrLine.push_back(';');
        }
        
        return 1;
    }
    
    return 0;
    
}


int snp_caller::snp_calling(const char* inputfile, const char* outputfile, const char* reffile, const char* chr)
{
    
    refreader Reader(reffile);
    
    refgenome=Reader.load_chr(chr);
    
    fstream fpread;
    fpread.open(inputfile, ios_base::in);
    
    if (!fpread)
    {
        
        std::cerr << "ERROR: Cannot access file: " << inputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
        
    }
    
    FILE *fwrite=fopen(outputfile, "w");
    
    if (fwrite==NULL)
    {
        
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
        
    }
    
    long long line_index=0;
    
    string inLine;
    pair<int,vector<string>> outLines;
    outLines.first=1;
    
    while (getline(fpread,inLine))
    {
        
        if (Var_Call(inLine))
        {
            
            fprintf(fwrite,"%s\n", inLine.c_str());
        }
        
        line_index++;
        
        if (!((line_index+1)%10000))
        {
            
            cout<<"<<<<<<<<<<< Finished Line: "<<line_index+1<<" >>>>>>>>>>>"<<endl;
            //if (line_index>500000) break;
            
        }
    }
    
    fclose(fwrite);
    
    Aprint(variants_scores, 1000)
    
    printf("Converting finished\n");
    
    return 0;
}

