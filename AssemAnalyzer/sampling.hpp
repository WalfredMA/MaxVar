//
//  stats.hpp
//  snp_caller
//
//  Created by Wangfei MA on 6/30/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef stats_hpp
#define stats_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "error_analyzer.hpp"
#include "coverage_analyzer.hpp"
#include "qual_analyzer.hpp"
#include "ref_reader.hpp"

#define cutoff 10000

#define  Vprint( n )  for ( auto i: n ) { cout<<i<<" ";} cout<<endl;
#define  Aprint( n , type , len ) for ( type i: vector<type> (n, n+ len) ) { cout<<i<<" ";} cout<<endl;


using namespace std;

class sampling
{
    
public:
    
    sampling();

    bool statistics();
    
    void output(string outputfile);
    
    void loadref(const char *reffile, const char *chr);
    
    size_t readsamples(string samplefile);
    
    long long totalnum;
    
    std::map<char, int> qual_cuts, qual_cuts_low;
    
    std::map<char, cover_curve> allele_coverages;
    
    std::map<std::string, error> type_errors;
    
private:
    
    int Read_alleles(string line, int line_data[], int &allele_sum, int &allele_max, bool ifmask);
    
    int get_consecutive(int current_index, char allele);
    
    vector<long long*> allele_counts, error_counts, coverage_errors, quals_counts, quals_errors;

    long long total_correct_moles[9]={1};
    long long total_correct_alleles[9]={1};
    
    Cover_anal cover_anal;
    Qual_anal qual_anal;
    Error_anal error_anal;
    
    int valide_start, valide_end, valide_size;
    
    std::string refgenome;
        
};




#endif /* stats_hpp */
