//
//  Samloader.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/26/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "Samloader.hpp"

using namespace std;

void Samloader::Locate_region(){
    
    std::string StrLine;
    
    Saminfile.open(inputfile, ios_base::in);
    
    if (!Saminfile) {
        
        std::cerr << "ERROR: Could not open " << inputfile << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    streampos Lastpos;
    while (getline(Saminfile,StrLine))
    {
        if (StrLine.length() ==0 || StrLine[0]=='@' || StrLine[0]=='#'){
            
            Lastpos=Saminfile.tellg();
            continue;
        }
    }
    
    Saminfile.seekg(Lastpos); //go back to the lastline end
}

int Samloader::GetNextAlignment(BamTools::BamAlignment &alignment, Cigar_Reader &cigar, int *barindex){ 
    
    std::string StrLine;
    
    if (!getline(Saminfile,StrLine)) return 0;
    
    int elements[20];
    const char *StrLine_c=StrLine.c_str();
    
    int ielements=0;
    for (int i=0;i<StrLine.length();++i){
        
        if (StrLine_c[i]!='\t') continue;
        
        elements[ielements++]=i+1;
        
        if (ielements==11) break;
        
    }
    
    if (ielements<11){
        
        return 0;
    }
    
    
    int AlignmentFlag=(int) atoi(StrLine_c+elements[0]);
    int Position=(int) atoi(StrLine_c+elements[2]);
    int MapQuality=(int) atoi(StrLine_c+elements[3]);
    string chrID(StrLine.begin()+elements[1], StrLine.begin()+elements[2]-1);
    
    int if_filter=Filtration(chrID, Position, MapQuality, AlignmentFlag);
    
    if (if_filter<=0) return if_filter;
    
    int alignment_size=cigar.Load(StrLine_c+elements[4]);
    
    if (alignment_size>=max_Base_Buffer-1000){
        
        printf("ERROR: %s Molecule size too large for current buffer size %d, please considering adjust buffer size in Config.h and recompile.\n", alignment.Name.c_str(), max_Base_Buffer );
        
        return 0;
        
    }
        
    int Qualities_length=elements[10]-elements[9]-1;
    int QueryBases_length=elements[9]-elements[8]-1;
    
    if (Qualities_length>2) {
        
        alignment.Qualities.assign(StrLine_c+elements[9],Qualities_length);
    }
    
    else if (QueryBases_length>2){
        
        alignment.Qualities.resize(QueryBases_length);
        
        fill_n(alignment.Qualities.begin(),QueryBases_length, '!'+20);
    }
    
    else {
        
        return 0;
    }
    
    
    if (QueryBases_length>2){
        
        alignment.QueryBases.assign(StrLine_c+elements[8],elements[9]-elements[8]-1);
    }
    
    else if (alignment_size>2) {
        
        alignment.QueryBases.resize(alignment_size);
        
        if (!MDtag.GrepMD(alignment.QueryBases, StrLine_c)) return 0;
        
    }
    
    alignment.Length=alignment_size;
    
    alignment.Name.assign(StrLine_c,elements[0]-1);
    
    alignment.Position=Position;
    
    long long barcode=Grep_Barcode(StrLine.c_str());
    
    if (!barcode && !read_nobarcode) return 0;
    
    
    if (barcode) {
        
        *barindex=Barcodes_Records.push_back(barcode);
    }
    
    else{
        
        *barindex=Barcodes_Records.push_back(alignment.Name.c_str());
        
    }
    
    if (Position+alignment_size-buff_start>max_Base_Buffer) return 2;
    
    return 1;
    
}

long long Samloader::Grep_Barcode(const char * Str){
    
    int Strlen=(int) strlen(Str);
    
    if (!Strlen) return 0;
    
    const char *bar_search=strstr(Str, "MI:i:");
    
    if (bar_search == NULL ) return 0;
    
    int bar_end;
    int bar_start=(int) (bar_search-Str)+5;
    unsigned long long barcode=0;
    
    for (bar_end=bar_start+1; bar_end< Strlen; ++bar_end){
        
        if (!(Str[bar_end]>='0'&&Str[bar_end]<='9')){
            break;
        }
        
        barcode*=10;
        barcode+=Str[bar_end]-'0';
    }
    
    return barcode;
}

void Samloader::Start_Read(){
    
    vector<thread> threads;
    
    vector<Processor*> Processors;
    
    void (*pFun)(Processor*, Samloader*, int) = &Processor::Initiate;
    
    for(int i = 0; i < nthreads; ++i)
    {
        
        Processors.push_back(new Processor);
    }
    
    for(int i=0; i< nthreads; ++i)
        
    {
        
        threads.push_back(thread(pFun,Processors[i], this, i));
    }
    
    for(int i=0; i< nthreads; ++i)
        
    {
        threads[i].join();
    }
    
}


void Samloader::Close(){
    
    if (Saminfile.is_open()) Saminfile.close();
    
}
