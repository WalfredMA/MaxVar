//
//  Bamloader.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/26/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "Bamloader.hpp"

using namespace BamTools;
using namespace std;


void Bamloader::Locate_region(){
    
    if (!Baminfile.Open(inputfile))
    {
        std::cerr << "ERROR: Could not open " << inputfile << " for reading.\n" <<endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    // Read bam.bai index.
    if (!Baminfile.OpenIndex(baifilename)) 
    {
        std::cerr << "ERROR: Could not read BAI index file, Please index " << baifilename << "\n";
        return ;
    }
    
    int chrID=0;
    if (phase_chrID.length()){
        // Translate from reference name to rID.
        chrID = Baminfile.GetReferenceID(phase_chrID);
        if (chrID<0)
        {
            std::cerr << "ERROR: Reference sequence named " << phase_chrID << " not known.\n";
            return ;
        }
        
        phase_chrID=to_string(chrID);
        
    }
    else{
        phase_chrID="";
    }
    
    // Jump the BGZF stream to this position.
    
    if (!Baminfile.Jump(chrID, max(1,phase_start)))
    {
        std::cerr << "ERROR: Could not find region " << phase_chrID <<":"<< phase_start << "\n";
        return ;
    }
}

int Bamloader::GetNextAlignment(BamTools::BamAlignment &alignment, Cigar_Reader &cigar, int *barindex){  
    
    if (!Baminfile.GetNextAlignment(alignment)) return -1;
    
    int iffilter=Filtration(to_string(alignment.RefID), alignment.Position, alignment.MapQuality, alignment.AlignmentFlag);
    
    if (iffilter<=0) return iffilter;
    
    alignment.Length=cigar.Load(&alignment.CigarData);
    
    if (alignment.Length>=max_Base_Buffer-1000){
        
        printf("ERROR: %s Molecule size too large for current buffer size %d, please considering adjust buffer size in Config.h and recompile.\n", alignment.Name.c_str(), max_Base_Buffer );
        
        return 0;
        
    }
    
    long long Barcode;
    bool iftag=alignment.GetTag("MI", Barcode);
    
    if (!iftag && !read_nobarcode) return 0;
    
    if (iftag) {
        
        *barindex=Barcodes_Records.push_back(Barcode);
    }
    
    else{
        
        *barindex=Barcodes_Records.push_back(alignment.Name.c_str());
        
    }
    
    ++alignment.Position;
    
    
    if (alignment.Length+alignment.Position-buff_start>max_Base_Buffer) return 2;
    
    return 1;
}


void Bamloader::Start_Read(){
    
    vector<thread> threads;
    
    vector<Processor*> Processors;
    
    void (*pFun)(Processor*, Bamloader*, int) = &Processor::Initiate;
    
    for(int i = 0; i < nthreads; ++i)
    {
        
        Processors.push_back(new Processor);
    }
    
    for(int i=0; i< nthreads; ++i)
        
    {
        
        threads.push_back(thread(pFun,Processors[i],this,i));
    }
    
    for(int i=0; i< nthreads; ++i)
        
    {
        threads[i].join();
    }
}


void Bamloader::Close(){
    
    if (Baminfile.IsOpen()) Baminfile.Close();
    
}
