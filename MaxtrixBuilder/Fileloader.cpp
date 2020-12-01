//
//  Fileloader.cpp
//  testcpp
//
//  Created by Wangfei MA on 8/19/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "Fileloader.hpp"

using namespace BamTools;
using namespace std;


int Check_Threads(){
    
    int writeend=last_read_start;
    int new_posi;
    for(int i=2;i<nthreads+2;i++)
    {
        
        new_posi=Threads_Posi[i];
        
        if (writeend<new_posi) {
            
            writeend=new_posi;
        }
    }
    
    last_read_start=writeend;
    
    return writeend;
    
}

int Fileloader::Filtration(string chrID, int Position, int MapQuality, int AlignmentFlag){
    
    if (Position>phase_end || ( phase_chrID.length() && phase_chrID.compare(chrID))) return -1;
    
    if (Position<last_read_start || MapQuality <qualcut || (ifcheckflag && issueflag(AlignmentFlag)) || Position-buff_start<0 ) return 0;
    
    return 1;
    
}

int Fileloader::issueflag(int flag){
    int binary[12]={0};
    int remainder0;
    int left=flag;
    int digit=0;
    
    
    while(left)
    {
        remainder0=left%2;
        binary[digit]=remainder0;
        
        left=left/2;
        ++digit;
    }
    
    
    if (binary[0] && binary[1] && !binary[8] && !binary[9] && !binary[11]){
        
        return 0;
        
        
    }
    else{
        
        return 1;
    }
}


bool Fileloader::Process_region(const char * chr, int start, int end){
        
    phase_chrID=string(chr);
    phase_start=start;
    phase_end=end;
    
    refgenome=reader.Load_Chr(chr,phase_start,phase_end);
    
    snps_records.set_Region(phase_start);
    
    Locate_region();
    
    printf("Start reading the region %s : %d \n", chr, buff_start);
    
    Start_Read();

    Close();
    
    Write(phase_end+100);

    Barcodes_Records.Outputall();
    
    printf("\n\nFinish Reading Total Reads: %lld \n", line_index);
    
    return 0;
    
}

int Fileloader::LoadAlignment(BamTools::BamAlignment &alignment, Cigar_Reader &cigar, int* barindex, int pid){
    
    std::lock_guard<std::mutex> lk(file_read_lock);
    
    int ifget=0;
    while (!ifget){
        
        ifget=GetNextAlignment(alignment, cigar, barindex);
        
    }
    
    if (ifget<0) {return ifget;}
    
    Threads_Posi[pid]=alignment.Position;
    
    Check_Threads();
        
    if (!(++line_index%100000) || ifget>1)
        
    {
        Wait_Write(last_read_start);
    }
    
        
    Barcodes_Records.Update_Infor(alignment.Length,alignment.Position,alignment.Length+alignment.Position);
    
    if (alignment.Position+alignment.Length>=buff_end) buff_end=alignment.Position+alignment.Length+1;
    
    return ifget;
}



bool Fileloader::Wait_Write(int New_start){
    
    if (New_start>last_write_start){
        
        Write(New_start);
        write_cd.notify_all();
        
    }
    
    return 0;
}




bool Fileloader::Write(int write_end){
    
    cout<<"<<<<<<<<<< Processing Position: "<<write_end<<" >>>>>>>>>>"<<endl;
    
    snps_records.WriteSnp( write_end);
    
    Barcodes_Records.Write_Barcodes(write_end);
    
    return true;
    
}
