//
//  Barcodes_Record.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/25/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "barcodes_record.hpp"

using namespace std;



Barcode_Records::Barcode_Records(const char* file){
    
    
    Allbarcodes.reserve(sizeof(Molecule)*1000000);
    
    Molecule* newbarcode=new Molecule;
    
    newbarcode->barname=(char*)calloc(1,sizeof(char));
    newbarcode->barcode=0;
    newbarcode->firstposi=0;
    newbarcode->lastposi=0;
    newbarcode->molesize=0;
    newbarcode->molecounter=0;
    
    Allbarcodes.push_back(newbarcode);
    ifvalide_barcodes.push_back(1);
    
    barcode_hashcells[0]=0;
    barcode_counter=1;
    singlemole_counter=1;
	barname_counter=0;
    deleted_barcode=0;
    
    barcodefile=std::string(file)+string("_barcode");
}

unsigned int Barcode_Records::push_back(long long barcode){
    
    if (!barcode){
        
        return 0;
    }
    
    lock_guard<mutex> guard(barcode_record_lock);
    
    current_barcode=barcode;
    
    current_barindex=Record_Search(current_barcode);
    
    if (current_barindex){
        
        return current_barindex;
    }
    
    string bar_string="MI:"+to_string(barcode);
    
    current_barname=(char*)malloc((bar_string.length()+1)*(sizeof(char)));
    strcpy(current_barname, bar_string.c_str());
    
    Addbarcode();
    
    return current_barindex;
}

unsigned int Barcode_Records::push_back(const char* molename){
        
    int namesize=(int)strlen(molename);
    
    if (!namesize){
        
        return 0;
    }
    
    lock_guard<mutex> guard(barcode_record_lock);
    
    current_barcode=-hash(molename, namesize);
    current_barindex=Record_Search(current_barcode);
    
    if (current_barindex){
        
        return current_barindex;
    }
    
    current_barname=(char*)malloc((namesize+1)*sizeof(char));
    strcpy(current_barname, molename);
    
    Addbarcode();
    
    ++singlemole_counter;
        
    return current_barindex;
}


void Barcode_Records::Addbarcode(){
    
    current_barindex=barcode_counter++;
    
    Molecule* newbarcode=new Molecule;

    newbarcode->barname=current_barname;
    newbarcode->barcode=current_barcode;
    newbarcode->firstposi=0;
    newbarcode->lastposi=0;
    newbarcode->molesize=0;
    newbarcode->molecounter=0;
    
    Allbarcodes.push_back(newbarcode);
    ifvalide_barcodes.push_back(1);
    
    barcode_hashcells[current_barcode]=current_barindex;
    
}

unsigned int Barcode_Records::Record_Search(long long barcode_hash){

    
    auto findindex = barcode_hashcells.find(barcode_hash);
    
    if(findindex==barcode_hashcells.end()){
        return 0;
    }
    
    else{
        
        return findindex->second;
    }
    
}

void Barcode_Records::Update_Infor(int molesize, int start, int end){
    
    lock_guard<mutex> guard(barcode_record_lock);
    
    Molecule *current_mole=Allbarcodes[current_barindex-deleted_barcode];
    
    current_mole->molesize+=molesize;
    current_mole->lastposi=end;
    ++current_mole->molecounter;
    
    if (!current_mole->firstposi){current_mole->firstposi=start;}
}


void Barcode_Records::Addbarname(){

	Allbarnames[barname_counter]=current_barname;

	Allbarnames_barcodes[barname_counter]=current_barcode;

    ++barname_counter;
    
	if (barname_counter==1000){

		Write_Barnames();

	}

}

void Barcode_Records::Write_Barnames(){

    FILE *f = fopen(barnamefile.c_str(), "a");

    if (f == NULL){

        std::cerr << "ERROR: Could not open barname file " << barnamefile << std::endl;

    }

    for (unsigned int i=0; i<barname_counter; ++i) {

        fprintf(f,"%s\t%lld\n",Allbarnames[i], Allbarnames_barcodes[i]);

        free(Allbarnames[i]);

        Allbarnames[i]=NULL;

    }

    fclose(f);

	barname_counter=0;

}



bool Barcode_Records::Clear_Barcode(long long index){
    
    ifvalide_barcodes[index]=0;
    
    Molecule *current_mole=Allbarcodes[index];
        
    barcode_hashcells.erase(current_mole->barcode);
    
    free(current_mole->barname);
    
    delete current_mole;
    current_mole=NULL;
    
    return 1;
    
}


void Barcode_Records::Write_Barcodes(long long current_posi){
    
    FILE *f = fopen(barcodefile.c_str(), "a");
    
    if (f == NULL){
        
        fclose(f);
        std::cerr << "ERROR: Could not open file " << barcodefile << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    bool iffindmin=0;
    unsigned int new_deleted_barcode=1;
    Molecule* current_mole;

    for (unsigned int i=1; i<Allbarcodes.size(); ++i){
        
        if (!ifvalide_barcodes[i]) continue;
        
        current_mole=Allbarcodes[i];
        
        if (current_posi>current_mole->lastposi+max_Base_Buffer) {
            
            if (!iffindmin) {
                iffindmin=1;
                new_deleted_barcode=i;
            }
            
            continue;
        }
        
        fprintf(f,"%u\t%s\t%u\t%u\t%zu\t%u\n",i+deleted_barcode, current_mole->barname, current_mole->firstposi, current_mole->lastposi, current_mole->molecounter,  current_mole->molesize);
        
        Clear_Barcode(i);
    }
    
    fclose(f);
    
    
    if (new_deleted_barcode>1) {
        
        ifvalide_barcodes.erase(ifvalide_barcodes.begin()+1, ifvalide_barcodes.begin()+new_deleted_barcode);
        Allbarcodes.erase(Allbarcodes.begin()+1,Allbarcodes.begin()+new_deleted_barcode);
    }
    
    
    deleted_barcode+=new_deleted_barcode-1;
}


void Barcode_Records::Outputall(){
    
    FILE *f = fopen(barcodefile.c_str(), "a");
    
    if (f == NULL){
        std::cerr << "ERROR: Could not open file " << barcodefile << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    Molecule* current_mole;
    for (unsigned int i=1; i<Allbarcodes.size() ; ++i){
        
        if (!ifvalide_barcodes[i]) continue;
        
        current_mole=Allbarcodes[i];
        
        fprintf(f,"%u\t%s\t%u\t%u\t%zu\t%u\n",i+deleted_barcode, current_mole->barname, current_mole->firstposi, current_mole->lastposi, current_mole->molecounter,  current_mole->molesize);
    }
    
    fclose(f);
}


long long Barcode_Records::hash(const char *s, int len){
    
    long long h = 1125899906842597LL;
    
    for (int i=0;i<len;++i) {
        h = 31*h + s[i];
    }
    return h;
}

