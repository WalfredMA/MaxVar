//
//  main.c
//  findbreaks
//
//  Created by Wangfei MA on 8/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>


#define allocsize 1000000
#define hash1 999
#define hash2 998

typedef struct{
    int size;
    int* values;
    int* indexs;
}hash;


typedef struct{
    
    int num_snp;
    int* snps;
    int num_blk;
    int** blks;
}cluster;

typedef struct{
    
    int *icluster;
}barcode;


int barcode_counter=0;
int cluster_counter=0;
int snp_counter=0;
barcode* allbarcodes;
cluster* allclusters;
int* coordinates;
char** allLines;
hash barcode_hash[1000][1000];
char StrLine[1000000];
int barcodes[1000000];
int barcodes_cluster[1000000];
int barcodes_index[1000000];

int compare_rev( const void* a, const void* b)
{
    int int_a = * ( (int*) a );
    int int_b = * ( (int*) b );
    
    return  (int_a < int_b)-(int_a > int_b);
}

int compare( const void* a, const void* b)
{
    int int_a = * ( (int*) a );
    int int_b = * ( (int*) b );
    
    return  (int_a > int_b)-(int_a < int_b);
}

int sum(int* list, int num){
    
    int thesum=0;
    for (int i=0;i<num;i++){
        
        thesum+=list[i];
    }
    return thesum;
}

int min(int *list, int size)
{
    int min=0;
    for (int i=0;i<size;i++)
    {
        if (list[i]<min) min=list[i];
    }
    return min;
    
}

int max(int *list, int size)
{
    int max=0;
    for (int i=0;i<size;i++)
    {
        if (list[i]>max) max=list[i];
    }
    return max;
    
}

int N50(int* list, int list_size){
    
    int iffind=0;
    int thesum=sum(list, list_size);
    
    qsort( list, list_size, sizeof(int), compare_rev);
    
    int left=1,right=list_size, current, current_sum, current_subone;
    int ifreachN50, ifpassN50;
    
    current=(left+right)/2;
    while (!iffind) {
        
        current_sum=sum(list, current);
        current_subone=current_sum-list[current-1];
        
        ifreachN50=2*current_sum-thesum;
        ifpassN50=2*current_subone-thesum;
        
        if (ifreachN50>=0 && ifpassN50<0){
            
            iffind=current;
            break;
            
        }
        
        else if (!ifpassN50) {
            
            iffind=current-1;
            break;
            
        }
        
        else if (ifpassN50>0){
            
            right=current;
            current=(right+left)/2;
        }
        
        else if (ifpassN50<0){
            
            left=current;
            if (right!=left+1){
                
                current=(current+right)/2;
            }
            
            else{
                
                iffind=right;
                break;
            }
        }
    }
    
    return list[iffind];
    
}



void newbarcode(){
    
    barcode_counter++;
    
    if ((barcode_counter-allocsize)>=0 && !((barcode_counter)%100)) {
        
        allbarcodes=(barcode*)realloc(allbarcodes, (barcode_counter+100)*sizeof(barcode));
        
    }
    
    allbarcodes[barcode_counter].icluster=NULL;
    
}

void newcluster(){
    
    cluster_counter++;
    
    if ((cluster_counter-allocsize)>=0 && !((cluster_counter)%100)) {
        
        allclusters=(cluster*)realloc(allclusters, (cluster_counter+100)*sizeof(cluster));
        
    }
    
    allclusters[cluster_counter].num_snp=0;
    allclusters[cluster_counter].num_blk=1;
    allclusters[cluster_counter].blks=(int**)malloc(sizeof(int**));
    allclusters[cluster_counter].blks[0]=(int*)malloc(sizeof(int*));
    *allclusters[cluster_counter].blks[0]=cluster_counter;
}


void blockaddsnp(int snp, int num_bar, int* bars, int cluster_index){
    
    cluster* cluster=&allclusters[cluster_index];
    
    if (!(cluster->num_snp))
    {
        
        cluster->snps= (int*)malloc(100*sizeof(int));
        
    }
    
    
    else if (!(cluster->num_snp%100))
    {
        
        cluster->snps= (int*)realloc(cluster->snps,  (cluster->num_snp+100)*sizeof(int));
        
    }
    
    cluster->snps[cluster->num_snp++]=snp;
    
    for (int i=0; i<num_bar; i++)
    {
        
        if (!allbarcodes[bars[i]].icluster)
        {
            allbarcodes[bars[i]].icluster=cluster->blks[0];
                        
        }
        
    }
    
}


void blockaddblock( int old_index,int new_index)
{
    
    cluster* old_cluster=&allclusters[old_index];
    cluster* new_cluster=&allclusters[new_index];
    
    new_cluster->blks= realloc(new_cluster->blks,sizeof(int*)*(new_cluster->num_blk+old_cluster->num_blk));
    
    for (int i=0; i<old_cluster->num_blk; i++)
    {
        *old_cluster->blks[i]=new_index;
        new_cluster->blks[new_cluster->num_blk++]=old_cluster->blks[i];
    }
    
    for (int i=0; i<old_cluster->num_snp; i++)
    {
        blockaddsnp(old_cluster->snps[i],0,0,new_index);
    }
    
    if (old_cluster->num_snp) {
        free(old_cluster->snps);
        old_cluster->num_snp=0;
    }
    
    free(old_cluster->blks);
    
}


int findindex(int barcode, hash* hash_cell){
    
    
    int cell_size=hash_cell->size;
    
    for (int i=0;i<cell_size;i++){
        
        if (barcode==hash_cell->values[i])
        {
            return hash_cell->indexs[i];
        }
    }
    
    newbarcode();
    
    hash_cell->size++;
    
    if (!cell_size){
        
        hash_cell->indexs=(int*)malloc(100*sizeof(int));
        hash_cell->values=(int*)malloc(100*sizeof(int));
    }
    
    else if (!(cell_size%100)){
        
        hash_cell->indexs=(int*)realloc(hash_cell->indexs, (cell_size+100)*sizeof(int));
        hash_cell->values=(int*)realloc(hash_cell->values, (cell_size+100)*sizeof(int));
        
    }
    
    hash_cell->values[cell_size]=barcode;
    hash_cell->indexs[cell_size]=barcode_counter;
    
    return 0;
}

int linkclusters(int num_counter)
{
    
    int barindex, barcode;
    int hash_1, hash_2;
    hash* hash_cell;
    
    int unique_icluster=0, ifunique=0;
    int *tempint=NULL;
    
    for (int k=0;k<num_counter;k++){
        
        barcode=barcodes[k];
        
        hash_1=(barcode+1000)%hash1;
        hash_2=(barcode+1000)%hash2;
        
        hash_cell=&barcode_hash[hash_1][hash_2];
        
        barindex=findindex(barcode, hash_cell);
        
        if (!barindex)
        {
            
            barcodes_index[k]=barcode_counter;
            
        }
        
        else {
            
            barcodes_index[k]=barindex;
            
            tempint=allbarcodes[barindex].icluster;
            
            if (!tempint) continue;
            
            ifunique=1;
            for (int i=0;i<unique_icluster;i++)
            {
                
                if (barcodes_cluster[i]==*tempint)
                {
                    
                    ifunique=0;
                    break;
                }
            }
            
            if (ifunique)
            {
                barcodes_cluster[unique_icluster++]=*tempint;
            }
        }
    }
    
    return unique_icluster;
}

void combine_clusters(int num_counter, int unique_icluster)
{
    
    int min_cluster=cluster_counter+1;
    
    for (int k=0; k<unique_icluster; k++)
    {
        
        if (barcodes_cluster[k]<min_cluster) min_cluster=barcodes_cluster[k];
        
    }
    
    
    if (!unique_icluster){
        
        
        newcluster();
        
        
        blockaddsnp(snp_counter,num_counter, barcodes_index, cluster_counter);
        
        
    }
    
    else {
        
        int index;
        for (int k=0; k<unique_icluster; k++){
            
            index=barcodes_cluster[k];
            
            if (index==min_cluster) continue;
            
            blockaddblock (index, min_cluster);
            
        }
        
        blockaddsnp(snp_counter,num_counter, barcodes_index, min_cluster);
        
    }
    
}

int blocking(const char *StrLine){
    
    if (!(snp_counter%10000)){
        
        printf("%d\n",snp_counter);
    }
    
    if (snp_counter>=allocsize && !(snp_counter%100)){
        
        coordinates=(int*) realloc(coordinates, (snp_counter+100)*sizeof(int));
        
    }
    
    int strlen0=(int) strlen(StrLine);
    
    int coordinate=0,curren_num=0,index=0,num_counter=0,barcode0=0;
    
    char current_chr=StrLine[index++];
    
    while (current_chr!='\t'){
        
        coordinate*=10;
        coordinate+=(int) atoi(&current_chr);
        current_chr=StrLine[index++];
    }
    
    coordinate+=1;
    coordinates[snp_counter]=coordinate;
    
    for (;index< strlen0;index++){
        
        current_chr=StrLine[index];
        
        if (current_chr>='0'&&current_chr<='9')
        {
            
            curren_num*=10;
            curren_num+=(int) atoi(&current_chr);
            
        }
        
        else
        {
            
            if (curren_num) barcodes[num_counter++]=curren_num;
            
            curren_num=0;
            
        }
        
    }
    
    int unique_icluster=linkclusters(num_counter);
    
    combine_clusters(num_counter, unique_icluster);
    
    return strlen0;
    
}


void indexfile(const char inputfile[]) {
    
    FILE *fpread=fopen(inputfile,"r");
    
    if (fpread==NULL ){
        
        printf("cannot open file\n");
        
        fclose(fpread);
        
        exit(1);
        
    }
    
    int strlen0;
    
    while (fgets(StrLine,1000000,fpread)!=NULL && !feof(fpread))
        
    {
        snp_counter++;
        
        if ((snp_counter>=1000000)&&!(snp_counter%100)){
            
            allLines=(char**)realloc(allLines,  (snp_counter+100)*sizeof(char*));
        }
        
        strlen0=blocking(StrLine);
        
        allLines[snp_counter]=(char*)malloc((strlen0+2)*sizeof(char));
        strncpy(allLines[snp_counter], StrLine, strlen0*sizeof(char));
        allLines[snp_counter][strlen0]='\0';
        
    }
    
    fclose(fpread);
    
}



void output(const char outputpath[]){
    
    printf("writing files\n");
    
    int pathlen=(int)(strlen(outputpath));
    char blockfile[pathlen+100];
    strncpy(blockfile, outputpath, pathlen);
    
    int block_sizes[cluster_counter];
    block_sizes[0]=0;
    
    cluster* current_cluster;
    int eachline=0;
    int valide_cluster=1;
    int cluster_size;
    
    int totalline=0;
    for (int icluster=1; icluster<cluster_counter; icluster++){
        
        block_sizes[icluster]=0;
        
        current_cluster=&allclusters[icluster];
        
        cluster_size=current_cluster->num_snp;
        
        totalline+=cluster_size;
        if (!cluster_size)
        {
            block_sizes[icluster]=0;
            continue;
        }
        
        qsort( current_cluster->snps, cluster_size, sizeof(int), compare );
        
        
        
        block_sizes[icluster]= coordinates[current_cluster->snps[cluster_size-1]]-coordinates[current_cluster->snps[0]];
        
        if (block_sizes[icluster]>block_sizes[0])
        {
            block_sizes[0]=block_sizes[icluster];
        }
        
        snprintf(blockfile+pathlen, 50, "block_%d",valide_cluster++);
        
        FILE *fwrite1=fopen(blockfile, "w");
        
        if (fwrite1!=NULL){
            
            char* Linestr;
            for (int i=0;i<cluster_size;i++){
                
                eachline=current_cluster->snps[i];
                Linestr=allLines[eachline];
                
                fprintf(fwrite1,"%s", Linestr);
                
            }
        }
        
        fclose(fwrite1);
    }
    
    
    //printf("N50 of all blocks: %d kb\nLargest block: %d kb\ntotal snp: %d kb\n", N50(block_sizes, cluster_counter)/1000, block_sizes[0]/1000, snp_counter/1000);
    
}




int main(int argc, const char ** argv) {
    
    printf("Indexing %s\n", argv[1]);
    
    const char *inputfile=argv[1];
    const char *outputpath=argv[2];
    
    
    allclusters=(cluster*)malloc((allocsize)*sizeof(cluster));
    allbarcodes=(barcode*)malloc((allocsize)*sizeof(barcode));
    allLines=(char**)malloc((allocsize)*sizeof(char*));
    coordinates=(int*)malloc((allocsize)*sizeof(int));
    
    for (int i=0;i<1000;i++){
        
        for (int j=0;j<1000;j++){
            
            barcode_hash[i][j].size=0;
            
        }
    }
    
    indexfile(inputfile);
    
    output(outputpath);
    
    printf("Converting finished\n");
    
    return 0;
}
