#!/usr/bin/python

import re 
import pandas as pd
import os
import argparse
import collections as cl
import numpy as np
import multiprocessing as mul
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
import lpa_init_new as slpa


def flatten_list_of_lists(lst_of_lsts):
	N = sum(map(len, lst_of_lsts))    
	starts = np.empty(len(lst_of_lsts)+1, dtype=np.int32)  
	values = np.empty(N, dtype=np.int32)

	starts[0], cnt = 0, 0
	for i,lst in enumerate(lst_of_lsts):
		for el in map(np.int32,lst):
			values[cnt] = el
			cnt += 1      
		starts[i+1] = cnt  

	return starts, np.int32(values)
	
def readdata(inputfile):
	
	table=pd.read_csv(inputfile,header=None, sep='\t',dtype='str',keep_default_na=False).values.tolist()

	
	data=[[snp[0], map(int,snp[1].split(',')[:-1]),map(int,snp[2].split(',')[:-1])] for snp in table]
		
	return data	
	
class network_matrix:
	
	def __init__(self):
		
		self.snp_=[[[],[]]]
		self.snpindex={0:0}
		self.snplist=[0]
		
		self.mole_=cl.defaultdict(lambda: [set(), set()])
		self.mole_[0]
		self.moleindex={0:0}
		self.molelist=[0]
		

	
	def add_row(self,data):
		
		coordi=int(data[0])
				
		plen0=len(data[1])
		nlen0=len(data[2])
		data=data[1:]
							
		if plen0>0 and nlen0>0:
			
			self.snpindex[coordi]=len(self.snplist)
			self.snp_.append(data)
			self.snplist.append(coordi)
			
			for mole in data[0]:
				
				(self.mole_[mole])[0].add(coordi)
			
			for mole in data[1]:	
				
				(self.mole_[mole])[1].add(coordi)
	
		return self
	
	
	def endload(self):
		
		
		self.molelist=sorted(self.mole_.keys())
		
		self.moleindex={x:i for i, x in enumerate(self.molelist)}
		
		self.mole_=map(lambda x: map(list,x),[self.mole_[x] for x in self.molelist])
				
		return self
	
def loaddata(data):
	
	thematrix=network_matrix()
		
	map(thematrix.add_row, data)
	
	thematrix.endload()
	
	snpindex=thematrix.snpindex
	mole_matrix=[[[snpindex[x] for x in thelist[0]], [snpindex[x] for x in thelist[1]]] for thelist in thematrix.mole_]
	
	moleindex=thematrix.moleindex
	snp_matrix=[[[moleindex[x] for x in thelist[0]], [moleindex[x] for x in thelist[1]]] for thelist in thematrix.snp_]

	lsnp=len(snp_matrix)	
		
	p_mole,n_mole=zip(*mole_matrix)
	
	p_snp,n_snp=zip(*snp_matrix)
	
	p_mole_starts,p_mole_list=flatten_list_of_lists(p_mole)
		
	n_mole_starts,n_mole_list=flatten_list_of_lists(n_mole)
		
	p_snp_starts,p_snp_list=flatten_list_of_lists(p_snp)
		
	n_snp_starts,n_snp_list=flatten_list_of_lists(n_snp)

	init_barcode=np.array(thematrix.molelist, dtype=np.int32)
	
	init_coordi=np.array(thematrix.snplist, dtype=np.int32)

	return (p_mole_starts,p_mole_list,n_mole_starts,n_mole_list, init_barcode),(p_snp_starts,p_snp_list,n_snp_starts,n_snp_list,init_coordi)

def phase(paths):
			
	inputfile,outputfile_snp,outputfile_mole =paths
		
	data=readdata(inputfile)
	
	mole_data,snp_data=loaddata(data)
		
	init_coordi=snp_data[-1]
	init_barcode=mole_data[-1]
	
	init_labels=np.array(range(len(init_coordi)), dtype=np.int32)
	
	init_signs=np.ones(len(init_coordi), dtype=np.int32)
	
	SLPA_network=slpa.SLPA_network()
	
	SLPA_network.load_snp(*snp_data)
	
	SLPA_network.load_mole(*mole_data[:-1])
		
	SLPA_network.load_blk(init_labels, init_signs)
	
	num_labels=0
	for cutoff in [16,10,8,6,4]:
						
		SLPA_network.set_para(cutoff)
		
		SLPA_network.cycle(10,10)
			
		labels, signs=SLPA_network.output_snp()
		
		num_labels=len(set(labels))
						
		if num_labels==1:
			
			break
	
	while num_labels>1:
		
		SLPA_network.cycle(10,10)
			
		labels, signs=SLPA_network.output_snp()
		
		num_labels=len(set(labels))
				
		if num_labels>=len(set(labels)):
			
			break
	
		num_labels=len(set(labels))
		
	snp_phases=[(c,a*b) for a,b,c in zip(labels, signs,init_coordi)[1:]]
		
	mole_phases=[(y,",".join(map(str,x))) for x,y in zip(SLPA_network.output_mole(),init_barcode[1:])]
						
	pd.DataFrame(snp_phases).to_csv(outputfile_snp,mode='w',header=None,index=False,sep='\t')
	
	pd.DataFrame(mole_phases).to_csv(outputfile_mole,mode='w',header=None,index=False,sep='\t')
	
	return 0
	
def phase_allfiles(args,threads):
		
	p=mul.Pool(processes=threads)
	
	p.map_async(phase,args)
	
	p.close()
	p.join()

def combine_allfiles(allfiles,outputfile_snp, outputfile_mole):
		
	block_counter=1
		
	for infile, snpfile, molefile in allfiles:
		
		snp_results=pd.read_csv(snpfile,header=None, sep='\t').fillna(0)
		
		unique_labels=set(map(abs,list(snp_results[1])))
	
		label_to_block={0:0}
		
		for x in unique_labels:
			
			label_to_block[x]=len(label_to_block)-1+block_counter
			
		block_counter+=len(label_to_block)-1
		
		out=[label_to_block[abs(x)] if x>0 else -label_to_block[abs(x)] for x in list(snp_results[1])]
		
		out=zip(snp_results[0],out)
		
		pd.DataFrame(out).to_csv(outputfile_snp,mode='a',header=None,index=False,sep='\t')
	
		mole_results=pd.read_csv(molefile,header=None, sep='\t',dtype='str').fillna('0')
		
		out=[','.join(map(lambda x: str(label_to_block[abs(x)]) if x>0 else str(-label_to_block[abs(x)]), map(int,str(y).split(',')))) for y in list(mole_results[1])]

		out=zip(mole_results[0],out)
		
		pd.DataFrame(out).to_csv(outputfile_mole,mode='a',header=None,index=False,sep='\t')
		

def run(args):
	
	inputfolder=args.input
	if inputfolder[-1]!='/':
		inputfolder+='/'
						
	outfolder=args.output
	if outfolder[-1]!='/':
		outfolder+='/'
	
	threads=args.threads
	
	allfiles=[[inputfolder+x,outfolder+x+'_phased_snp',outfolder+x+'_phased_mole'] for x in sorted([x for x in os.listdir(inputfolder) if "block_" in x], key=lambda x: int(x.split('_')[-1]))]
	
	phase_allfiles(allfiles,threads)
	
	combine_allfiles(allfiles,outfolder+'phased_snps',outfolder+'phased_moles')
	

def main():
		
	parser=argparse.ArgumentParser(description="algthorithm to phase all blocks")
	parser.add_argument("-i","--inputfolder",help="input data folder" ,dest="input",type=str,required=True)
	parser.add_argument("-o","--outputfolder",help="output path" ,dest="output",type=str,required=True)
	parser.add_argument("-t","--threads",help="number of threads" ,dest="threads",type=int,default=1)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)	

	
if __name__=='__main__':
	
	main()

