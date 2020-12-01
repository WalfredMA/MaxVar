#!/usr/bin/python

import re 
import pandas as pd
import os
import argparse
import collections as cl
import numpy as np
import math
import multiprocessing as mul
import ctypes
import gc
from scipy.stats import ranksums
from scipy.stats import fisher_exact


manager=mul.Manager()


file1="/Users/walfred/Desktop/NA19238_9_hetero.txt"
file1="/Users/walfred/Desktop/NA19238_9_zygo.txt"

parents={}
with open(file1,mode='r') as f:
	
	counter=0
	for line in f:
		counter+=1
		if counter>10000000:
			break
		split=line.split()
		
		alleles=split[1:4]
		#if alleles[1]==alleles[2]:
		parents[int(split[0])]=split[1:]
f.close()

global confi_cut,processor,refdata,mole_data,snpsites


def read_mole(molefile, snpfile):
	global snpsites
	label_snpcount=cl.defaultdict(list)
	snpsites=set([])
	with open(snpfile,mode='r') as f:
		
		for line in f:
			
			split=line.split()
			coordinate=int(split[0])
			labels=map(int,split[1].split(','))
			
			if len(labels) and 0 not in labels:
				snpsites.add(coordinate)
			for label in labels:
				label_snpcount[abs(label)].append(coordinate)
	f.close()
	
	label_snpcount=set([0]+[abs(label) for label,coordis in label_snpcount.items() if len(coordis)>10 and max(coordis)-min(coordis)>10000])
	
	data={}
	with open(molefile,mode='r') as f:
		
		for line in f:
			
			split=line.split()
			data[int(split[0])]=np.array([x for x in map(int,split[1].split(',')) if abs(x) in label_snpcount],dtype=np.int32)
			
	f.close()	
	
	return data

def read_ref(reffile, chrom):
	
	with open(reffile,mode='r') as f:
		read=f.read().split('>')[1:]
	f.close()
	
	out=""
	for x in read:
		
		if x.splitlines()[0].split()[0]==chrom:
			
			out=''.join(x.splitlines()[1:])
	
	return out

def allele_pair(x1,x2,ifphase,phase):
			
	if phase[0]<phase[1]:
		genotype=x1[0]+ifphase+x2[0]
		bx=x1[-1]+";"+x2[-1]
	elif phase[0]>phase[1]:
		genotype=x2[0]+ifphase+x1[0]
		bx=x2[-1]+";"+x1[-1]
	else:
		genotype=x1[0]+"/"+x2[0]
		bx=x1[-1]+";"+x2[-1]
	
	return ":".join([genotype]+x1[1:4]+map(lambda x:x[0]+','+x[1],zip(x1,x2)[4:6])+x1[6:-1]+[bx])

def mean(x):
	
	if len(x):
		
		return sum(x)/len(x)
	else:
		return 0
	

class bionomial:
	
	def __init__(self,error=0.5):
		
		self.log2_add = [100, 100, 99, 99, 98, 98, 97, 97, 96, 96, 95, 95, 94, 94, 93, 93, 92, 92, 91, 91, 90, 90, 89, 89, 88, 88, 88, 87, 87, 86, 86, 85, 85, 84, 84, 84, 83, 83, 82, 82, 81, 81, 81, 80, 80, 79, 79, 78, 78, 78, 77, 77, 76, 76, 76, 75, 75, 74, 74, 73, 73, 73, 72, 72, 72, 71, 71, 70, 70, 70, 69, 69, 68, 68, 68, 67, 67, 67, 66, 66, 65, 65, 65, 64, 64, 64, 63, 63, 63, 62, 62, 62, 61, 61, 61, 60, 60, 60, 59, 59, 58, 58, 58, 58, 57, 57, 57, 56, 56, 56, 55, 55, 55, 54, 54, 54, 53, 53, 53, 52, 52, 52, 52, 51, 51, 51, 50, 50, 50, 49, 49, 49, 49, 48, 48, 48, 47, 47, 47, 47, 46, 46, 46, 46, 45, 45, 45, 44, 44, 44, 44, 43, 43, 43, 43, 42, 42, 42, 42, 41, 41, 41, 41, 40, 40, 40, 40, 39, 39, 39, 39, 38, 38, 38, 38, 38, 37, 37, 37, 37, 36, 36, 36, 36, 36, 35, 35, 35, 35, 34, 34, 34, 34, 34, 33, 33, 33, 33, 33, 32, 32, 32, 32, 32, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29, 29, 29, 29, 29, 29, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 26, 26, 26, 26, 26, 26, 25, 25, 25, 25, 25, 25, 24, 24, 24, 24, 24, 24, 23, 23, 23, 23, 23, 23, 23, 22, 22, 22, 22, 22, 22, 22, 21, 21, 21, 21, 21, 21, 21, 21, 20, 20, 20, 20, 20, 20, 20, 19, 19, 19, 19, 19, 19, 19, 19, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
		
		self.cdf_storage=manager.dict({})
						
		self.pdf_storage=manager.dict({})
		
		self.exp_error=np.log2(error)
		
		self.exp_noerror=np.log2(1-error)

	def pdf(self,total,obs,factor=1):
					
		if (total,obs) in self.pdf_storage:
			
			return self.pdf_storage[(total,obs)]
	
		else:
			
			value=0
			for x in xrange(total-obs+1,total+1):
				
				value+=np.log2(x)
			
			for x in xrange(1,obs+1):
				
				value-=np.log2(x)
			
			value+=self.exp_error*obs
			value+=self.exp_noerror*(total-obs)
			
			value=int(value*factor+0.5)
		
			self.pdf_storage[(total,obs)]=value
		
			return value
		
	def cdf(self, total, obs):
		
		if (obs, total) in self.cdf_storage:
			
			return self.cdf_storage[(obs, total)]
	
		else:
			if obs<total:
				
				results=sorted([self.pdf(total,i,100) for i in xrange(obs,total+1)],reverse=True)
				
				result=results[0]
				for x in results[1:]:
					
					dis=result-x
					if dis>=len(self.log2_add):
						break
					
					result+=self.log2_add[dis]
				
				result*=0.01
			
			else:
				
				result=self.exp_noerror*(total-obs)+self.exp_error*obs
						
			self.cdf_storage[(obs, total)]=result
		
			return result


global gcounter
gcounter=[0,0,0]

class line_processor:
	
	def __init__(self):
		
		self.splitcomma=re.compile(r"(?<!,),")
		self.last_coordi=0
		
		self.ifvcf=mul.Value("i",1)
		self.coverage=mul.Value("i",0)
		self.chrom=mul.Value(ctypes.c_char_p,"")
		self.sample=mul.Value(ctypes.c_char_p,"")
		self.errorrates=manager.dict({})
		self.cutoffs=manager.dict({})
		self.qualcuts=manager.dict({})
		
		global confi_cut
		
		self.confi_cut=confi_cut
		self.indel_cut=2*confi_cut
				
	def load_phasedata(self, ifvcf, chrom, sample, infofile, reffile, molefile):
			
		self.ifvcf.value=ifvcf
			
		self.chrom.value=chrom
		if chrom[:3]=='chr':
			self.chrom.value=chrom[3:]
			
		self.sample.value=sample
									
		with open(infofile,mode='r') as f:
			readtext=f.read().splitlines()
		f.close()
		

		
		cutoffs={}
		errorrates={}
		qualcuts={}
				
		self.coverage.value=int(readtext[2])/int(readtext[1])
		
		for line in readtext:
			
			splits=line.split()
			
			findcutoff=[x for x in splits if "lowconfi_cutoff:" in x]
			
			if len(findcutoff):
				
				cutoffs[splits[0]]=int(findcutoff[0].split('lowconfi_cutoff:')[1])
			
			qualcutoff=[x for x in splits if "lowcut:" in x]
			
			if len(qualcutoff):
				
				qualcuts[splits[0]]=float(qualcutoff[0].split('lowcut:')[1])
				
			error_rate=[x for x in splits if "error_rate:" in x]
			
			if len(error_rate):
				
				errorrates[splits[0]]=float(error_rate[0].split("error_rate:")[1])
				
		self.errorrates.update(dict(errorrates))
		self.cutoffs.update(dict(cutoffs))
		self.qualcuts.update(qualcuts)
		
	
		
	
	def score(self,coordinate,valide_alleles,phased_alleles,barcodes,barcode_phases,del_sizes):
			
		if phased_alleles>1:
			phase_score=sum([20*abs(barcode_phases[allele][0]-barcode_phases[allele][1]) for allele in phased_alleles])
		else:
			phase_score=0
		
		sv_score=0
		allele_scores=[]
		for allele in valide_alleles[1:]:
			
			confi=self.errorrates.get("ATCGDNMI"[min(7,valide_alleles[0])]+'-'+"ATCGDNMI"[min(7,allele)],[0,0,0,0,0,0,0,0,0,0])[self.consect_index(coordinate,"ATCGDNMI"[min(7,valide_alleles[0])])]
			
			if confi:
				confi=-math.log10(confi)
			else:
				confi=-math.log10(0.01)
			
			
			
			cutoff=self.cutoffs.get("ATCGDNMI"[min(7,valide_alleles[0])]+'-'+"ATCGDNMI"[min(7,allele)],[0,0,0,0,0,0,0,0,0,0])[self.consect_index(coordinate,"ATCGDNMI"[min(7,valide_alleles[0])])]
			
		
			allele_score=(len([confi for x in barcodes[allele] if x>0])-cutoff)*confi+20
		
			allele_scores.append(allele_score)
			
			if allele==4:
								
				sv_score=sum([min(100,10*(del_size-10)) for del_size in del_sizes if del_size>10])
		
		return round((max(allele_scores)if allele_scores else 0)+sv_score+phase_score,2)	
	
	def phase_label(self,columns):
		
		global refdata,mole_data
		
		label_counter=cl.defaultdict(list)
		label_counter[0]=[]
		for i,column in enumerate(columns):
			
			for mole in column:
			
				labels=mole_data.get(abs(mole),[])
									
				for label in labels:
					
					label_counter[label].append(i)					
				
		max_label=0
		max_count=0
		for label in set(map(abs,label_counter.keys())):
			
			count=len(label_counter[label])+len(label_counter[-label])
			
			if count > max_count:
				
				max_label=label
				max_count=count
		
		allele_counter=cl.defaultdict(int)
		for allele in label_counter[max_label]:
			
			allele_counter[allele]+=1
			
		posi_accuracy=[max(allele_counter.values()+[0]),sum(allele_counter.values())]
		
		allele_counter=cl.defaultdict(int)
		for allele in label_counter[-max_label]:
			
			allele_counter[allele]+=1
		
		nega_accuracy=[max(allele_counter.values()+[0]),sum(allele_counter.values())]
		
		barcode_phases=[[0,0] for x in columns]
		if max_label:
			
			for i,column in enumerate(columns):
				
				for mole in column:
				
					labels=mole_data.get(abs(mole),[])
					
					for label in labels:
						
						if max_label==abs(label):	
							
							if label>0:
								barcode_phases[i][0]+=1
							else:
								barcode_phases[i][1]+=1
		
		alleles_phases=[0 for column in columns]
		for i in xrange(len(barcode_phases)):
			
			if i==4:
				phase_cutoff=self.indel_cut/2
			else:
				phase_cutoff=self.confi_cut/2
			
			if barcode_phases[i][0]>max(phase_cutoff,2*barcode_phases[i][1]):
				
				alleles_phases[i]=1
				
			elif barcode_phases[i][1]>max(phase_cutoff,2*barcode_phases[i][0]):
				alleles_phases[i]=-1
			
			#elif min(barcode_phases[i][0],barcode_phases[i][1])>max(phase_cutoff,(barcode_phases[i][0]+barcode_phases[i][1])/4):
				#alleles_phases[i]=0
			
		
		return max_label,barcode_phases,alleles_phases,[posi_accuracy, nega_accuracy]
	
	def consect_index(self, coordinate, base):
		
		global refdata
		if base not in ['A','T','C','G']:
			return 0
		
		front_consect=0		
		for ref_coordi in range(coordinate-1,max(0,coordinate-4),-1):
			
			ref_base=refdata[ref_coordi-1]
			if ref_base==base.upper() or ref_base==base.lower():
				front_consect+=1
			else:
				break
		
		back_consect=0
		for ref_coordi in range(coordinate+1,coordinate+4):
			
			ref_base=refdata[ref_coordi-1]
			if ref_base==base.upper() or ref_base==base.lower():
				back_consect+=1
			else:
				break
		
						
		if front_consect>2 or back_consect>2:
			return 9
		
		return front_consect*3+back_consect
	
	def read_alleles(self, splits):
			
		barcodes, quals, del_sizes=[],[],[]
		
		inserts=cl.defaultdict(list)
		inserts_qual=cl.defaultdict(list)
		for i,column in enumerate(splits[:7]):
							
			barcodes.append([])
			quals.append([])
			
			if not len(column):
				continue
			
			for mole in column:
				
				if i==5:
					continue
				
				barcodes[-1].append(int(mole[2:].split(":")[0]))
				
				quals[-1].append(ord(mole[0])-32)
					
				if i==4:
										
					del_sizes.append(mole[2:].split(':')[1])
					
		if len(splits)>7:
		
			for mole in column:
									
				genotype=mole[2:].split(":")[1]
				
				barcode=mole[2:].split(":")[0]
				
				if len(genotype)>2:
					genotype=genotype[1]+genotype[-1]
				
				if int(barcode) and len(genotype):
					
					inserts[genotype].append(int(barcode))
					inserts_qual[genotype].append(ord(mole[0])-32)	
				
		alleles=['A','T','C','G','-','N','M']+inserts.keys()
		
		barcodes.extend(inserts.values())
		
		quals.extend([inserts_qual[geno] for geno in inserts.keys()])
		
		return alleles, barcodes, quals, del_sizes
		
	def clean_redundent(self, columns, quals, del_sizes):
		
		barcodes=cl.defaultdict(lambda: [0 for i in xrange(len(columns))]) 
		barcodes_quals=cl.defaultdict(lambda: [0 for i in xrange(len(columns))]) 
				
		for i,(column,qual) in enumerate(zip(columns,quals)):
			
			for barcode,bar_qual in zip(column,qual):
				
				barcodes[abs(barcode)][i]+=1
				barcodes_quals[abs(barcode)][i]+=bar_qual
		
		barcodes_uniq={bar:allele_count.index(max(allele_count)) for bar,allele_count in barcodes.items() if max(allele_count)>=0.9*sum(allele_count)}
		
		for i,(column,qual) in enumerate(zip(columns,quals)):
			
			if i==4:
												
				filtered=[[bar, barcodes_quals[abs(bar)][i]/barcodes[abs(bar)][i],sum(map(int,size.split('-')))] for bar,size in zip(column,del_sizes) if barcodes_uniq.get(abs(bar),-1)==4]
								
				if len(filtered):
					columns[i],quals[i],del_sizes=map(list,zip(*filtered))
										
				else:
					columns[i],quals[i],del_sizes=[],[],[]
						
			else:
				
				filtered=[[bar,barcodes_quals[abs(bar)][i]/barcodes[abs(bar)][i]] for bar in column if barcodes_uniq.get(abs(bar),-1)==i]
								
				if len(filtered):
					columns[i],quals[i]=map(list,zip(*filtered))
										
				else:
					columns[i],quals[i]=[],[] 
		
		
		return columns, quals, del_sizes
	
	def readline(self,splits):		 
		
		alleles, barcodes, quals, del_sizes=self.read_alleles(splits)
		
		barcodes, quals, del_sizes=self.clean_redundent(barcodes, quals,del_sizes)
						
		return alleles, barcodes,quals, del_sizes
	
	def findalleles(self, columns, mean_quals, coordinate, alleles_counts,alleles_cutoff):
		
		global refdata,mole_data
		
		total_count=sum(alleles_counts)
		
		ref_base=refdata[abs(coordinate)-1].upper()
		
		alleles_index=sorted([i for i in xrange(len(columns)) if len(columns[i])>alleles_cutoff and mean_quals[i]>self.qualcuts["ATCGDNMI"[min(7,i)]]], key=lambda i:len(columns[i]),reverse=True)
				
		if coordinate<0 or len(columns[4])<self.indel_cut:
			
			alleles_index=[x for x in alleles_index if x!=4]
		
		if len(alleles_index)==1 and "ATCGDNMI"[min(7,alleles_index[0])]==ref_base:
						
			return []
				
		return alleles_index
		
	def valide_alleles(self, coordinate,columns,alleles_index,alleles_cutoff, mean_quals):
		
		global refdata,mole_data
		
		ref_base_f=refdata[abs(coordinate)-2].upper()
		ref_base_b=refdata[abs(coordinate)].upper() 
		
		alleles_counts=map(len,columns)
		 
		columns_posi=[int(sum([1 if x>0 else 0.5 for x in column ])) if (i==4 or (i<4 and "ATCGDNMI"[i]==ref_base_f and "ATCGDNMI"[i]==ref_base_b)) else len(column) for i,column in enumerate(columns)]
		
		total_posi=sum(columns_posi)
					
		valide_alleles=[]
				
			
		cutoff1=self.cutoffs.get("ATCGDNMI"[min(7,alleles_index[0])]+'-'+"ATCGDNMI"[min(7,alleles_index[1])],[0,0,0,0,0,0,0,0,0,0])[self.consect_index(abs(coordinate),"ATCGDNMI"[min(7,alleles_index[0])])]
		
		
		if columns_posi[alleles_index[1]]>=max(self.confi_cut,cutoff1) and alleles_counts[alleles_index[1]]>alleles_cutoff and mean_quals[alleles_index[1]]>=self.qualcuts["ATCGDNMI"[min(7,alleles_index[1])]]:
			
			valide_alleles.append(alleles_index[1])
		
				
		cutoff2=self.cutoffs.get("ATCGDNMI"[min(7,alleles_index[1])]+'-'+"ATCGDNMI"[min(7,alleles_index[0])],[0,0,0,0,0,0,0,0,0,0])[self.consect_index(abs(coordinate),"ATCGDNMI"[min(7,alleles_index[1])])]
		
		if columns_posi[alleles_index[0]]>=max(self.confi_cut,cutoff2) and alleles_counts[alleles_index[0]]>alleles_cutoff and mean_quals[alleles_index[0]]>=self.qualcuts["ATCGDNMI"[min(7,alleles_index[0])]]:
			
			valide_alleles.append(alleles_index[0])
			
					
		return valide_alleles
		
	
	def process_data(self, line):
		
		global bionom,refdata,mole_data,gcounter
		
		splits=line.split("\t")
		
		coordinate=int(splits[0])
		
		splits=[re.split(self.splitcomma,","+column)[1:-1] for column in splits[1:-1]]

		alleles_counts=map(len, splits)
		
		alleles_sum=sum(alleles_counts)
		if coordinate<0:
			alleles_sum-=alleles_counts[4]
		
		alleles_cutoff=max(self.confi_cut,max(alleles_counts)/3,self.coverage.value/4)		
								
		alleles, barcodes, quals, del_sizes=self.readline(splits)
		
		alleles_counts=map(len, barcodes)
		
		mean_quals=[2.0*sum(sorted(qual[len(qual)/2:]))/(len(qual)+1) if qual else 0 for qual in quals]			
				
		alleles_found=self.findalleles(barcodes, mean_quals, coordinate, alleles_counts,alleles_cutoff)
		
		if not len(alleles_found) or len(alleles_found)>2:

			return 0,[]			
		
		allele_data=[barcodes[i] for i in alleles_found]
		
		barcodes=[x if i in alleles_found else [] for i,x in enumerate(barcodes)]
		
		max_label,barcode_phases,alleles_phase,allele_accuracy=self.phase_label(barcodes)
		
		allele_quals=[quals[i] for i in alleles_found]
		#if min([-error_bionom.cdf(x[1], x[1]-x[0]) for x in allele_accuracy if x[1]>5 and x[0]<x[1]*0.5]+[0])>7:
		
							
		if len(alleles_found)==1:
			
			if alleles_counts[alleles_found[0]]>max(self.confi_cut, sum(alleles_counts)*0.8):
				confidence=1
			else:
				return 0,[]	
			
			valide_alleles=alleles_found
			
			score=self.score(abs(coordinate),alleles_found,alleles_found,barcodes,barcode_phases,del_sizes)
			
			return confidence,[coordinate,score,max_label,zip(*[[alleles[index] for index in valide_alleles],[barcodes[index] for index in valide_alleles],[quals[index] for index in valide_alleles],[barcode_phases[index] for index in valide_alleles],[alleles_phase[allele] for allele in valide_alleles],valide_alleles])]
			
		reverse_phase_total=max(sum([barcode_phases[x][i%2] for i,x in enumerate(alleles_found)]), sum([barcode_phases[x][(i+1)%2] for i,x in enumerate(alleles_found)]))
		
		phase_total=[abs(barcode_phases[x][0]-barcode_phases[x][1]) for i,x in enumerate(alleles_found)]
		
		phase_bias=(min(phase_total)<max(phase_total)/3)
		
		same_phase_total=sum([x[0]+x[1] for x in barcode_phases])-reverse_phase_total
				
		#if reverse_phase_total<max(3*same_phase_total,self.confi_cut+same_phase_total) or alleles_phase.count(-1)!=1 or alleles_phase.count(1)!=1 or mean_quals[alleles_found[1]]<=self.qualcuts["ATCGDNMI"[min(7,alleles_found[1])]] or mean_quals[alleles_found[0]]<=self.qualcuts["ATCGDNMI"[min(7,alleles_found[0])]] or sum([alleles_counts[x] for x in alleles_found])<alleles_sum*0.8: 
		if reverse_phase_total<max(2*same_phase_total,self.confi_cut+same_phase_total) or alleles_phase.count(-1)!=1 or alleles_phase.count(1)!=1 or phase_bias or sum([alleles_counts[x] for x in alleles_found])<alleles_sum*0.70: 
			
			if coordinate in parents and 4 not in alleles_found and parents[coordinate][1]!=parents[coordinate][2]:
				
				print "nega",coordinate,max([alleles_counts[x] for x in alleles_found]),min([alleles_counts[x] for x in alleles_found]), parents[coordinate],max_label,barcode_phases,alleles_phase,allele_accuracy
		
				
				print "nega",reverse_phase_total, same_phase_total,quals
				
				gcounter[1]+=1
					
				print gcounter
				
			return 0,[]
		
		if coordinate in parents and 4 not in alleles_found:
		
			if parents[coordinate][1]==parents[coordinate][2]:
				print coordinate,max([alleles_counts[x] for x in alleles_found]),min([alleles_counts[x] for x in alleles_found]), parents[coordinate],max_label,barcode_phases,alleles_phase,allele_accuracy
		
				
				print reverse_phase_total, same_phase_total,quals
				gcounter[0]+=1
					
			else:
				gcounter[2]+=1
					
		
				print gcounter
		
		return 0,[]
		if alleles_phase.count(1)>1 or alleles_phase.count(-1)>1 or (4 in alleles_found and max([-error_bionom.cdf(sum(barcode_phase), min(barcode_phase)) for barcode_phase in barcode_phases if min(barcode_phase)<0.8*sum(barcode_phase)]+[0])>5) or (max([-error_bionom.cdf(sum(barcode_phase), min(barcode_phase)) for barcode_phase in barcode_phases if min(barcode_phase)<0.8*sum(barcode_phase)]+[0])>7):
			
			#print coordinate, barcode_phases, [-bionom.cdf(sum(barcode_phase), min(barcode_phase)) for barcode_phase in barcode_phases]

			return 0,[]
		
		
		alleles_max=max(alleles_counts)
		confi_alleles=[x for x,phase in enumerate(alleles_phase) if phase and alleles_counts[x]>alleles_max/2]
		
		valide_alleles=[]
		
		if len(confi_alleles)==2 and len(alleles_found)==2:
			
			valide_alleles=sorted(list(confi_alleles),key=lambda x:alleles_counts[x])
						
			score=self.score(abs(coordinate),confi_alleles,confi_alleles, barcodes,barcode_phases,del_sizes)	
						
			return 2,[coordinate,score,max_label,zip(*[[alleles[index] for index in valide_alleles],[barcodes[index] for index in valide_alleles],[quals[index] for index in valide_alleles],[barcode_phases[index] for index in valide_alleles],[alleles_phase[allele] for allele in valide_alleles],valide_alleles])]
	
		phased_alleles=[x for x,phase in enumerate(alleles_phase) if phase]
				
		
		valide_alleles=self.valide_alleles(coordinate,barcodes, alleles_found,alleles_cutoff, mean_quals)
						
		valide_alleles=set(valide_alleles+phased_alleles)
		
		if len([x for x in del_sizes if x>10])>min(sum(alleles_counts)/4,self.coverage.value/4):
			
			valide_alleles.add(4)
			
			self.last_coordi=abs(coordinate)
		
		valide_alleles=sorted(list(valide_alleles),key=lambda x:alleles_counts[x])
		
		if len(valide_alleles)<=1:
						
			return 0,[]
			
		score=self.score(abs(coordinate),valide_alleles, phased_alleles,barcodes,barcode_phases,del_sizes)	
		
		return -len(valide_alleles),[coordinate,score,max_label,zip(*[[alleles[index] for index in valide_alleles],[barcodes[index] for index in valide_alleles],[quals[index] for index in valide_alleles],[barcode_phases[index] for index in valide_alleles],[alleles_phase[allele] for allele in valide_alleles],valide_alleles])]
		
	def tovcf(self, line):	
		
		global refdata,mole_data
		
		number_alleles, out=self.process_data(line)
		
		if not number_alleles or abs(number_alleles)>2:
			return ""
		
		elif not self.ifvcf.value:
			if abs(number_alleles)>1:
				return line
			else:
				return ""
		
		coordinate, score,max_label,data=out
			
		ref_base=refdata[abs(coordinate)-1].upper()
				
		if number_alleles<0 or (abs(number_alleles)>1 and [1 for i,allele_data in enumerate(data) if max(allele_data[3])<max(self.confi_cut if i !=4 else self.indel_cut, len(allele_data[1])/4)]):
			Filter="LowQual"
			ifphase='/'
		else:
			Filter="PASS"
			ifphase='|'
		
		an=len([x for x in data if x[0]!=ref_base])
		
		if abs(number_alleles)==1:
			ac=2
		else:
			ac=1
			
		ref_index=[index for index,allele_data in enumerate(data) if allele_data[0]==ref_base]
		posi_var_index=[index for index,allele_data in enumerate(data) if allele_data[0]!=ref_base and allele_data[4]>=0]
		nega_var_index=[index for index,allele_data in enumerate(data) if allele_data[0]!=ref_base and allele_data[4]<0]
		min_len=min(len(posi_var_index),len(nega_var_index))
		
		sort_index=[x for pair in zip(posi_var_index,nega_var_index) for x in pair]+posi_var_index[min_len:]+nega_var_index[min_len:]
		data=[data[index] for index in ref_index+sort_index]
		
		if ref_index:
			ref_index=[0]

		ref_qual=[]
		ro=0
		rq=0.00
		ref_phase=0
		if ref_index:
			ref_phase=data[ref_index[0]][4]	
			ref_qual=data[ref_index[0]][1]		
			ro=len(ref_qual)
			if ref_qual:
				rq=str(round(sum(ref_qual)/len(ref_qual),2))
		
		dp=str(sum(map(lambda x:len(x[1]),data)))
		pq=str(round(20.00*sum(map(lambda x:abs(x[3][1]-x[3][0]),data)),2))
		
		infors,fields,alts,phases=[],[],[],[]
		
		infor_proto=["AN=%d"%an, "AC=%d"%ac, "AF=%s"%str(round(1.0/(an+len(ref_index)),2)), "BaseQRankSum=0.00"]
		
		if max_label:
			field_proto=["0", dp, "%d"%ro, "%s"%rq, "0", "%s"%str(0.0), "%d"%abs(max_label), "%s"%pq, ""]	
		else:
			field_proto=["0", dp, "%d"%ro, "%s"%rq, "0", "%s"%str(0.0),""]		
		
		for index,allele_data in enumerate(data):
			
			infor=[x for x in infor_proto]
			field=[x for x in field_proto]
			
			allele,barcodes,quals,barcode_phases,alleles_phase,alleles_index=allele_data
						
			if ref_index and len(ref_qual)>1 and len(quals)>1:
				
				BaseQRankSum=str(round(-10*math.log10(ranksums(ref_qual,quals)[1]+1.0e-10),2))
				
				infor[-1]="BaseQRankSum=%s"%BaseQRankSum
			
			ao=len(quals)
			aq=sum(quals)/ao if ao else 0.00
			
			field[0]="%d"%(min(2,index+1-len(ref_index)))
			field[-1]=",".join(["%d_%d"%(abs(bar),qual) for bar, qual in zip(barcodes, quals)])
			field[4]="%s"%str(ao)
			field[5]="%d"%abs(max_label)
			
			infors.append(infor)
			fields.append(field)
			alts.append(allele)
			phases.append(alleles_phase)
					
		pair_infor=[";".join(map(lambda x:x[0]+','+x[1].split('=')[1],zip(infors[2*index],infors[2*index+1]))) for index in xrange((len(infors))/2)]
		
		pair_field=[allele_pair(fields[2*index],fields[2*index+1],ifphase, [phases[2*index],phases[2*index+1]]) for index in xrange((len(infors))/2)]
		
		pair_alt=[alts[2*index]+','+alts[2*index+1] for index in xrange((len(infors))/2)]
		
		if len(infors)%2:
			
			fields[-1][0]=fields[-1][0]+'/'+fields[-1][0]

			pair_infor.append(";".join(infors[-1]))
			pair_field.append(":".join(fields[-1]))
			pair_alt.append(alts[-1])
		
					
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BC04103_longranger
		vcfline="\n".join(["\t".join([self.chrom.value, str(abs(coordinate)), ".", ref_base, ",".join(alts[len(ref_index):]), str(score), Filter, infor, "GT:DP:RO:QR:AO:QA:PS:PQ:BX", field]) for infor, field, alt in zip(pair_infor, pair_field, pair_alt)] )
				
		return vcfline		
		



	
def proxy(lines,counter,results):
	
	global processor, bionom
	
	for counter,line in lines:
				
		results[counter%10000]=processor.tovcf(line)
	
class block:
	
	def __init__(self,numthreads):
		self.results=manager.list(["" for x in xrange(10000)])
		self.threads_input=[[] for x in xrange(numthreads)]
		self.processes=[]
		self.counter=0
		self.numthreads=numthreads
				
	
	def __del__(self):
		
		del self.results[:]
		
		for x in self.threads_input:
			del x[:]
		del self.processes[:]
		
		del self.results,self.threads_input, self.processes, self.pool
		

	def takeline(self,line):
		
		self.threads_input[self.counter%self.numthreads].append((self.counter%10000,line))
		self.counter+=1
	
	def run(self,fwrite):
		
		self.pool=mul.Pool(processes=self.numthreads)

		
		for thread in xrange(self.numthreads):
						
			self.pool.apply_async(proxy,args=(self.threads_input[thread], self.counter,self.results))
			
			
		self.pool.close()
		self.pool.join()
			
		fwrite.write("\n".join([x for x in self.results if len(x)])+'\n')
				
	
		

def run(args):
	
	inputfile=args.input
	outputfile=args.output
	infofile=args.info
	reffile=args.ref
	molefile=args.mole
	snpfile=args.snp
	sample=args.sample
	ifvcf=args.ifvcf
	chrom=args.chr
	numthreads=args.threads
	confi=args.confidence
	
	
	global confi_cut, processor,refdata,mole_data, error_bionom, allele_bionom
	
	confi_cut=20/(confi)
	
	processor=line_processor() 
	error_bionom=bionomial(0.05)
	allele_bionom=bionomial(0.5)
	
	#snp_data=pd.read_csv(snpfile,header=None,sep='\t',dtype=int).set_index(0)[1].to_dict()
		
	if ifvcf:
		fwrite=open(inputfile+'_valid.vcf',mode='w')
	else:	
		fwrite=open(inputfile+'_valide',mode='w')
	
	counter=0
	out=[]
	fread=open(inputfile,mode='r')	
	
	processor.load_phasedata(ifvcf, chrom, sample, infofile, reffile, molefile)
	
	mole_data=read_mole(molefile,snpfile)
	
	refdata=read_ref(reffile, chrom)
	
	
	eachblock=block(numthreads)
	
	for counter,line in enumerate(fread):
		
		
		if 1==0:	
			if abs(int(line.split()[0])) not in [874766]:
				continue
			print line
			print processor.process_data(line)
		
		processor.process_data(line)
		
		continue
		eachblock.takeline(line)
		
		if not (counter+1)%10000:
			
			eachblock.run(fwrite)
			
			del eachblock
			
			gc.collect()

			eachblock=block(numthreads)
			
			print counter+1
			
	eachblock.run(fwrite)
						
	fwrite.close()
	fread.close()
	


def main():
		
	parser=argparse.ArgumentParser(description="algthorithm to validate all snps using phased data")
	parser.add_argument("-i","--input",help="input file" ,dest="input",type=str,required=True)
	parser.add_argument("-f","--info",help="info file" ,dest="info",type=str,required=True)
	parser.add_argument("-o","--output",help="output file" ,dest="output",type=str,required=True)
	parser.add_argument("-r","--ref",help="reference" ,dest="ref",type=str,required=True)
	parser.add_argument("-c","--chr",help="chromosome" ,dest="chr",type=str,required=True)
	parser.add_argument("-m","--mole",help="referenced phased molecules" ,dest="mole",type=str,required=True)
	parser.add_argument("-n","--snp",help="referenced phased snps" ,dest="snp",type=str,required=True)
	parser.add_argument("-s","--sample",help="sample name" ,dest="sample",type=str,default="sample")
	parser.add_argument("-v","--outvcf",help="if output as vcf" ,dest="ifvcf",type=int,default=1)
	parser.add_argument("-t","--threads",help="number of threads" ,dest="threads",type=int,default=1)
	parser.add_argument("-d","--confidence",help="phased confidence" ,dest="confidence",type=int,default=4)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)	

	
if __name__=='__main__':
	
	main()