#!/usr/bin/python

import re 
import pandas as pd
import os
import argparse
import collections as cl
import numpy as np
import random as rd	

global splitcomma
splitcomma=re.compile(r"(?<!,),")

global allbarcodes
allbarcodes=set([])

def group_inserts(elements):
	
	if not len(elements):
		
		return {}
	
	inserts=cl.defaultdict(list)
	
	for insert in re.split(splitcomma,","+elements)[2:-1]:
		
		insert=insert[2:]
				
		barcode=insert.split(":")[0]
		
		genotype=insert.split(":")[1]
				
		if len(genotype)>3:
			genotype=genotype[:2]+genotype[-1]
		
		if int(barcode) and len(genotype):
			
			inserts[genotype].append(abs(int(barcode)))
	
	return inserts


def process_line(line,refgenotypes):
	
	splits=line.split("\t")
	
	coordinate=int(splits[0])
		
	refgenotype=refgenotypes.get(coordinate,[])
	
	columns=[map(lambda x: abs(int(x.split(":")[0])),[a[2:] for a in re.split(splitcomma,","+x)[1:-1]]) for i,x in enumerate(splits[1:]) if i!=5 and i<7]
	
	inserts={}	
	if len(splits)>8:
		inserts=group_inserts(splits[8])
		
	
	alleles=['A','T','C','G','D','M']+inserts.keys()
	
	columns.extend(inserts.values())
			
	barcodes=cl.defaultdict(lambda: [0 for i in xrange(len(alleles))])
	
	col_notused=[]
	
	for i,column in enumerate(columns):
		
		for barcode in column:
			
			barcodes[barcode][i]+=1;
			
	barcodes={(bar if max(allele_counts)>=0.9*sum(allele_counts) else col_notused.append(bar)):allele_counts.index(max(allele_counts)) for bar,allele_counts in barcodes.items()}
	
	if None in barcodes:
		del barcodes[None]
		
	allbarcodes.update(barcodes.keys())	
	
	columns=[[] for i in xrange(len(alleles))]
	
	for barcode,allele in barcodes.items():
		
		columns[allele].append(barcode)
	
	if refgenotype:
		
		row_sortindex=[i for i,allele in enumerate(alleles) if allele in refgenotype]
		
		if len(row_sortindex)!=2:
			return -1,[],[]
			
	elif splits[-1][-2]=='-' or splits[-1][-1]=='-':
				
		refgenotype=[x.upper()  for x in splits[-1].split('-')[:-1]]
		
		if 'I' in refgenotype and not len(inserts):
			
			return 0, [], []
				
		row_sortindex=[i for i,allele in enumerate('ATCGDMI') if allele in refgenotype]
		
		if len(row_sortindex)<2:
			
			col_counts=map(len,columns)
			
			row_sortindex+=sorted([x for x in range(len(alleles)) if x not in row_sortindex], key=lambda i:col_counts[i],reverse=True)[:1]
			
	else:
		
		col_counts=map(len,columns)
		
		sv_dele=len([x for x in re.split(splitcomma,","+splits[5])[1:-1] if sum(map(int, x[1:].split(":")[1].split('-')))>10])
		
		if sv_dele>sum(col_counts)/4:
						
			row_sortindex=sorted(range(len(alleles)), key=lambda i:col_counts[i],reverse=True)[:2]
			
		else:
			
			row_sortindex=sorted(range(len(alleles)), key=lambda i:col_counts[i],reverse=True)[:2]
			
			if rd.getrandbits(1):
				
				row_sortindex=row_sortindex[::-1]
				
		
		if 4 in row_sortindex:
			
			row_sortindex=[4]+[x for x in row_sortindex if x!=4]
		
			
	return coordinate,[columns[i] for i in row_sortindex], [alleles[i] for i in row_sortindex]


		
			

def run(args):
	
	inputfile=args.inputfile
	
	outputfile=args.outputfile
	
	refvcf=args.refvcf
			
	refgenotypes={}
	last_coordi=0
	if len(refvcf):
		
		with open(refvcf, mode="r") as fvcf:
			
			for line in fvcf:
				
				line=line.split()
				
				coordinate=int(line[1])
				
				if last_coordi==coordinate:
					continue
					
				last_coordi=coordinate
								
				current_allele=[line[3]]+line[4].split(',')

				genotype=line[-1].split(":")[0]
				genotype=[current_allele[int(genotype[0])],current_allele[int(genotype[-1])]]
					
				refgenotypes[coordinate]=genotype
	
		fvcf.close()
	

	fdata=open(outputfile,mode='w')
	readfile=open(inputfile,mode='r')
	
	line_index=0
	for line in readfile:
		
		if len(line.strip())<2:
			continue
		
		line_index+=1
		
		if not line_index%1000:
			print line_index
		
		coordinate,data,alleles=process_line(line,refgenotypes)
		
		if not len(alleles):
			continue
				
		fdata.write("%d\t%s,\t%s,\t%s,%s\n"%(abs(coordinate),','.join(map(str,data[0])),','.join(map(str,data[1])),alleles[0],alleles[1]))
		
	fdata.close()
	readfile.close()
	
	totalbar=len(allbarcodes)
	
	
	
	

def main():
	parser=argparse.ArgumentParser(description="script to convert all snp data to binary")
	parser.add_argument("-i","--input",help="input snp data" ,dest="inputfile",type=str)
	parser.add_argument("-o","--output",help="output file" ,dest="outputfile",type=str)
	parser.add_argument("-r","--refvcf",help="reference allele" ,dest="refvcf",default="",type=str)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)	
	

if __name__=='__main__':
	
	main()