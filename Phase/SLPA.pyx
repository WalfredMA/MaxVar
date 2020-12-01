#!/usr/bin/python
# coding=utf-8
#distutils: language = c++

import cython
import numpy as np 
cimport cython 
cimport numpy as np
from libc.string cimport memset
from libc.stdlib cimport malloc as Malloc, free as Free, realloc as Realloc
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

DTYPE=np.int32

ctypedef np.int32_t DTYPE_t	


global debug 
debug=0

global allcoordi
allcoordi=[]

cdef struct label:
	
	int index
	int total
	int posi
	int nega
	char sign
	int value

cdef struct mole:
	
	label* labels
	unsigned short num_labels
	
	int* pconnected
	unsigned short num_pconnected
		
	int* nconnected
	unsigned short num_nconnected


cdef struct snp:
	
	label* labels
	unsigned short num_labels
	unsigned short size_labels
		
	int detered_label
	char detered_sign
	
	int* pconnected
	unsigned short num_pconnected
	
	int* nconnected
	unsigned short num_nconnected
	
	int init_label
	char init_sign
	int coordi


cdef struct blk:
	
	label* labels
	unsigned short num_labels
	unsigned short size_labels
	
	int init_label
	
	int detered_label
	char detered_sign
	unsigned short detered_value
	
	int* connected
	char* connected_sign
	int num_connected	
	int size_connected
	


cdef inline label *push_back_connect(blk* current_blk, int new_arr_int, char new_arr_char):

	cdef int arr_l=current_blk.num_connected
	cdef int arr_size=current_blk.size_connected
	

	if arr_l ==arr_size or arr_size-arr_l>100:
		
		current_blk.connected=<int*> Realloc(current_blk.connected, (arr_l+100)*sizeof(int))
		current_blk.connected_sign=<char*> Realloc(current_blk.connected_sign, (arr_l+100)*sizeof(char))
		current_blk.size_connected=arr_l+100

	current_blk.connected[arr_l]=new_arr_int
	current_blk.connected_sign[arr_l]=new_arr_char
	current_blk.num_connected+=1

cdef inline label *push_back_snp(snp *current_snp):

	cdef int arr_l=current_snp.num_labels,arr_size=current_snp.size_labels

	if arr_l==arr_size or arr_size-arr_l>100:
		current_snp.labels=<label*> Realloc(current_snp.labels,(arr_l+100)*sizeof(label))
		current_snp.size_labels=arr_l+100
	
	current_snp.num_labels+=1
	
	return &current_snp.labels[arr_l]

cdef inline label *push_back_blk(blk *current_blk):

	cdef int arr_l=current_blk.num_labels,arr_size=current_blk.size_labels

	if arr_l==arr_size or arr_size-arr_l>100:
		
		current_blk.labels=<label*> Realloc(current_blk.labels, (arr_l+100)*sizeof(label))
		current_blk.size_labels=arr_l+100
	
	current_blk.num_labels+=1
	
	return &current_blk.labels[arr_l]

cdef inline void tops(np.int32_t[:] tops_values, np.int64_t[:] tops_indexs, int new_value, int new_index, int top_num):

	for i in xrange(top_num):
		
		if abs(new_value)>abs(tops_values[i]):
			
			for j in xrange(i+1,top_num):
			
				tops_values[j]=tops_values[j-1]
				tops_indexs[j]=tops_indexs[j-1]
				
			tops_values[i]=new_value
			tops_indexs[i]=new_index
			
			break

cdef inline int getlabelsign(snp* current_snp, int findlabel):

	for index in xrange(current_snp.num_labels):
				
		if current_snp.labels[index].index==findlabel:
		
			return current_snp.labels[index].sign
	
	return 0
	
		



cdef allvalidelabel_num(np.int32_t[:] sorted_coordi, set allLabels, int cutoff):

	cdef np.ndarray[np.int32_t, ndim=1,mode="c"] result = np.zeros([len(sorted_coordi)], dtype = np.int32)
		
	cdef np.ndarray[np.int32_t, ndim=1,mode="c"] allvalide_snp = np.zeros([len(sorted_coordi)], dtype = np.int32)
	cdef np.ndarray[np.int32_t, ndim=1,mode="c"] allvalide_snpindex = np.zeros([len(sorted_coordi)], dtype = np.int32)
	
	cdef int coordi, num_valide_blk=0, num_valide_snp=0
	
	for i in xrange(len(sorted_coordi)):
	
		coordi=sorted_coordi[i]
		
		for i in xrange(num_valide_snp):	
			
			if coordi-allvalide_snp[i]>cutoff:
				num_valide_snp-=1
				if allvalide_snpindex[i] in allLabels:
					num_valide_blk-=1
			else:
				break
		
		result[i]=num_valide_blk
		num_valide_snp+=1
		allvalide_snp[num_valide_snp]=coordi
		allvalide_snpindex[num_valide_snp]=i
		
		if coordi in allLabels:
		
			for i in xrange(num_valide_snp):	
		
				result[allvalide_snpindex[num_valide_snp]]+=1
		
			num_valide_blk+=1
					
	
	return result
	
#find the max abs value index
cdef filter_sort(dict associations, np.int64_t[:] associations_sortindex,np.int32_t[:] associations_sortvalue,  int cutoff):

	cdef int size=len(associations)
	
	cdef np.ndarray[np.int32_t, ndim=1,mode="c"] valide_values = np.empty([size], dtype=np.int32)
	cdef np.ndarray[np.int32_t, ndim=1,mode="c"] valide_indexs = np.empty([size], dtype=np.int32)
	
	cdef int newiter
	cdef int num_pass=0
	cdef int index,sign, val
	for i,newiter in associations.items():
						
		if abs(newiter)>=cutoff:
			
			if newiter>0:
				sign=1
			else:
				sign=-1
			
			valide_values[num_pass]=abs(newiter)
			valide_indexs[num_pass]=i*sign
			num_pass+=1	
	
	cdef np.ndarray[np.int64_t, ndim=1] sortindex=valide_values[:num_pass].argsort()[::-1]
	
	 
	for i in xrange(num_pass):
	
		index=valide_indexs[sortindex[i]]
		val=valide_values[sortindex[i]]
		sign=index/abs(index)
		index*=sign
		val*=sign
		
		associations_sortindex[i]=index
		associations_sortvalue[i]=val
		
	return num_pass


#pass the label from snp layer to molecule layer
cdef pass_label_tomole(snp *snps, mole *current_mole, dict current_mole_valide_labels, int sign, int snp_max):

	cdef int current_index, index, dict_index
	cdef int num_current_mole_labels=current_mole.num_labels
	cdef label* current_des_label
	cdef label* current_des_labels=current_mole.labels
	cdef char current_sign
	cdef int current_val=1
	
	if sign>0:
		connected=current_mole.pconnected
		connected_size=current_mole.num_pconnected
		
	elif sign<0:
		connected=current_mole.nconnected
		connected_size=current_mole.num_nconnected
	
	for index in xrange(connected_size):
									
		current_snp=snps[connected[index]]
		
		current_index=current_snp.detered_label
		
		if not current_index:
			continue
				
		current_sign=sign*current_snp.detered_sign
		
		dict_index=current_mole_valide_labels.get(current_index,-1)
				
		if dict_index<0:
		
			current_mole_valide_labels[current_index]=current_mole.num_labels
			
			current_des_label=&current_des_labels[current_mole.num_labels]
			
			current_mole.num_labels+=1
			
			current_des_label.index=current_index
					
			current_des_label.total=current_val
			
			if current_sign>0:
			
				current_des_label.posi=current_val
				
				current_des_label.nega=0
			
			elif current_sign<0:
			
				current_des_label.posi=0
				
				current_des_label.nega=current_val
				
			else:
				
				current_des_label.posi=0
				
				current_des_label.nega=0
												
		else:
		
			current_des_label=&current_des_labels[dict_index]
			
			current_des_label.total+=current_val
		
			if current_sign>0:
			
				current_des_label.posi+=current_val

			elif current_sign<0:
				
				current_des_label.nega+=current_val
				
	


#analyze and phase each molecule
cdef summary_mole(mole *current_mole, dict current_mole_valide_labels, int snp_cutoff, float phase_cutoff):

	cdef int current_index,current_sign,total_count, nega_count, posi_count
	
	cdef label* current_des_labels=current_mole.labels
	cdef label* current_des_label
	cdef float the_cutoff
	
	#phase all labels it has
	for current_index in xrange(current_mole.num_labels):
		
		current_des_label=&current_des_labels[current_index]
		
		total_count, posi_count, nega_count=current_des_label.total, current_des_label.posi, current_des_label.nega
		
		the_cutoff=max(phase_cutoff*(posi_count+nega_count),total_count/2+1,snp_cutoff)
		
		if total_count<the_cutoff:
			current_des_label.sign=0
			current_des_label.value=1
			continue
		
		elif nega_count>=the_cutoff:
		
			true_count=nega_count-posi_count
			current_sign=-1
		
		elif posi_count>=the_cutoff:
			
			true_count=posi_count-nega_count
			current_sign=1
		
		else:
			true_count=total_count
			current_sign=0
			
		current_des_label.sign=current_sign	
		current_des_label.value=true_count
	
	

#pass the label from molecule layer to snp layer
cdef pass_label_tosnp(snp *current_snp, dict current_snp_valide_labels, mole *moles, int* label_coordi, int sign, unsigned short mole_max, int max_dis):
	
	cdef mole* current_mole
	cdef label* current_ori_labels
	cdef label* current_ori_label
	
	cdef label *current_des_labels=current_snp.labels
	cdef label *current_des_label
	
	cdef int* connected
	cdef unsigned short connected_size
	if sign>0:
		connected=current_snp.pconnected
		connected_size=current_snp.num_pconnected
		
	elif sign<0:
		connected=current_snp.nconnected
		connected_size=current_snp.num_nconnected
	
	cdef int index, dict_index, label_index, current_index, mole_index1
	cdef char current_sign
	cdef int current_val
	cdef int current_coordi=current_snp.coordi
	
	#iter all connected molecules 
	for index in xrange(connected_size):
	
		mole_index1=connected[index]
	
		current_mole=&moles[mole_index1]
															
		if not current_mole.num_labels:
				
			continue
	
		current_ori_labels=current_mole.labels
			
		 
		for current_index in xrange(current_mole.num_labels):
		
			current_ori_label=&current_ori_labels[current_index]
			
			label_index=current_ori_label.index
			
			current_val=current_ori_label.value
						
			if not current_val or abs(label_coordi[label_index]-current_coordi)>max_dis:
				continue
			
			current_val=min(mole_max,current_val)
			
			current_sign=sign*current_ori_label.sign
			
			dict_index=current_snp_valide_labels.get(label_index,-1)
			
			if dict_index<0:
											
				current_snp_valide_labels[label_index]=current_snp.num_labels
				
				current_des_label=push_back_snp(current_snp)
				
				current_des_labels=current_snp.labels
				
				current_des_label.index=label_index

				current_des_label.total=current_val
				
				if current_sign>0:
				
					current_des_label.posi=current_val
					
					current_des_label.nega=0
				
				elif current_sign<0:
				
					current_des_label.posi=0
					
					current_des_label.nega=current_val
					
				else:
					
					current_des_label.posi=0
					
					current_des_label.nega=0
			
														
			else:
				
				current_des_label=&current_des_labels[dict_index]
				
				current_des_label.total=current_des_label.total+current_val			
				
				if current_sign>0:
				
					current_des_label.posi=current_des_label.posi+current_val

				elif current_sign<0:
					
					current_des_label.nega=current_des_label.nega+current_val
	

#analyze and phase each snp
cdef summary_snp(snp *current_snp, dict current_snp_valide_labels,int mole_cutoff, float phase_cutoff):
	
	cdef int index, current_index, current_sign, total_count, nega_count, posi_count, true_count
	
	cdef label *current_des_labels=current_snp.labels
	cdef label *current_des_label
	cdef float the_cutoff

	#filter and analyze all labels this snp has 
	for current_index in xrange(current_snp.num_labels):
		
		current_des_label=&current_des_labels[current_index]
		
		total_count, posi_count, nega_count=current_des_label.total, current_des_label.posi, current_des_label.nega
		
		#print current_snp.coordi, current_des_label.index, total_count, posi_count, nega_count
		
		the_cutoff=max(phase_cutoff*(posi_count+nega_count),total_count/2+1,mole_cutoff)
		
		if total_count<the_cutoff:
			current_des_label.value=0
			current_des_label.sign=0
			continue
		
		elif nega_count>=the_cutoff:
		
			current_sign=-1
			true_count=nega_count-posi_count
		
		elif posi_count>=the_cutoff:
			
			current_sign=1
			true_count=posi_count-nega_count
		
		else:
			
			true_count=total_count
			current_sign=0
				
					
		current_des_label.value=true_count
		current_des_label.sign=	current_sign
		
		
		
#pass the label from snp layer to blocks
cdef pass_label_toblk(snp *snps, blk *current_blk,  dict current_blk_valide_labels, int old_index, int snp_max):
		
	cdef int label_index, dict_index, num_current_ori_labels, current_val, current_index, snp_index1, temp_max, total_count, nega_count, posi_count
	
	cdef char current_sign
	
	cdef label	*current_ori_labels, *current_ori_label
	cdef label	*current_des_labels=current_blk.labels, *current_des_label
	
	cdef int* connected=current_blk.connected
	cdef char* connected_sign=current_blk.connected_sign
	cdef unsigned short connected_size=current_blk.num_connected
			
	cdef snp* current_snp

	for index in xrange(connected_size):
	
		snp_index1=connected[index]
		
		current_sign=connected_sign[index]
		
		#not use the snp already nulled 
		if not current_sign:
			continue
		
		current_snp=&snps[snp_index1]
		
		current_ori_labels=current_snp.labels
		
		num_current_ori_labels=current_snp.num_labels
		
		#storage all labels this snp has 
		for current_index in xrange(num_current_ori_labels):
						
			current_ori_label=&current_ori_labels[current_index]
			
			label_index=current_ori_label.index
			
			if label_index==old_index:

				continue	
				
			total_count=current_ori_label.total
			
			#snp can be positively connect or negatively connected
			if current_sign>0:
				posi_count, nega_count=current_ori_label.posi, current_ori_label.nega
			
			elif current_sign<0:	
				posi_count, nega_count=current_ori_label.nega, current_ori_label.posi
			
			#print "blk",old_index,label_index, total_count, posi_count, nega_count
			
			#if this label already passed by other snps
			dict_index=current_blk_valide_labels.get(label_index,-1)
					
			if dict_index<0:
						
				current_blk_valide_labels[label_index]=current_blk.num_labels
								
				current_des_label=push_back_blk(current_blk)
				
				current_des_labels=current_blk.labels
								
				current_des_label.index=label_index
				
				#if over the max value of snp propagation, then shrink by proportion
				if total_count>snp_max:
				
					current_des_label.posi=posi_count*snp_max/total_count
				
					current_des_label.nega=nega_count*snp_max/total_count
														
					current_des_label.total=snp_max
				
				else:
				
					current_des_label.posi=posi_count
				
					current_des_label.nega=nega_count
													
					current_des_label.total=total_count
					
								
			else:
							
				current_des_label=&current_des_labels[dict_index]
			
			
				if total_count>snp_max:
				
					current_des_label.posi+=posi_count*snp_max/total_count
				
					current_des_label.nega+=nega_count*snp_max/total_count
														
					current_des_label.total+=snp_max
				
				else:
				
					current_des_label.posi+=posi_count
				
					current_des_label.nega+=nega_count
													
					current_des_label.total+=total_count
	
	
	

#analyze and phase each block
cdef summary_blk(blk *current_blk, dict current_blk_valide_labels, dict associations,dict associations_full,  int blk_index, int lblk, int blk_cutoff, float phase_cutoff):

	cdef int  index, current_val, current_sign, total_count, posi_count, nega_count, true_count
	
	cdef long long current_index
	cdef label* current_label
	cdef float the_cutoff
	
	cdef int top_num=10	 #only use 10 highest labels 

	cdef np.ndarray[np.int32_t, ndim=1] tops_values =  np.zeros(top_num, dtype=np.int32)
	cdef np.ndarray[np.int64_t, ndim=1] tops_indexs =  np.zeros(top_num, dtype=np.int64)
		
	#filter and analyze all labels it has 	
	for index in xrange(current_blk.num_labels):
		
		current_label=&current_blk.labels[index]
		
		current_index=current_label.index
		
		total_count, posi_count, nega_count=current_label.total, current_label.posi, current_label.nega
		
		#print current_blk.detered_label, current_label.index, total_count, posi_count, nega_count
						
		the_cutoff=max(phase_cutoff*(posi_count+nega_count),total_count/2+1, blk_cutoff)
		
		if nega_count>=the_cutoff:
			
			true_count=nega_count-posi_count
		
			current_sign=-1
		
		elif posi_count>=the_cutoff:
			
			true_count=posi_count-nega_count
			current_sign=1
		else:
			
			current_label.value=0
			current_label.sign=0
			continue
		
		current_label.value=true_count
		current_label.sign=current_sign
		
		#print (blk_index,current_index,true_count,current_sign,the_cutoff, total_count,posi_count,nega_count)
		
		if current_index<blk_index:
			
			current_index=lblk*current_index+blk_index
		
		else:
			current_index=lblk*blk_index+current_index
			
		tops(tops_values, tops_indexs, current_sign*true_count, current_index, top_num)
		
		#print "check0", current_index, true_count,current_sign
		
		if associations_full.get(current_index,0):
			associations_full[current_index]+=true_count*current_sign
			
		else:
			associations_full[current_index]=true_count*current_sign
				
	#only top labels are considered  
	for index in xrange(top_num):
	
		current_val=tops_values[index]
		
		if not current_val:
			
			continue
		
		current_index=tops_indexs[index]
		
		if associations.get(current_index,0):
			associations[current_index]+=current_val
		else:
			associations[current_index]=current_val
			



cdef class SLPA_network:
	
	cdef int snp_cutoff, mole_cutoff, snp_max, mole_max, blk_cutoff, max_dis
	cdef float phase_cutoff
	cdef int lsnp, lmole, lblk
	cdef int *popularity #storage all popularities of each block
	cdef int *label_coordi # all labels initiate coordinate, used to build sandboxs
	cdef blk* blks #data of all blocks
	cdef snp* snps #data of all snps
	cdef mole* moles #data of all molecules
	
	cdef dict associations #the association between labels used for switch
	cdef dict associations_full #full association between labels 
	cdef dict origin_labels #find label's original value (the label inputed and storaged in allLabels)
	
	
	def __cinit__(self):
		
		self.mole_max=3 #max confidence each molecule can provide
		self.mole_cutoff=3 #min mole confidence needed to phase
		self.snp_max=5 #max confidence each snp can provide
		self.snp_cutoff=1 #min snp confidence to phase
		self.blk_cutoff=5 #min value to phase block
		self.max_dis=1000000 #the max distance to propagate
		self.phase_cutoff=0.6
		
		self.origin_labels = {}	
		self.associations = {}
		self.associations_full={}
		self.lmole=0
		self.lsnp=0
		self.lblk=0
		
	def set_para(self, int blk_cutoff=5, int mole_max=3, int mole_cutoff=3, int snp_max=5, phase_cutoff=0.6):
		
		self.blk_cutoff=blk_cutoff
		self.mole_cutoff=mole_cutoff
		
		self.mole_max=max(mole_cutoff,mole_max)
		self.snp_max=max(snp_max,blk_cutoff)
		self.phase_cutoff=phase_cutoff
	
	def __dealloc__(self):
		
		for blk_index1 in xrange(1,self.lblk):
			
			Free(self.blks[blk_index1].labels)
			Free(self.blks[blk_index1].connected)
		
		for snp_index1 in xrange(1,self.lsnp):
			
			Free(self.snps[snp_index1].labels)
			Free(self.snps[snp_index1].pconnected)
			Free(self.snps[snp_index1].nconnected)
		
		for mole_index1 in xrange(1,self.lmole):
			
			Free(self.moles[mole_index1].labels)
			Free(self.moles[mole_index1].pconnected)
			Free(self.moles[mole_index1].nconnected)
		
		Free(self.blks)
		Free(self.snps)
		Free(self.moles)
		Free(self.popularity)
		
	#output snp data
	def output_snp(self):
		
		cdef int snp_index1
		
		cdef np.ndarray[DTYPE_t, ndim=1] snp_labels =  np.zeros(self.lsnp, dtype=DTYPE)
		cdef np.ndarray[DTYPE_t, ndim=1] snp_signs =  np.zeros(self.lsnp, dtype=DTYPE)
		
		for snp_index1 in xrange(1,self.lsnp):
			
			snp_labels[snp_index1]=self.snps[snp_index1].detered_label
			#snp_labels[snp_index1]=self.origin_labels[self.snps[snp_index1].detered_label]
			snp_signs[snp_index1]=self.snps[snp_index1].detered_sign
		
		return snp_labels, snp_signs
	
	def output_mole(self):
			
		cdef int mole_index1,label_index1,num_labels,sign
		cdef label* mole_labels, current_label
		
		allout_labels=[]
		out_labels=[]
		for mole_index1 in xrange(1,self.lmole):
			
			num_labels=self.moles[mole_index1].num_labels
			
			mole_labels=self.moles[mole_index1].labels
			
			out_labels=[]
			for label_index1 in xrange(0,num_labels):
				
				current_label=mole_labels[label_index1]
				
				if current_label.posi>current_label.total*0.6:
					
					out_labels.append(int(current_label.index))
				
				elif current_label.nega>current_label.total*0.6:
				
					out_labels.append(-int(current_label.index))
				
				else:
					continue
			
			allout_labels.append(out_labels)
		
		return allout_labels
	
	#phase all snps into blocks
	def phase(self):
	
		cdef int current_index, current_sign, blk_index1, snp_index1
		cdef blk* current_blk
		cdef snp* current_snp
		
		#reset all blocks
		for blk_index1 in xrange(1,self.lblk):
			
			self.blks[blk_index1].num_connected=0
		
		#phase each snp by adding their index into blocks
		for snp_index1 in xrange(1,self.lsnp):
			
			current_snp=&self.snps[snp_index1]
			#print current_snp.coordi,current_snp.detered_label, current_snp.detered_sign
			
			#low confidence snp, switch back to original phase
			if not current_snp.detered_sign:
			
				current_snp.detered_label=current_snp.init_label
				current_snp.detered_sign=current_snp.init_sign
			
			else:
				current_blk=&self.blks[current_snp.detered_label]	
				push_back_connect(current_blk, snp_index1, current_snp.detered_sign)
						
	#load molecular data and their connected snp
	def load_mole(self,np.ndarray[np.int32_t, ndim=1, mode="c"] starts_mole_p not None, np.ndarray[np.int32_t, ndim=1, mode="c"] matrix_mole_p not None, np.ndarray[np.int32_t, ndim=1, mode="c"] starts_mole_n not None, np.ndarray[np.int32_t, ndim=1, mode="c"] matrix_mole_n not None):
		
		cdef int start, end, snp_index1, mole_index1
		
		self.lmole=len(starts_mole_p)-1
		
		self.moles = <mole *> Malloc((self.lmole+1) * sizeof(mole))
		
		cdef mole *current_mole
		
		#initiate all molecule nodes 
		cdef int alloc_size
		for mole_index1 in xrange(1,self.lmole):
					
			current_mole=&self.moles[mole_index1]

			#storage positively connect snp index
			start, end= starts_mole_p[mole_index1:(mole_index1+2)]
			
			alloc_size=(end-start)
			current_mole.pconnected=<int *> Malloc((end-start+1) * sizeof(int))
			current_mole.num_pconnected=(end-start)
			if end>start:
				i=0
				for snp_index1 in np.nditer(matrix_mole_p[start:end]):
					current_mole.pconnected[i]=snp_index1
					i+=1
			
			#storage negatively connect snp index
			start, end= starts_mole_n[mole_index1:(mole_index1+2)]
			
			alloc_size+=(end-start)
			current_mole.nconnected=<int *> Malloc((end-start+1) * sizeof(int))
			current_mole.num_nconnected=(end-start)
			if end>start:
				i=0
				for snp_index1 in np.nditer(matrix_mole_n[start:end]):
					current_mole.nconnected[i]=snp_index1
					i+=1

			#the maximum label storaged=total connected snp
			current_mole.labels=<label *> Malloc((alloc_size+1) * sizeof(label))
			current_mole.num_labels=0	
			
	
	#load snp data and their connected molecules 
	def load_snp(self,np.ndarray[np.int32_t, ndim=1, mode="c"] starts_snp_p not None, np.ndarray[np.int32_t, ndim=1, mode="c"] matrix_snp_p not None, np.ndarray[np.int32_t, ndim=1, mode="c"] starts_snp_n not None,np.ndarray[np.int32_t, ndim=1, mode="c"] matrix_snp_n not None, np.ndarray[np.int32_t, ndim=1, mode="c"] snp_coordi not None):
	
		cdef int start, end, snp_index1, mole_index1i, i, j
		
		self.lsnp=len(starts_snp_p)-1
		self.snps = <snp *> Malloc((self.lsnp+1) * sizeof(snp))
		
		global allcoordi
		allcoordi=snp_coordi
		
		cdef snp* current_snp
				
		for snp_index1 in xrange(1,self.lsnp):
			
			#number of labels storaged is dynammic, use 100 as starting allocated size, but able to extend after
			current_snp=&self.snps[snp_index1]
			current_snp.coordi=snp_coordi[snp_index1]
						
			current_snp.num_labels=0
			current_snp.size_labels=100
			current_snp.labels=<label *> Malloc(100 * sizeof(label))
			
			#storage positively connect snp index
			start, end= starts_snp_p[snp_index1:(snp_index1+2)]
			current_snp.pconnected=<int *> Malloc((end-start+1) * sizeof(int))
			current_snp.num_pconnected=(end-start)
			if end>start:
				i=0
				for mole_index1 in np.nditer(matrix_snp_p[start:end]):
					current_snp.pconnected[i]=mole_index1
					i+=1
			
			#storage negatively connect snp index
			start, end= starts_snp_n[snp_index1:(snp_index1+2)]	
			current_snp.nconnected=<int *> Malloc((end-start+1) * sizeof(int))
			current_snp.num_nconnected=(end-start)		
			if end>start:
				i=0
				for mole_index1 in np.nditer(matrix_snp_n[start:end]):
					current_snp.nconnected[i]=mole_index1
					i+=1

	#load initiate block and labels on them
	def load_blk(self, np.ndarray[np.int32_t, ndim=1, mode="c"] init_blks not None, np.ndarray[np.int32_t, ndim=1, mode="c"] init_signs not None):
	
		cdef int start, end, snp_index1, mole_index1, blk_index1
		
		cdef list allLabels=list(set([0]+init_blks))
		
		self.lblk=len(allLabels)
		
		self.label_coordi=<int *> Malloc(self.lblk * sizeof(int))
		self.label_coordi[0]=0
				
		self.blks = <blk *> Malloc((self.lblk+1) * sizeof(blk))
		cdef blk *current_blk
		cdef mole *current_mole
		cdef snp *current_snp
						
		self.popularity =  <int *> Malloc(self.lblk * sizeof(int))
			
		#initiate all block nodes
		for blk_index1 in xrange(1,self.lblk):
		
			self.origin_labels[allLabels[blk_index1]]=blk_index1
			self.label_coordi[blk_index1]=self.snps[allLabels[blk_index1]].coordi
					
			current_blk=&self.blks[blk_index1]			
			current_blk.detered_label=blk_index1	
			current_blk.detered_sign=1
			
			current_blk.num_labels=0
			current_blk.size_labels=100
			current_blk.labels=<label *> Malloc(100 * sizeof(label))
			
			current_blk.num_connected=0
			current_blk.size_connected=100
			current_blk.connected=<int *> Malloc(100 * sizeof(int))
			current_blk.connected_sign=<char *> Malloc(100 * sizeof(char))
			
		#initiate initial labels for all snp nodes
		for snp_index1 in xrange(1,self.lsnp):
			
			current_snp=&self.snps[snp_index1]
			current_snp.init_label=init_blks[snp_index1]
			current_snp.init_sign=init_signs[snp_index1]
			current_snp.detered_label=init_blks[snp_index1]
			current_snp.detered_sign=init_signs[snp_index1]
				
		#initial phase 
		self.phase()
		
	#passing all labels to molecular_layer
	def molecular_layer(self):
		
		cdef int* connected
		cdef mole* current_mole
		cdef int mole_index1, snp_index1
		
		cdef dict current_valide_labels={}
		
		for mole_index1 in xrange(1,self.lmole):
		
			current_mole=&self.moles[mole_index1]
			
			current_mole.num_labels=0
								
			current_valide_labels.clear()
			
			#get labels from positively connected snp													
			pass_label_tomole(self.snps, current_mole, current_valide_labels, 1, self.snp_max)
			
			#get labels from negatively connected snp			
			pass_label_tomole(self.snps, current_mole, current_valide_labels, -1, self.snp_max)
						
			if 	current_mole.num_labels:
			
				#filter labels 
				summary_mole(current_mole, current_valide_labels, 1, self.phase_cutoff)
		
	
	#passing all labels to snp_layer
	def snp_layer(self):
		
		cdef int* connected		
		cdef snp* current_snp
		
		cdef int snp_index1,mole_index1
		
		cdef dict current_valide_labels={}
		
		for snp_index1 in xrange(1,self.lsnp):
					
			current_snp=&self.snps[snp_index1]
						
			current_snp.num_labels=0
				
			current_valide_labels.clear()	
			
			#get labels from positively connected molecules		
			pass_label_tosnp(current_snp, current_valide_labels,self.moles,self.label_coordi, 1, self.mole_max, self.max_dis)	
			
			#get labels from negatively connected molecules	
			pass_label_tosnp(current_snp, current_valide_labels, self.moles, self.label_coordi, -1,self.mole_max, self.max_dis)	
			
			if len(current_valide_labels):
			
				#filter labels 
				summary_snp(current_snp, current_valide_labels, self.mole_cutoff,self.phase_cutoff)	
			
	
	#passing all labels to block layer
	def blk_layer(self):
		
		cdef int* connected, valide_labels, valide_des_labels
		cdef int num_valide_label, num_des_labels
		
		cdef blk* current_blk
		
		cdef int index,label_index,old_label,current_index, current_sign
		cdef int total_count, posi_count, nega_count
		
		memset(self.popularity,0,sizeof(int)*self.lblk)
		
		cdef dict current_valide_labels={}
		
		for blk_index1 in xrange(1,self.lblk):
									
			current_blk=&self.blks[blk_index1]
			
			current_des_labels=current_blk.labels
									
			self.popularity[old_label]+=1
			
			current_blk.num_labels=0
			
			current_valide_labels.clear()
			
			#get labels from all connected snps				
			pass_label_toblk(self.snps, current_blk, current_valide_labels, blk_index1, self.snp_max)

			if len(current_valide_labels):
				
				#filter labels and determine all associations
				summary_blk(current_blk,current_valide_labels, self.associations, self.associations_full, blk_index1,self.lblk, self.snp_max,self.phase_cutoff)	
				
				
			
	#propagate and update labels of each block, also passing back to snp layer.
	def propagate(self):
			
		cdef int themax, index, index2, maxindex, firstblk, secondblk, switch_ori, switch_des, labelchangeto, current_sign, current_index,signchangeto, oldchangeblk, oldchangesign,pair_index, ifconflict
		
		cdef snp* current_snp
		cdef blk* switch_blk
		cdef label* current_label
		#find the connected block pair with the highest association value
			
		cdef list switch_ori_dic=[[] for i in xrange(self.lblk)] #original label to be switched to each label
		
		cdef np.ndarray[np.int64_t, ndim=1] switch_des_dic=np.zeros(self.lblk, dtype = np.int64) #each label's switch destination, one label can either be switch origin or switch destination, not both.
		
		#sort associations by their absolute values and storaged into sortindex and sortvalue by order
		cdef np.ndarray[np.int32_t, ndim=1] associations_sortvalue = np.zeros(len(self.associations), dtype = np.int32)
		cdef np.ndarray[np.int64_t, ndim=1] associations_sortindex = np.zeros(len(self.associations), dtype = np.int64)
		
		cdef int pass_num=filter_sort(self.associations, associations_sortindex, associations_sortvalue, self.blk_cutoff)
		
		#print self.associations
		#print [(x/self.lblk,x%self.lblk,associations_sortvalue[i]) for i,x in enumerate(associations_sortindex)]
			
		#print [(i/self.lblk,i%self.lblk,x) for i,x in self.associations_full.items()]
	
		#iter each association and determine the corrsponding switch 
		for index in xrange(pass_num):
		
			maxindex=associations_sortindex[index]
			themax=associations_sortvalue[index]
			
			if not themax:
				continue
			
			#two labels coupled
			firstblk=maxindex/self.lblk
			secondblk=maxindex%self.lblk	
			
			#print "maxindex",maxindex, self.lblk, firstblk, secondblk
			#if both already used as destinations cannot paired, too complex
			if (switch_des_dic[firstblk]) and (switch_des_dic[secondblk]):
			
				continue
						
			#determine which change to which when it is possible
			if self.popularity[firstblk]>=self.popularity[secondblk]:
			
				if not switch_des_dic[firstblk]:
					switch_des=firstblk
					switch_ori=secondblk
					
				else: 
					switch_des=secondblk
					switch_ori=firstblk
					
			else:
			
				if not switch_des_dic[secondblk]:
					switch_des=secondblk
					switch_ori=firstblk
					
				else:
					switch_des=firstblk
					switch_ori=secondblk
			
			
			#considering destination label might already switched, find the destination of destination
			if themax>0:
				signchangeto=1
			else:
				signchangeto=-1
				
			labelchangeto=switch_des
			while switch_des_dic[labelchangeto]:
				#redirect to final destination label
				labelchangeto=switch_des_dic[labelchangeto]

				if labelchangeto<0:
					signchangeto*=-1
					labelchangeto*=-1
					
			#check if any conflicts in associations, all have to be concordant
			ifconflict=0
			for oldchangeblk in switch_ori_dic[labelchangeto]:
			
				oldchangesign=1
				if oldchangeblk<0:
					oldchangesign=-1
					oldchangeblk*=-1
				
				if oldchangeblk<switch_ori:
				
					index2=self.lblk*oldchangeblk+switch_ori
				
				else:
				
					index2=self.lblk*switch_ori+oldchangeblk
				
				if self.associations_full.get(index2,0)*signchangeto*oldchangesign<0:
				
					#print "conflict", switch_ori, switch_des, labelchangeto, signchangeto,"\t",oldchangeblk, labelchangeto,oldchangesign,"\t",self.associations.get(index2,0),self.associations_full.get(index2,0)
					ifconflict=1
					break
			
			if ifconflict:
				continue
			
			#storage to switch
			switch_des_dic[switch_ori]=labelchangeto*signchangeto
			switch_ori_dic[labelchangeto].append(switch_ori*signchangeto)
			
			switch_blk=&self.blks[switch_ori]
			switch_blk.detered_label=labelchangeto
			switch_blk.detered_sign=signchangeto
		
		#swich all block's inner snp sites
		for switch_ori in xrange(self.lblk):
		
			if not switch_des_dic[switch_ori]:
				continue
					
			switch_blk=&self.blks[switch_ori]
			labelchangeto=switch_blk.detered_label
			signchangeto=switch_blk.detered_sign
						
			
			for index in xrange(switch_blk.num_connected):
				
				current_index=switch_blk.connected[index]
				current_sign=switch_blk.connected_sign[index]
				current_snp=&self.snps[current_index]
				
				current_snp.detered_label=labelchangeto				
				
				#if conflicts with snp itself's phase, then null this snp		
										
				if getlabelsign(current_snp, labelchangeto)!=-current_sign*signchangeto:	
					
					current_snp.detered_sign=current_sign*signchangeto					
					
				else:
				
					current_snp.detered_sign=0
					
					
		self.associations_full.clear()
		self.associations.clear()
			
	def run(self):
			
		self.molecular_layer()
		self.snp_layer()
		self.blk_layer()
		self.propagate()
		
	#running each cycle of network
	def cycle(self, int cycle1, int cycle2):
					
		for icycle in xrange(cycle1):
		
			#print "turn: ", icycle
			for i in xrange(cycle2):
			
				self.molecular_layer()
				self.snp_layer()
				self.blk_layer()
				self.propagate()
			
			self.phase()
