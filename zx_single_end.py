#!/usr/bin/env python2.7
'''
@auther: xiao zhang
@contact:zhangxiao@big.ac.cn
@version:1.0.0
'''

#python module
import argparse
import collections
import sys
import pysam
import time

from zx_EM import *
		
#parameters
def parameters():

	desc = "Identifying transcription of transposable elements."

	exmp = "Example: Tran --input RNAseq1.bam RNAseq2.bam  --TE TE_annotation.gtf "

	parser = argparse.ArgumentParser(prog='zx',description=desc, epilog=exmp) 

	parser.add_argument('--input', metavar='sample file', dest='files',nargs='+', required=True,
						help='Sample files')
	parser.add_argument('--exp_level', metavar='TE expression level', dest='exp_level',type=str, nargs='?', default='instance',choices=['instance','element','family','class'],
						help='estimate TE expression on ? level:(instance, element/subfamily, family or class). DEFAULT: instance.')
	parser.add_argument('--TE', metavar='TE-GTF-file', dest='tefile', type=str, required=True,
						help='GTF file for transposable element annotations')
	parser.add_argument('--format', metavar='input file format', dest='format', type=str, nargs='?', default='SAM', choices=['BAM','SAM'],
						help='Input file format: BAM or SAM. DEFAULT: BAM')
	parser.add_argument('--type', metavar='reads type', dest='reads', type=str, nargs='?', default='multi' , choices=['multi','uniq'],
						help='estimate TE expression by ? reads:(uniq or multi). DEFAULT: multi.')
	parser.add_argument('--stranded', metavar='option', dest='stranded', nargs='?', type=str, default="no", choices=['yes','no','reverse'],
						help='Is this a stranded library? (yes, no, or reverse). DEFAULT: yes.')
	parser.add_argument('-o','--output', metavar='output directory', dest='output_dir', nargs='?', default='ZX_copy_out',
						help='Name of output directory. DEFAULT: ZX_copy_out')	   
	parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
	
	args = parser.parse_args()

	return args

def calulate_abundance(samfile,TE,reads_type,maxL=500,numItr=50):
	references = samfile.references

	nonunique = 0
	uniq_reads = 0

	prev_read_name =''

	alignments_per_read = []	
	avgReadLength =	 0
	tmp_cnt = 0
	
	multi_reads = []
	uniq_counts=[0.0]*TE.numInstances()
	multi_counts = [0.0]*TE.numInstances()
	readMappability = []
	TE_position_list =[]
	TE_position_ratio_list=[]
	Mappability_list=[]
	uniq_counts_com={}
	
	for aligned_read in samfile.fetch():
		
		if aligned_read.is_unmapped or aligned_read.is_duplicate or aligned_read.is_qcfail :				
			continue		   
		
		cur_read_name = aligned_read.query_name

		if aligned_read.is_paired :	 
			break
			
		else : #single end read
			if prev_read_name == "":
				if aligned_read.get_tag("NH")==1:
					alignments_per_read.append((aligned_read,None))
					continue
				else:
					alignments_per_read.append((aligned_read,None))
					prev_read_name = cur_read_name
					continue
					
			elif cur_read_name == prev_read_name :
				alignments_per_read.append((aligned_read,None))
				#print alignments_per_read
				prev_read_name = cur_read_name
				continue
				
			else : # a new read					
				if tmp_cnt < 10000 :
					avgReadLength += aligned_read.query_length
					tmp_cnt += 1
					
				#print alignments_per_read			
				#print aligned_read.query_name
				(TE_position_ratio,Mappability,TE_position)= reads_ovp_TE(references,alignments_per_read,TE,aligned_read.query_length)
				
				if	Mappability == 1 :
					for i in TE_position_ratio:
						if type(i)==int:
							uniq_counts[i] += 1
						else:
							if i not in uniq_counts_com:
								uniq_counts_com[i]=1
							else:
								uniq_counts_com[i]+=1
								
				TE_position_ratio_list.append(TE_position_ratio)
				Mappability_list.append(Mappability)

		
		prev_read_name = cur_read_name
		alignments_per_read = []
		alignments_per_read.append((aligned_read,None))
	
	#for c in uniq_counts_com:
		#print c,uniq_counts_com[c]
		
	TE_combine_trackback=TE_ovp_reads(TE,TE_position_ratio_list,Mappability_list,uniq_counts,uniq_counts_com,multi_counts,multi_reads,readMappability)	
	
	te_counts=[0] * len(uniq_counts)	
	te_multi_counts=[0] * len(uniq_counts)
	
	if avgReadLength > 0 :
		avgReadLength = int(avgReadLength/tmp_cnt)
		
	if len (multi_reads) > 0 :
		multi_reads_within_uniq_reads(TE,multi_reads,uniq_counts,te_multi_counts,TE_combine_trackback,uniq_counts_com)
		#te_multi_counts = EM(TE,multi_reads,numItr,avgReadLength,multi_counts,readMappability,TE_position_list)
		#print uniq_counts
		#for i in range(len(te_multi_counts)):
		#	if i <TE.numInstances():
		#		te_counts[i]=uniq_counts[i]+te_multi_counts[i]
		#	else:
		#		for t in TE_combine_trackback[i]:
		#			if TE_combine_trackback[i] in uniq_counts_com:
		#				te_counts[t]=uniq_counts[t]+te_multi_counts[i]+uniq_counts_com[TE_combine_trackback[i]]
		#			else:
		#				te_counts[t]=uniq_counts[t]+te_multi_counts[i]
		
	te_counts = map(operator.add,uniq_counts,te_multi_counts)
	return te_counts		

class TEindex():

	def __init__ (self,filename):
		self.__namelist = collections.defaultdict(dict)
		self.__IDmap = collections.defaultdict(dict)
		self._length = []
		self._nameIDmap = []
		self._copynameIDmap = []
		self._elements = []
		self._start = []
		self._end= []
		
		f = open(filename,'r')

		name_idx = 0
		idx = 0
		linenum = 0
		cur_chrom = ''
		chr_name_idx = 0
		cur_name_idx = 0
		cur_start = 0

		for line in f :
			line = line.strip()
			items = line.split('\t')
			chrom = items[0]
			start = int(items[3])
			end = int(items[4])
			strand = items[6]
			items[8] = items[8].replace("; ",";")
			desc = items[8].split(';')
			name = ""
			family_id = ""
			ele_id = ""
			class_id = ""
			tlen = end - start + 1
			linenum += 1

			for i in range(len(desc)) :
				desc[i] = desc[i].replace("\"","")
				pos = desc[i].find(" ")
				tid = desc[i][:pos]
				val = desc[i][pos+1:len(desc[i])]

				if tid == "gene_id" :
					ele_id = val
				if tid == "transcript_id" :
					name = val
				if tid == "family_id" :
					family_id = val
				if tid == "class_id" :
					class_id = val				  

			full_name = name+':'+ele_id+':'+family_id+':'+class_id+':'+strand
			ele_name = ele_id+':'+family_id+':'+class_id
			full_copy_name = name+':'+ele_id+':'+family_id+':'+class_id+':'+str(start)+'-'+str(end)
			
			if ele_name not in self._elements :
					self._elements.append(ele_name)
												   
			self._length.append(tlen)
			self._nameIDmap.append(full_name)
			self._copynameIDmap.append(full_copy_name)	

			if cur_chrom!=chrom :
				name_idx=0
				
			self.__namelist[chrom][name_idx] = [(start,end)]
			self.__IDmap [chrom][name_idx] = idx
			cur_chrom=chrom
			name_idx += 1
			idx +=1
				
		f.close()
	def getID(self,chrom,name_idx):
		return self.__IDmap [chrom][name_idx]
		
	def getNames(self) :
		return self.__namelist	
		
	def getElements(self) :
		return self._elements
		 
	def getFullCopyName(self) :
		return self._copynameIDmap
		
	def getStrand(self,idx) :
		f_name = self._nameIDmap[idx]
		return f_name[len(f_name)-1]
		
	def numInstances(self) :
		return len(self._nameIDmap)

	def getEleName(self,idx) :
		full_name = None
		if idx >= len(self._nameIDmap) or idx < 0 :
			return None
		else :
			full_name =	 self._nameIDmap[idx]
		if full_name is not None:
			pos = full_name.find(':')
			val = full_name[pos+1:(len(full_name)-2)]
			return val
		else :
			return None
			

	def getFullName(self,idx) :
		if idx >= len(self._nameIDmap) or idx < 0 :
			return None
		else :
			return self._nameIDmap[idx]

	def getLength(self,idx) :
		if idx < len(self._length) :
			return self._length[idx]
		else :
			return -1
			
	def groupByEle(self,te_inst_counts) :
		TEs = self.getElements()
		te_ele_counts = dict(zip(TEs,[0]*len(TEs)))

		for i in range(len(te_inst_counts)) :
			ele_name = self.getEleName(i)
			if ele_name in te_ele_counts :
				te_ele_counts[ele_name] += te_inst_counts[i]

		return te_ele_counts
		
	def groupByCopy(self,te_inst_counts) :
		names= self.getFullCopyName()
		te_copy_counts = dict(zip(names,[0]*len(names)))

		for i in range(len(te_inst_counts)) :
			copy_name = self._copynameIDmap[i]
			te_copy_counts[copy_name] = te_inst_counts[i]
			#print te_copy_counts
		return te_copy_counts

def multi_reads_within_uniq_reads(TE,multi_reads,uniq_counts,te_multi_counts,TE_combine_trackback,uniq_counts_com):
	for i in multi_reads:
		te_num=[]
		for te in i:
			if te <TE.numInstances():
				if uniq_counts[te]!=0:
					te_num.append(te)
					
			else:
				tee = TE_combine_trackback[te]
				if tee in uniq_counts_com:
					te_num.append(te)
				
				
		if len(te_num)==1:
			if te_num[0]<TE.numInstances():
				te_multi_counts[int(te_num[0])]+=1	
			else:
				for t in TE_combine_trackback[te_num[0]]:
					te_multi_counts[int(t)]+=1	
	
			
def reads_ovp_TE(references,alignments_per_read,TE,length):
	TE_position = 0
	Mappability=0
	TE_position_ratio={}
	#loop over every alignment
	for r1,r2 in alignments_per_read :
		itv_list = []
		TE_count=0
		direction = "+"
		if r1 is not None and r1.is_reverse :
			direction = "-"
		if r2 is not None and not r2.is_reverse :
			direction = "-"
		chrom1 = ""
		if r1 is not None :
			chrom1 = references[r1.tid]
		chrom2 = ""
		if r2 is not None :
			chrom2 = references[r2.tid]

		#fetch all mapping intervals
		itv_list = []
		
		if r1 is not None :
			itv_list.append([chrom1,r1.reference_start,r1.reference_end,direction,r1.get_tag("NH")])###########################atention
   
		if chrom2 != "" : #paired-end, both mates are mapped
			itv_list2.append([chrom2,r2.reference_start,r2.reference_end,direction,r2.get_tag("NH")])############################
			itv_list.extend(itv_list2)
	
		TEnamelist = []
		TE_ratio=0
		for iv in itv_list :
			#print iv
			low = 0 
			high=len(TE.getNames()[iv[0]])-1
			middle = 0
					
			while(low <= high):
				middle = int((low + high)/2)
				[(s,e)] = TE.getNames()[iv[0]][middle]
				if iv[1]+1 <= e and iv[2] >= s :
					TEnamelist.append(middle)
					cur_position=middle
					[(s,e)] = TE.getNames()[iv[0]][cur_position]
					while (e >= iv[1]+1):
						if cur_position==0:
							break
						cur_position =cur_position-1
						[(s,e)] = TE.getNames()[iv[0]][cur_position]
						if iv[1]+1 <= e and iv[2] >= s :
							TEnamelist.append(cur_position)
													
					cur_position=middle
					[(s,e)] = TE.getNames()[iv[0]][cur_position]
					while (s <= iv[2]):
						if cur_position==len(TE.getNames()[iv[0]])-1:
							break
						cur_position =cur_position+1
						[(s,e)] = TE.getNames()[iv[0]][cur_position]
						if iv[1]+1 <= e and iv[2] >= s :
							TEnamelist.append(cur_position)
						
					break
					
				elif e < iv[1]+1:
					low = middle+1

				elif s > iv[2]:
					high = middle-1
								
			#print len(TEnamelist)	
			if len(TEnamelist)>0:				
				for t in TEnamelist :
					TE_count=1
					Mappability=iv[4]	
					TE_ratio=1/float(Mappability)
					s=TE.getID(iv[0],t)
					if s not in TE_position_ratio:		
						TE_position_ratio[s] =TE_ratio
					else :
						TE_position_ratio[s] +=TE_ratio
						
				if len(TEnamelist)>1:
					ss=[]
					for t in TEnamelist :
						s=TE.getID(iv[0],t)
						ss.append(s)
						
					ss=tuple(ss)	
					TE_position_ratio[ss] =TE_ratio
			
			TEnamelist = []										
			TE_ratio = 0
			
		TE_position +=TE_count	
	#print TE_position_ratio	
	
	return (TE_position_ratio,Mappability,TE_position)

	
def TE_ovp_reads(TE,TE_position_ratio_list,Mappability_list,uniq_counts,uniq_counts_com,multi_counts,multi_reads,readMappability):
	TE_combine={}
	TE_combine_trackback={}
	num=int(TE.numInstances()-1)
	for i in range(len(Mappability_list)):
		if	Mappability_list[i] > 1 :
			multi_algn=[]
			mark=0
			for te in TE_position_ratio_list[i]:				
				if type(te)==int:					
					if uniq_counts[te]!=0:
						multi_counts[int(te)] +=1
						mark=1
						multi_algn.append(te)
						
				else :
					if te not in TE_combine:
						num +=1
						TE_combine[te]=num
						TE_combine_trackback[num]=te
						multi_counts.append([0.0])
						multi_counts[TE_combine[te]] =0.0
						if te in uniq_counts_com:						
							multi_counts[TE_combine[te]] +=1
							mark=1
							tee=TE_combine[te]
							te=tee
							multi_algn.append(te)	
					
				
			
			if mark==0:
				for te in TE_position_ratio_list[i]:
					if type(te)==int:
						multi_counts[int(te)] +=TE_position_ratio_list[i][te]
					else:
						multi_counts[TE_combine[te]] +=TE_position_ratio_list[i][te]						
						tee=TE_combine[te]
						te=tee
						
					multi_algn.append(te)
			#print multi_algn		
			multi_reads.append(multi_algn)
			readMappability.append(Mappability_list[i])	
			#print multi_algn
	return TE_combine_trackback
			
def output_count_tbl(t_tbl,fname):
	try:
		f = open(fname, 'w')
	except IOError:
		error("Cannot create report file %s !\n" % (fname))
		sys.exit(1)
	else:
		cnt_tbl = {}
		header = "TE"
		keys = set([])
		for tsmp in t_tbl.keys():
			keys = keys.union(t_tbl[tsmp].keys())
			header +="\t"+tsmp
			
		for tsmp in t_tbl.keys():
			cnts = t_tbl[tsmp]
			for k in keys:
				val = 0
				if k in cnts :
					val = cnts[k]
				if cnt_tbl.has_key(k) :
					cnt_tbl[k].append(int(val))
				else :
					vallist = []
					vallist.append(int(val))
					cnt_tbl[k] = vallist

		#output
		f.write(header + "\n")
		for te in sorted(cnt_tbl.keys()):
			vals = cnt_tbl[te]
			vals_str =te
			for i in range(len(vals)) :
				vals_str +="\t"+str(vals[i])
			f.write(vals_str + "\n")

		f.close()

	return			
			
	
def main():
	#ISOTIMEFORMAT='%Y-%m-%d %X'
	#sys.stderr.write("BEGIN\n")
	#print time.strftime(ISOTIMEFORMAT,time.localtime(time.time()))
	
	args=parameters()
	TE=TEindex(args.tefile)
	reads_type=args.reads
	
	#sys.stderr.write("TE index building done\n")
	#print time.strftime(ISOTIMEFORMAT,time.localtime(time.time()))
	
	for filename in args.files :
		cnt_tbl = {}
		if format == "BAM" :
			samfile = pysam.AlignmentFile(filename,'rb')
		else :
			samfile = pysam.AlignmentFile(filename,'r')
		te_instance_counts=calulate_abundance(samfile,TE,reads_type)
		te_lev_counts = TE.groupByCopy(te_instance_counts)
		#print te_lev_counts
		cnt_tbl[filename]=dict(te_lev_counts.items())		
		
	f_cnt_tbl = args.output_dir + ".cntTable"
	output_count_tbl(cnt_tbl, f_cnt_tbl)	
	
	
if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt !\n")
		sys.exit(0)				