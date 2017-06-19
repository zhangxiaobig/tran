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

	exmp = "Example: Tran --input RNAseq1.bam RNAseq2.bam  --TE TE_annotation.gtf --sortByPos"

	parser = argparse.ArgumentParser(prog='zx',description=desc, epilog=exmp) 

	parser.add_argument('--input', metavar='sample file', dest='files',nargs='+', required=True,
						help='Sample files')
	parser.add_argument('--exp_level', metavar='TE expression level', dest='exp_level',type=str, nargs='?', default='element',choices=['instance','element','family','class'],
						help='Calculate TE expression on ? level:(instance, element/subfamily, family or class). DEFAULT: element.')
	parser.add_argument('--TE', metavar='TE-GTF-file', dest='tefile', type=str, required=True,
						help='GTF file for transposable element annotations')
	parser.add_argument('--format', metavar='input file format', dest='format', type=str, nargs='?', default='BAM', choices=['BAM','SAM'],
						help='Input file format: BAM or SAM. DEFAULT: BAM')
	parser.add_argument('--stranded', metavar='option', dest='stranded', nargs='?', type=str, default="yes", choices=['yes','no','reverse'],
						help='Is this a stranded library? (yes, no, or reverse). DEFAULT: yes.')
	parser.add_argument('-o','--output', metavar='output directory', dest='output_dir', nargs='?', default='ZX_uniq_out',
						help='Name of output directory. DEFAULT: tran_out')	   
	parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
	
	args = parser.parse_args()

	return args
	
#def parameter_log():

def calulate_abundance(samfile,TE,maxL=500,numItr=10):
	references = samfile.references
	
	empty = 0
	nonunique = 0
	uniq_reads = 0

	i = 0
	prev_read_name =''

	alignments_per_read = []	
	avgReadLength =	 0
	tmp_cnt = 0
	multi_read1 = []
	multi_read2 = []
	multi_reads = []
	uniq_counts=[0.0]*TE.numInstances()
	multi_counts = [0.0]*TE.numInstances()
	readMappability = []
	TE_position_list =[]
	
	for aligned_read in samfile.fetch():
		
		if aligned_read.is_unmapped or aligned_read.is_duplicate or aligned_read.is_qcfail :				
			continue		   
		
		cur_read_name = aligned_read.query_name

		if aligned_read.is_paired :	 # this is not reliable if read mate is unmapped
			if cur_read_name == prev_read_name or prev_read_name == "":
				prev_read_name = cur_read_name
				if aligned_read.is_read1 :
						multi_read1.append(aligned_read)
				if aligned_read.is_read2 :
						multi_read2.append(aligned_read)
				continue
				
			else :
				if len(multi_read1) <= 1 and len(multi_read2) <= 1 :
					uniq_reads += 1
					read1 = None
					read2 = None
					if len(multi_read1) == 1:
						read1 = multi_read1[0]
					if len(multi_read2) == 1:
						read2 = multi_read2[0]
					if read1 is not None and read1.is_proper_pair :
						if read2 is None :
							sys.stderr.write("******NOT COMPLETE*******\n")
							sys.stderr.write("If the BAM file is sorted by coordinates, please use unsorted file and re-run!\n")
							sys.exit(0)
						if tmp_cnt < 10000 :
							pos1 = read1.reference_start
							pos2 = read2.reference_start
							if abs(pos1-pos2) <= maxL :
								avgReadLength += abs(pos1-pos2)+read2.query_length
							tmp_cnt += 1

					alignments_per_read.append((read1,read2))
					
						
				else:
					nonunique += 1
					#singlton
					if len(multi_read2) == 0 :
						for r in multi_read1 :
							alignments_per_read.append((r,None))
					if len(multi_read1) == 0:
						for r in multi_read2 :
							alignments_per_read.append((None,r))
					if len(multi_read2) == len(multi_read1) :
						for i in range(len(multi_read1)) :
							read1 = multi_read1[i]
							read2 = multi_read2[i]
							alignments_per_read.append((read1,read2))
							
		else : #single end read
			if cur_read_name == prev_read_name or prev_read_name == "":
				alignments_per_read.append((aligned_read,None))
				prev_read_name = cur_read_name
				continue
			else : # a new read
				if tmp_cnt < 10000 :
					avgReadLength += aligned_read.query_length
					tmp_cnt += 1

				if len(alignments_per_read) == 1 :
					uniq_reads += 1
				else :
					nonunique += 1
		
		(TE_position_ratio,Mappability,TE_position)= reads_ovp_TE(references,alignments_per_read,TE)
		TE_ovp_reads(TE_position_ratio,Mappability,TE_position,uniq_counts,multi_counts,multi_reads,readMappability,TE_position_list)
		
		alignments_per_read = []
		multi_read1 = []
		multi_read2 = []
		prev_read_name = cur_read_name
		if not aligned_read.is_paired :
			alignments_per_read.append((aligned_read,None))
		else :
			if aligned_read.is_read1 :
				multi_read1.append(aligned_read)
			if aligned_read.is_read2 :
				multi_read2.append(aligned_read)
	
	sys.stderr.write("processing reads done!\n")
	ISOTIMEFORMAT='%Y-%m-%d %X'
	print time.strftime(ISOTIMEFORMAT,time.localtime(time.time()))
	'''
	te_multi_counts=[0] * len(uniq_counts)	
	
	if avgReadLength > 0 :
		avgReadLength = int(avgReadLength/tmp_cnt)
		
	if len (multi_reads) > 0 :
		new_te_multi_counts = EM(TE,multi_reads,uniq_counts,te_multi_counts,numItr,avgReadLength,multi_counts,readMappability,TE_position_list)
										
	te_counts = map(operator.add,uniq_counts,new_te_multi_counts)
	sys.stderr.write("EM done!\n")
	print time.strftime(ISOTIMEFORMAT,time.localtime(time.time()))
	'''
	return uniq_counts
	
class TEindex():

	def __init__ (self,filename):
		self.__namelist = collections.defaultdict(dict)
		self.__window = collections.defaultdict(dict)
		self._length = []
		self._nameIDmap = []
		self._elements = []
		self._start = []
		self._end= []
		
		f = open(filename,'r')

		name_idx = 0
		linenum = 0
		cur_chrom = 'chr1'
		chr_name_idx = 0
		cur_name_idx = 0
		cur_start = 0
		cur_end = 0
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
			if ele_name not in self._elements :
					self._elements.append(ele_name)
												   
			self._length.append(tlen)
			self._nameIDmap.append(full_name)	
				
			self.__namelist[chrom][name_idx] = [(start,end)]		  
				
			if chrom == cur_chrom:
				if (name_idx-chr_name_idx)%500 == 0 :						 
					#print chrom,cur_name_idx,cur_start,end
					cur_start=start
					cur_name_idx=name_idx
			else :
				chr_name_idx=name_idx
			
			self.__window[chrom][cur_name_idx] = [(cur_start,end)]	
				
			cur_chrom = chrom  
			cur_end = end				
			name_idx += 1			   
				
		f.close()

	def getNames(self) :
		return self.__namelist
		
	def getElements(self) :
		 return self._elements
		 
	def getStrand(self,idx) :
		f_name = self._nameIDmap[idx]
		return f_name[len(f_name)-1]
		
	def numInstances(self) :
		return len(self._nameIDmap)
		
	def getWindow(self) :
		return self.__window

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
														
			
def reads_ovp_TE(references,alignments_per_read,TE):
	TE_position = 0
	Mappability=0
	TE_position_ratio={}
	#loop over every alignment
	for r1,r2 in alignments_per_read :
		itv_list = []
		TE_count=0
		direction = 1
		if r1 is not None and r1.is_reverse :
			direction = -1
		if r2 is not None and not r2.is_reverse :
			direction = -1
		chrom1 = ""
		if r1 is not None :
			chrom1 = references[r1.tid]
		chrom2 = ""
		if r2 is not None :
			chrom2 = references[r2.tid]

		#fetch all mapping intervals
		itv_list = []
		if r1 is not None :
			itv_list = fetch_exon(chrom1,r1.pos,r1.cigartuples,direction,format,r1.get_tag("NH"))###########################atention
   
		if chrom2 != "" : #paired-end, both mates are mapped
			itv_list2 = fetch_exon(chrom2,r2.pos,r2.cigar,direction,format,r2.get_tag("NH"))############################
			itv_list.extend(itv_list2)
		
		TEnamelist = []
		TE_ratio={}		
		for iv in itv_list :		
			for keys in TE.getWindow()[iv[0]].keys() :
				[(s,e)] = TE.getWindow()[iv[0]][keys]						  
				if iv[1] <= e and iv[2] >= s :
					try :
						for i in range(499):
							[(s,e)] = TE.getNames()[iv[0]][keys+i]
							if iv[1] <= e and iv[2] >= s :
								TEnamelist.append(keys+i)								
								break	
					except KeyError :
						pass
					
					break
					
			if len(TEnamelist)>0:				
				for t in TEnamelist :
					if iv[3] == TE.getStrand(t) :
						TE_count=1
						TE_ratio[t]=1/float(iv[4])
						Mappability=iv[4]									
						
		TE_position +=TE_count	
	
	if len(TE_ratio.keys())> 0 :
		for k in TE_ratio:
			if k not in TE_position_ratio:
				TE_position_ratio[k]=0.0			
				TE_position_ratio[k] +=TE_ratio[k]
			else :
				TE_position_ratio[k] +=TE_ratio[k]

		#print r1.qname,TE_list,TE_position_ratio_list,TE_position,Mappability
					
	return (TE_position_ratio,Mappability,TE_position)


def fetch_exon(chrom, st, cigar,direction,format,readMappability):
	
	chrom_st = st
	if format == "BAM" :
		chrom_st += 1
	exon_bound =[]

	for c,s in cigar:	 #code and size
		if c==0:			#match			
			chrom_st += s
		elif c==1:		  #insertion to ref
			continue
		elif c==2:		  #deletion to ref
			chrom_st += s
		elif c==3:		  #gap or intron
			chrom_st += s
		elif c==4:		  #soft clipping. We do NOT include soft clip as part of exon
			chrom_st += s
		else:
			continue
	
	if direction == 1 :
		exon_bound.append([chrom, chrom_st,chrom_st + s-1,"+",readMappability])
	if direction == -1 :
		exon_bound.append([chrom, chrom_st,chrom_st + s-1,"-",readMappability])
			
	return exon_bound
	
def TE_ovp_reads(TE_position_ratio,Mappability,TE_position,uniq_counts,multi_counts,multi_reads,readMappability,TE_position_list):	
	if	Mappability == 1 :
		for i in TE_position_ratio.keys():
			uniq_counts[i] += 1

	if	Mappability > 1 :
		multi_algn=[]
		for te in TE_position_ratio:
			multi_counts[int(te)] +=TE_position_ratio[te]
			multi_algn.append(te)
		multi_reads.append(multi_algn)
		readMappability.append(Mappability)
		TE_position_list.append(TE_position)
		
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
		for te in sorted(cnt_tbl.keys()) :
			vals = cnt_tbl[te]
			vals_str = te
			for i in range(len(vals)) :
				vals_str +="\t"+str(vals[i])
			f.write(vals_str + "\n")

		f.close()

	return
				
def main():
	ISOTIMEFORMAT='%Y-%m-%d %X'
	sys.stderr.write("BEGIN\n")
	print time.strftime(ISOTIMEFORMAT,time.localtime(time.time()))
	
	args=parameters()
	TE=TEindex(args.tefile)
	
	sys.stderr.write("TE index building done\n")
	print time.strftime(ISOTIMEFORMAT,time.localtime(time.time()))
	
	for filename in args.files :
		cnt_tbl = {}
		if format == "BAM" :
			samfile = pysam.AlignmentFile(filename,'rb')
		else :
			samfile = pysam.AlignmentFile(filename,'r')
		te_instance_counts=calulate_abundance(samfile,TE)
		te_ele_counts = TE.groupByEle(te_instance_counts)
		cnt_tbl[filename]=dict(te_ele_counts.items())		
		
	f_cnt_tbl = args.output_dir + ".cntTable"
	output_count_tbl(cnt_tbl, f_cnt_tbl)	
if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt !\n")
		sys.exit(0)				