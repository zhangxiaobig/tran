#!/usr/bin/env python2.7
'''
@auther: xiao zhang
@contact:zhangxiao@big.ac.cn
@version:
'''

# python module
import pdb
import sys
import os.path
import math
import operator
import argparse
import traceback
import collections
import subprocess
from time import time
import pysam
import logging

from constants import *
from EMAlgorithm import *
from Index import *
from Input import *

#parameters
def prepare_parser():
    desc = "Identifying differential transcription of gene and transposable elements."

    exmp = "Example: Tran -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --GTF gene_annotation.gtf --TE TE_annotation.gtf --sortByPos"

    parser = argparse.ArgumentParser(prog='Tran',description=desc, epilog=exmp) 

    parser.add_argument('-t','--treatment', metavar='treatment sample', dest='tfiles',nargs='+', required=True,
                        help='Sample files in group 1 (e.g. treatment/mutant)')
    parser.add_argument('-c','--control', metavar='control sample', dest='cfiles',nargs='+', required=True,
                        help='Sample files in group 2 (e.g. control/wildtype)')
    parser.add_argument('--GTF', metavar='genic-GTF-file', dest='gtffile', type=str, required=True,
                        help='GTF file for gene annotations')
    parser.add_argument('--TE', metavar='TE-GTF-file', dest='tefile', type=str, required=True,
                        help='GTF file for transposable element annotations')
    parser.add_argument('--format', metavar='input file format', dest='format', type=str, nargs='?', default='BAM', choices=['BAM','SAM'],
                        help='Input file format: BAM or SAM. DEFAULT: BAM')
    parser.add_argument('--stranded', metavar='option', dest='stranded', nargs='?', type=str, default="yes", choices=['yes','no','reverse'],
                        help='Is this a stranded library? (yes, no, or reverse). DEFAULT: yes.')
    parser.add_argument('--project', metavar='name', dest='prj_name', nargs='?', default='tran_out',
                        help='Name of this project. DEFAULT: tran_out')
    parser.add_argument('--sortByPos', dest = 'sortByPos', action="store_true",
                        help='Alignment files are sorted by chromosome position.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

    return parser


class UnknownChrom(Exception):
    pass

def count_reads(samples, format, geneIdx, teIdx,stranded, sortByPos, prj_name,numItr,fragLength,maxL):
    cnt_tbl = {}
    libsize = []
    # check input files exist or notb
    for filename in samples :
        if not os.path.isfile(filename):
            sys.stderr.write("File %s does not exist or is not a file.\n" % (filename))
            sys.exit(1)
    try:
        for filename in samples :
            (gene_counts,te_instance_counts) = count_transcript_abundance(filename,format,geneIdx,teIdx,stranded,sortByPos, prj_name,numItr,fragLength,maxL)
            
            #summarize on elements
            te_ele_counts = teIdx.groupByEle(te_instance_counts)
            # save gene counts and TE counts into count table
            cnt_tbl[filename] = dict(gene_counts.items() + te_ele_counts.items())
            num_reads = sum(gene_counts.values())+sum(te_ele_counts.values())
            libsize.append(num_reads)
    except:
        sys.stderr.write("Error: %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" %
                          ( sys.exc_info()[1].__class__.__name__,
                            os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]),
                            traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
        sys.exit(1)

    return (cnt_tbl,libsize)

def fetch_exon(chrom, st, cigar,direction,format,readMappability):
    ''' fetch exon regions defined by cigar. st must be zero based
    return list of tuple of (chrom,st, end)
    '''
    #match = re.compile(r'(\d+)(\D)')
    chrom_st = st
    if format == "BAM" :
        chrom_st += 1
    exon_bound =[]

    for c,s in cigar:    #code and size
        if c==0:        #match
            if direction == 0 :
                exon_bound.append([chrom, chrom_st,chrom_st + s-1,".",readMappability])
            if direction == 1 :
                exon_bound.append([chrom, chrom_st,chrom_st + s-1,"+",readMappability])
            if direction == -1 :
                exon_bound.append([chrom, chrom_st,chrom_st + s-1,"-",readMappability])
            chrom_st += s
        elif c==1:        #insertion to ref
            continue
        elif c==2:        #deletion to ref
            chrom_st += s
        elif c==3:        #gap or intron
            chrom_st += s
        elif c==4:        #soft clipping. We do NOT include soft clip as part of exon
            chrom_st += s
        else:
            continue
    return exon_bound

def ovp_annotation(references,multi_reads,geneIdx,teIdx,stranded,format) :
    annot_gene = []
    annot_TE = []
    readMappability = []
    #loop over every alignment
    for r1,r2 in multi_reads :
        itv_list = []
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

        if stranded == "no":
            direction = 0
        if stranded == "reverse":
            direction = direction *(-1)

        #fetch all mapping intervals
        itv_list = []
        if r1 is not None :
            itv_list = fetch_exon(chrom1,r1.pos,r1.cigartuples,direction,format,r1.get_tag("NH"))###########################atention
   
        if chrom2 != "" : #paired-end, both mates are mapped
            itv_list2 = fetch_exon(chrom2,r2.pos,r2.cigar,direction,format,r2.get_tag("NH"))############################
            itv_list.extend(itv_list2)
        try:
            (TEs,Mappability) = teIdx.TE_annotation(itv_list)
            genes = geneIdx.Gene_annotation(itv_list)

            if len(TEs) > 0   :
                annot_TE.append(TEs)
                #print TEs
                
            if len(Mappability) > 0:
                readMappability.append(Mappability)
                #print Mappability
            if len(genes) > 0 :
                annot_gene.append(list(set(genes)))
        except:
            sys.stderr.write("Error occurred during read assignments\n")
            raise
    #print annot_TE,readMappability
    return (annot_gene,annot_TE,readMappability)
    
def readInAlignment(filename, format, sortByPos, prj_name):
    try:
        if format == "BAM" :
            if not sortByPos  :
                samfile = pysam.AlignmentFile(filename,'rb')
                if len(samfile.header) ==0:
                    print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                    sys.exit(1)

            else :
                cur_time = time.time()
                bam_tmpfile = '.'+str(cur_time)
                pysam.sort("-n",filename,bam_tmpfile)
                samfile = pysam.AlignmentFile(bam_tmpfile+".bam",'r')
                subprocess.call(['rm -f ' + bam_tmpfile+'.bam' ],shell=True)

        else :
                samfile = pysam.AlignmentFile(filename,'r')

    except:
        sys.stderr.write("Error occured when reading first line of sample file %s.\n" % filename)
        raise
    return samfile  

def count_transcript_abundance(filename, format, geneIdx, teIdx, stranded, sortByPos, prj_name,numItr,fragLength,maxL):
    #read in BAM/SAM file
    samfile = readInAlignment(filename, format, sortByPos,prj_name)

    references = samfile.references 
    
    empty = 0
    nonunique = 0
    uniq_reads = 0

    i = 0
    prev_read_name =''

    alignments_per_read = []
    leftOver_gene = []
    leftOver_te = []
    avgReadLength =  0
    tmp_cnt = 0
    multi_read1 = []
    multi_read2 = []
    paired = False
    gene_counts = dict(zip(geneIdx.getFeatures(),[0.0]*len(geneIdx.getFeatures()))) 
    te_counts = [0.0] * teIdx.numInstances()
    te_multi_counts = [0.0]*len(te_counts)
    multi_reads = []
    annot_TE_position = [[0 for col in range(1)] for row in range(len(te_counts))]
    readMappability = [0.0]*len(te_counts)
    cc = 0
    try:
        while(1):
            i += 1
            cc += 1
            aligned_read = samfile.next()
            if aligned_read.is_unmapped or aligned_read.is_duplicate or aligned_read.is_qcfail :                
                continue
           
            cur_read_name = aligned_read.query_name

            if aligned_read.is_paired :  # this is not reliable if read mate is unmapped
                added_flag_pos = aligned_read.query_name.find('/')
                
                if added_flag_pos == -1 :
                    added_flag_pos = len(aligned_read.query_name)
                if aligned_read.query_name.endswith(".1") or aligned_read.query_name.endswith(".2") :
                    added_flag_pos = len(aligned_read.query_name) - 2
                
                cur_read_name = aligned_read.query_name[:added_flag_pos]
                
                paired = True
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
                                sys.stderr.write("If the BAM file is sorted by coordinates, please specify --sortByPos and re-run!\n")
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
                            
        
               
            try :
                
                (annot_gene,annot_TE,Mappability) = ovp_annotation(references,alignments_per_read, geneIdx, teIdx,stranded,format)
                
                if len(alignments_per_read) > 1 : #multi read, prior to TE
                    no_annot_te = parse_annotations_TE(multi_reads,annot_TE,annot_TE_position,Mappability, te_counts, te_multi_counts,readMappability, leftOver_te)

                    if no_annot_te :
                        no_annot_gene = parse_annotations_gene(annot_gene,gene_counts,leftOver_gene)
                        if no_annot_gene :
                            empty += 1

                else : #uniq read , prior to gene
                    no_annot_gene = parse_annotations_gene(annot_gene,gene_counts,leftOver_gene)
                    if no_annot_gene :
                        no_annot_te = parse_annotations_TE(multi_reads,annot_TE,annot_TE_position,readMappability,te_counts, te_multi_counts, readMappability,leftOver_te)
                        if no_annot_te :
                            empty += 1
                
            except:
                sys.stderr.write("Error occurred when processing annotations of %s in file %s.\n" % (prev_read_name, filename))
                raise

            if i % 1000000 == 0 :
                    sys.stderr.write("%d %s processed.\n" % (i, "alignments " ))

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
      
            
    except StopIteration:
            pass
                #the last read  
          
    try:
            #resolve ambiguity
            if len(leftOver_gene) > 0 :
                resolve_annotation_ambiguity(gene_counts,leftOver_gene)
            if len(leftOver_te) > 0 :
                resolve_annotation_ambiguity(te_counts,leftOver_te)

            ss = sum(te_counts)
            sys.stderr.write("uniq te counts = %s \n" % (str(ss)))

            te_tmp_counts = [0]*len(te_counts)

            if numItr > 0 :
              
              try :
                '''iterative optimization on TE reads '''
                sys.stderr.write(".......start iterative optimization ..........\n")
                if not paired and fragLength > 0:
                    avgReadLength = fragLength

                elif avgReadLength > 0 :
                    avgReadLength = int(avgReadLength/tmp_cnt)
                else :
                    sys.stderr.write("There are not enough reads to estimate fragment length. \n" )
                    raise 
                new_te_multi_counts = [0] *len(te_counts)
                                
                if len(multi_reads) > 0 :
                    new_te_multi_counts = EMestimate(teIdx,multi_reads,te_tmp_counts,te_multi_counts,avgReadLength,readMappability,annot_TE_position)
                    #print '######################'teIdx,'\t',multi_reads,'\t',te_tmp_counts,'\t',te_multi_counts,'\t',numItr,'\t',avgReadLength ###########
              except :
                sys.stderr.write("Error in optimization\n")
                raise
              te_counts = map(operator.add,te_counts,new_te_multi_counts)
              
            else :
              te_counts = map(operator.add,te_counts,te_multi_counts)


    except:
            sys.stderr.write("Error occurred when assigning read gene/TE annotations in file %s.\n" % (filename))
            raise
    st = sum(te_counts)
    sg = sum(gene_counts.values())
    num_reads = st + sg

    sys.stderr.write("TE counts total %s\n" % (st))
    sys.stderr.write("Gene counts total %s\n" % (sg))
    sys.stderr.write("\nIn library %s:\n" % (filename))
    sys.stderr.write("Total annotated reads = %s\n" % ( str(num_reads) ) )
    sys.stderr.write("Total non-uniquely mapped reads = %s\n" % (str(int(nonunique))))
    sys.stderr.write("Total unannotated reads = %s\n\n" %(str(int(empty))))
    
    return (gene_counts,te_counts)
  
    
def parse_annotations_gene(annot_gene, gene_counts,leftOver_gene):
        no_annot = True

        if len(annot_gene) > 1 :
            leftOver_gene.append((annot_gene,1.0))
        elif len(annot_gene) == 1 :
            genes = annot_gene[0]
            if len(genes) == 1 :
                gene_counts[genes[0]] += 1
            else :
                if genes[0] == genes[1] :
                    gene_counts[genes[0]] += 1
                else :
                    gene_counts[genes[0]] += 0.5
                    gene_counts[genes[1]] += 0.5
        else :
            return True

        return False

# Assign ambiguous genic reads mapped to a location with multiple annotations
def resolve_annotation_ambiguity(counts, leftOvers ) :
    for (annlist,w) in leftOvers :
        readslist = {}
        total = 0.0
        size = len(annlist)
        ww = 1.0 * w
        if size > 1:
            ww = ww/size

        for ann in annlist :
            for a in ann:
                if a not in readslist :
                    readslist[a] = counts[a]
                    total += counts[a]

        if total > 0.0 :
            for a in readslist :
                v = ww * readslist[a] / total
                counts[a] += v
        else :
            for a in readslist :
                counts[a] = ww/len(readslist)


def parse_annotations_TE(multi_reads,annot_TE,annot_TE_position,Mappability,uniq_counts,multi_counts, readMappability,leftOver_list):
    if len(annot_TE) == 0 :
        return True

    if len(annot_TE) == 1 and len(annot_TE[0]) == 1 :
        uniq_counts[annot_TE[0][0]] += 1

    if len(annot_TE) == 1 and len(annot_TE[0]) > 1 :
        leftOver_list.append((annot_TE,1.0))

    if len(annot_TE) > 1:
        multi_algn = []
        for i in range(len(annot_TE)):
            for te in annot_TE[i] :
                multi_counts[te] += 1.0 / (len(annot_TE) * len(annot_TE[i]))
                readMappability[te] +=1.0 / Mappability[0][0]
                if annot_TE_position[te] == 0:
                    annot_TE_position[te]=len(annot_TE)
                else :
                    annot_TE_position[te].append(len(annot_TE))
                multi_algn.append(te)                   
        multi_reads.append(multi_algn)
        #print len(annot_TE_position)

    return False

def output_count_tbl(t_tbl, c_tbl, fname):
    try:
        f = open(fname, 'w')
    except IOError:
        error("Cannot create report file %s !\n" % (fname))
        sys.exit(1)
    else:
        cnt_tbl = {}
        header = "gene/TE"
        keys = set([])
        for tsmp in t_tbl.keys():
            keys = keys.union(t_tbl[tsmp].keys())
            header +="\t"+tsmp+".T"
        for csmp in c_tbl.keys():
            keys = keys.union(c_tbl[csmp].keys())
            header +="\t"+csmp+".C"

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

        for csmp in c_tbl.keys():
            cnts = c_tbl[csmp]
            for k in keys:
                val = 0
                if k in cnts :
                   val = cnts[k]

                if cnt_tbl.has_key(k) :
                    cnt_tbl[k].append(int(val))
                else :
                    vallist = []
                    vallist.append(int(val))
                    cnt_tbl[gene] = vallist
        #output
        f.write(header + "\n")
        for gene in sorted(cnt_tbl.keys()) :
           vals = cnt_tbl[gene]
           vals_str = gene
           for i in range(len(vals)) :
              vals_str +="\t"+str(vals[i])
           f.write(vals_str + "\n")

        f.close()

    return
    
def main():

    args=read_input(prepare_parser())
    
    info = args.info
    warn = args.warn
    debug = args.debug
    error = args.error

    # Output arguments used for program
    info("\n" + args.argtxt + "\n")

    info("Processing GTF files ...\n")

    geneIdx = GeneFeatures(args.gtffile,args.stranded,"exon","gene_id")
    info("Done building gene index ......\n")

    
    teIdx = TEfeatures()
    cur_time = time.time()
    te_tmpfile = '.'+str(cur_time)+'.te.gtf'
    subprocess.call(['sortBed -i '+ args.tefile+ ' >'+ te_tmpfile],shell=True )
    info("\nBuilding TE index .......\n")
    teIdx.build(te_tmpfile)
    subprocess.call(['rm -f ' + te_tmpfile ],shell=True)
    info("Done building TE index ......\n")
   
    # Read sample files make count table
    info("\nReading sample files ...\n")
    (tsamples_tbl, tsamples_rpm) = count_reads(args.tfiles, args.parser, geneIdx, teIdx, args.stranded, args.sortByPos, args.prj_name, numItr, fragLength, maxL)
    (csamples_tbl, csamples_rpm) = count_reads(args.cfiles, args.parser, geneIdx, teIdx, args.stranded, args.sortByPos, args.prj_name, numItr, fragLength, maxL)
    
    info("Finished processing sample files")

    info("Generating counts table")

    f_cnt_tbl = args.prj_name + ".cntTable"
    output_count_tbl(tsamples_tbl, csamples_tbl, f_cnt_tbl)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt !\n")
        sys.exit(0) 
    
