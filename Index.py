#!/usr/bin/env python2.7

import pdb
import sys, time, re
import logging
import operator
import collections
from math import ceil,floor

class TEfeatures :
    
    def __init__ (self):

        self.__namelist = collections.defaultdict(dict)
        self.__window = collections.defaultdict(dict)
        self.indexlist = {}
        self._length = []
        self._nameIDmap = []
        self._elements = []
        self._start = []
        self._end= []

    def getNames(self) :
        names = []
        return self._nameIDmap
    def numInstances(self) :
        return len(self._nameIDmap)

    def getElements(self) :
         return self._elements

    def getStrand(self,idx) :
        f_name = self._nameIDmap[idx]
        return f_name[len(f_name)-1]

    def getEleName(self,idx) :
        full_name = None
        if idx >= len(self._nameIDmap) or idx < 0 :
            return None
        else :
            full_name =  self._nameIDmap[idx]
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

    def getLength(self,TE_name_idx) :
        if TE_name_idx < len(self._length) :
            return self._length[TE_name_idx]
        else :
            return -1
                    
    def findOvpTE(self,chrom,start,end,Mappability): #attention when test
        TEnamelist = []
        readMappability = []
        for keys in self.__window[chrom].keys() :
            [(s,e)] = self.__window[chrom][keys]                        
            if start <= e and end >= s :
                #print s,e
                try :
                    for i in range(499):
                        [(s,e)] = self.__namelist[chrom][keys+i]
                        if start <= e and end >= s :
                            TEnamelist.append(keys)
                            break   
                except KeyError :
                    pass
                    
                break
        if len(TEnamelist)>0:
            readMappability.append(Mappability)
          
        return (TEnamelist,readMappability)
    
    def TE_annotation(self,iv_seq): #iv_seq store reads information
        TEs = []
        readMappability =[]
        for iv in iv_seq :
            chromo = iv[0]
            start = iv[1]
            end = iv[2]
            strand = iv[3]
            Mappability = iv[4]
            (name_idx_list, Mappability) = self.findOvpTE(chromo,start,end,Mappability)

            if name_idx_list is not None :
                for t in name_idx_list :
                    if strand != "." : #stranded
                        if strand == self.getStrand(t)  :
                            if t not in TEs :
                                TEs.append(t)
                                readMappability=Mappability                             
                    else :#not stranded
                        if t not in TEs :
                            TEs.append(t)
                            readMappability=Mappability              
        
        return (TEs,readMappability)

    #group by element
    def groupByEle(self,te_inst_counts) :

        TEs = self.getElements()
        te_ele_counts = dict(zip(TEs,[0]*len(TEs)))

        for i in range(len(te_inst_counts)) :
            ele_name = self.getEleName(i)

            if ele_name is None:
                sys.stderr.write("TE out of index boundary!\n")
                sys.exit(1)

            if ele_name in te_ele_counts :
                te_ele_counts[ele_name] += te_inst_counts[i]
            else :
                sys.stderr.write("TE inconsistency! "+ele_name+"\n")
                sys.exit(1)

        return te_ele_counts

    def build (self,filename):
            self.__srcfile = filename             
            try:
                f = open(self.__srcfile,'r')
            except:
                logging.error("cannot open such file %s !\n" %(self.__srcfile))
                sys.exit(1)

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

                if ele_id == "" or name == "" or family_id == "" or class_id == "" :
                    sys.stderr.write(line+"\n")
                    sys.stderr.write("TE GTF format error! There is no annotation at line %s.\n" % (linenum))
                    raise

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
                
         

            

#    def __init__ (self):
#
#        self._qname = {}
#        self._flag = []
#        self._rname = []
#        self._pos = []
#       self._mapq = []
#       self._cigar = []
#       self._rnext = []
#       self._pnext =[]
#       self._seq = []
#       self._nh =[]
        


class GFF_Reader( ):


   def __init__( self, filename, id_attribute):
      self.line_no = None
      self.filename = filename
      self.id_attribute = id_attribute
      self._re_attr_main = re.compile( "\s*([^\s\=]+)[\s=]+(.*)" )

   def __iter__( self ):
      self.line_no = 0
      if self.filename.lower().endswith( ( ".gz" , ".gzip" ) ):
            lines = gzip.open( self.filename )
      else:
            lines = open( self.filename )

      for line in lines:
          self.line_no += 1
          if line == "\n" or line.startswith('#'):
              continue
          ( chr, source, feature, start, end, score, strand, frame, attributeStr ) = line.split("\t")
          id = self.__parse_GFF_attr_string(attributeStr,self.id_attribute)#return transcript_id

          yield (id, chr, strand, int(start), int(end), feature)

      lines.close()
      self.line_no = None

   def __parse_GFF_attr_string(self,attributeStr,id_interested) :


       for pairs in attributeStr.split(';') :
           if pairs.count('"') not in [0,2] :
               raise ValueError, "The attribute string seems to contain mismatched quotes."
           nv = self._re_attr_main.match(pairs)
           if not nv :
               raise ValueError, "Failure parsing GFF attribute line."
           val = nv.group(2)
           name = nv.group(1)
           if name == id_interested :
               return val
       return None

   def get_line_number_string( self ):
      if self.line_no is None:
            return "file %s closed" % self.filename

      else:
         return "line %d of file %s" % ( self.line_no, self.filename )


class GeneFeatures:
    """index of Gene annotations.
        """
    def __init__ (self,GTFfilename,stranded,feature_type,id_attribute):

        #dict of dicts since the builtin type doesn't support it for some reason
        self.featureIdxs_plus = collections.defaultdict(dict)
        self.featureIdxs_minus = collections.defaultdict(dict)
        self.featureIdxs_nostrand = collections.defaultdict(dict)
        self.features = []

        self.read_features(GTFfilename,stranded,feature_type,id_attribute)


    # Reading & processing annotation files
    def read_features(self,gff_filename, stranded, feature_type, id_attribute) :

        # read count of features in GTF file
        gff = GFF_Reader(gff_filename,id_attribute)  # (id, chr, strand, int(start), int(end), feature)recode gene feature
        i = 0
        counts = 0
        try:
            for f in gff:
                if f[0] is None :
                    continue
                if f[5] == feature_type:
                    counts += 1
                    if stranded != "no" and f[2] == "." :
                        sys.stderr.write("Feature %s does not have strand information." % (f[0]))
                    try:
                        if f[2] == "."  :
                            self.featureIdxs_nostrand[f[1]][f[0]].append((f[3],f[4]))
                    except:
                        self.featureIdxs_nostrand[f[1]][f[0]] = [(f[3],f[4])]

                    try:
                        if f[2] == "+"  :
                            self.featureIdxs_plus[f[1]][f[0]].append((f[3],f[4]))
                    except:
                        self.featureIdxs_plus[f[1]][f[0]] = [(f[3],f[4])]

                    try:
                        if f[2] == "-" :
                            self.featureIdxs_minus[f[1]][f[0]].append((f[3],f[4]))
                    except KeyError:
                        self.featureIdxs_minus[f[1]][f[0]] = [(f[3],f[4])]

                    #save gene id
                    if f[0] not in self.features :
                        self.features.append(f[0])

                    i += 1
                    if i % 100000 == 0 :
                        sys.stderr.write("%d GTF lines processed.\n" % i)
        except:
            sys.stderr.write("Error occured in %s.\n" % gff.get_line_number_string())
            raise

        if counts == 0 :
            sys.stderr.write("Warning: No features of type '%s' found in gene GTF file.\n" % feature_type)
            
    def getFeatures(self) :
        return self.features
    
    def findOvpgene(self,start,end): #attention when test
        gene_namelist=[]
        for each_gene in self.keys() :
            for (s,e) in self[each_gene]:
                if start <= e and end >= s :                    
                    gene_namelist.extend(each_gene)
                                            
            gene_namelist +=gene_namelist
        
        return gene_namelist
                                           

    def Gene_annotation(self,itv_list): #itv_list = exon_bound.append([chrom, chrom_st,chrom_st + s-1,"."])
        genes = []
        fs = []
        for itv in itv_list :
            try:
                if itv[3] == "+" :
                    if itv[0] in self.featureIdxs_plus :
                        fs = self.featureIdxs_plus[itv[0]].findOvpgene(itv[1],itv[2])


                if itv[3] == "-" :
                        if itv[0] in self.featureIdxs_minus:
                            fs = self.featureIdxs_minus[itv[0]].findOvpgene(itv[1], itv[2])


                if itv[3] == "." :
                        if itv[0] in self.featureIdxs_minus:
                            fs = self.featureIdxs_minus[itv[0]].findOvpgene(itv[1], itv[2])

                        if itv[0] in self.featureIdxs_plus :
                            fs += self.featureIdxs_plus[itv[0]].findOvpgene(itv[1],itv[2])
                        if itv[0] in self.featureIdxs_nostrand :
                            fs += self.featureIdxs_nostrand[itv[0]].findOvpgene(itv[1],itv[2])

                if len(fs) > 0:
                        genes = genes + fs

            except AttributeError:
                    pass


        return genes

