# python modules

import math 
import sys
import array
import operator

from constants import *

def normalizeMeans(meansIn):
	total = sum(meansIn)
	meansOut = [0]*len(meansIn)
	sys.stderr.write("total means = " + str(total) +"\n")	 
	if total > 0 :
		meansOut = map(lambda x: 1.0*x/total, meansIn)
	return meansOut
	
def M_step(meansIn,TE,uniq_counts,multi_reads,avgReadLength,readMappability,TE_position_list):
	meansOut = [0] *len(meansIn)
	multi_counts=E_step(meansIn,TE,multi_reads,avgReadLength,readMappability,TE_position_list)
	for tid in range(len(meansIn)) :
		meansOut[tid] = multi_counts[tid]						  
	meansOut = normalizeMeans(meansOut)
	#sys.stderr.write("after normalization toatl means = "+str(sum(meansOut))+"\n") ##################
   
	return meansOut

	
def EM(TE,multi_reads,uniq_counts,te_multi_counts,numItr,avgReadLength,multi_counts,readMappability,TE_position_list,OPT_TOL=0.0001):
	#initializaiton
	means0 = []
	for tid in range(len(uniq_counts)) :
			means0.append(multi_counts[tid])
	means0 = normalizeMeans(means0) 
	#################
	cur_iter = 0
	t_size = len(uniq_counts)
	r = [0]*t_size
	outerIteration = 1
	while cur_iter < numItr :
		cur_iter += 1
		
		means1 = M_step(means0,TE,uniq_counts,multi_reads,avgReadLength,readMappability,TE_position_list)					  
		
		for tid in range(len(means0)) :
			r[tid] = means1[tid] - means0[tid]
			
		rNorm = math.sqrt(sum(map(operator.mul,r,r)))
		
		if rNorm < OPT_TOL :
			means1 = means0
			#sys.stderr.write("rNome = OPT_TOL \n")
			break
		
		means1, means0 = means0, means1
			
	if cur_iter >= numItr :
		sys.stderr.write("not converge.....\n")
	else :
		sys.stderr.write("converge at iteration " + str(cur_iter)+"\n")	   
	new_multi_counts = E_step(means0,TE,multi_reads,avgReadLength,readMappability,TE_position_list)
	return new_multi_counts
			
def E_step(meansIn,TE,multi_reads,avgReadLength,readMappability,TE_position_list):
	
	multi_counts = [0] *len(meansIn)
	
	sys.stderr.write("num of multi reads = "+str(len(multi_reads))+"\n")
	
	for f in range(len(multi_reads)) :
		totalMass = 0.0		
		for te in multi_reads[f] :
			tlen = TE.getLength(te)
			effectiveLength = (tlen - avgReadLength + 1)
			totalMass += meansIn[te]*effectiveLength    #E_step
		
		if totalMass > 0.0 :
				  norm = 1.0 / totalMass
		else :
				  norm = 0.0
				  
		for te in multi_reads[f] :
			  multi_counts[te] += meansIn[te] * norm * TE_position_list[f]/readMappability[f]    #M_step
	
	#sys.stderr.write("total multi counts = "+ str(sum(multi_counts))+"\n")
	return multi_counts
		
		