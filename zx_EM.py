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
	sys.stderr.write("after normalization toatl means = "+str(sum(meansOut))+"\n") ##################
   
	return meansOut

	
def EM(TE,multi_reads,uniq_counts,te_multi_counts,numItr,avgReadLength,multi_counts,readMappability,TE_position_list,OPT_TOL=0.0001):
	#initializaiton
	means0 = []
	for tid in range(len(uniq_counts)) :
			means0.append(multi_counts[tid])
	means0 = normalizeMeans(means0) 
	#################
	minStep0 = 1.0
	minStep = 1.0
	maxStep0 = 1.0
	maxStep = 1.0
	mStep = 4.0
	cur_iter = 0
	t_size = len(uniq_counts)
	r = [0]*t_size
	r2 = [0] * t_size
	v = [0]*t_size
	meansPrime = [0.0] * t_size
	outerIteration = 1
	while cur_iter < numItr :
		cur_iter += 1
		
		means1 = M_step(means0,TE,uniq_counts,multi_reads,avgReadLength,readMappability,TE_position_list)					  
		means2 = M_step(means1,TE,uniq_counts,multi_reads,avgReadLength,readMappability,TE_position_list)
		
		for tid in range(len(means0)) :
			r[tid] = means1[tid] - means0[tid]
			r2[tid] = means2[tid] - means1[tid]
			v[tid] = (means2[tid] - means1[tid]) - r[tid]
			
		rNorm = math.sqrt(sum(map(operator.mul,r,r)))
		r2Norm = math.sqrt(sum(map(operator.mul,r2,r2))) 
		vNorm = math.sqrt(sum(map(operator.mul,v,v)))
		rr = sum(map(operator.mul,r,v))
		rvNorm = math.sqrt(abs(rr))
		
		if vNorm == 0 :
			means0 = means1
			break
		alphaS = rNorm / rvNorm
		alphaS = max(minStep, min(maxStep, alphaS))
		
		if rNorm < OPT_TOL :
			sys.stderr.write("rNome = OPT_TOL \n")
			break
		if r2Norm < OPT_TOL :
			sys.stderr.write("r2Nome = OPT_TOL \n")
			means0 = means2
			break
			
		for tid in range(len(means0)) :
			meansPrime[tid] = max(0.0, means0[tid] + 2*alphaS*r[tid] + (alphaS*alphaS)*v[tid])	
		
		# Stabilization step
		if abs(alphaS - 1.0) > 0.01 :
			try :
				meansPrime = M_step(meansPrime,TE,uniq_counts,multi_reads,avgReadLength,readMappability,TE_position_list)
			except :
				sys.stderr.write("Error in EMupdate\n")
				raise		
			
			if alphaS == maxStep :
				maxStep = max(maxStep0, maxStep/mStep)
				alphaS = 1.0
		
		if alphaS == maxStep :
			maxStep = mStep * maxStep
		
		if minStep < 0 and alphaS == minStep :
			minStep = mStep * minStep
		
		meansPrime, means0 = means0, meansPrime
			
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
			totalMass += meansIn[te]*effectiveLength
		
		if totalMass > 0.0 :
				  norm = 1.0 / totalMass
		else :
				  norm = 0.0
				  
		for te in multi_reads[f] :
			  multi_counts[te] += meansIn[te] * norm * TE_position_list[f]/readMappability[f]
	
	sys.stderr.write("total multi counts = "+ str(sum(multi_counts))+"\n")
	return multi_counts
		
		