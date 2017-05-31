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
    
def ReadLengthDistribution( estimatedReadLength ): #reads distribution
    return 1
    
def EMUpdate(meansIn, te_features,multi_reads,estimatedReadLength,readMappability,annot_TE_position): #readMappability :list of reads map position
    #M-step 2nd equation
    meansOut = [0] *len(meansIn)
    multi_counts = computeAbundances(te_features,meansIn,multi_reads,estimatedReadLength)    
    sys.stderr.write("total multi counts = "+ str(sum(multi_counts))+"\n")   
    #print len(meansIn),len(multi_reads),len(multi_counts),len(readMappability)
	#########################################
    for tid in range(len(annot_TE_position)) :
        for i in range(len(annot_TE_position[tid])) :
			tlen = te_features.getLength(tid)
			#sys.stderr.write("#####################"+str(tid)+"\t"+str(tlen)+"\n" )#################
			if tlen <0 :
				sys.stderr.write("Error in optimization: the TE does not exist!\n")
				raise
			effectiveLength = tlen - estimatedReadLength + 1
			#sys.stderr.write("#####################"+str(tid)+"\t"+str(effectiveLength)+"\t"+str(tlen)+"\t"+str(estimatedReadLength)+"\n" )#################
			if  readMappability[tid]> 0 and  effectiveLength >0 : 
				meansOut[tid]+=  multi_counts[tid]* annot_TE_position[tid][i] / readMappability[tid] /effectiveLength
			else :
			#   sys.stderr.write("effective length is less than read length\n")
				meansOut[tid] = 0.0
	#########################################
    '''
    for tid in range(len(annot_TE_position)) :
        for i in range(len(annot_TE_position[tid])) :
            if readMappability[tid]>0 :
                meansOut[tid]+=  multi_counts[tid]* annot_TE_position[tid][i] / readMappability[tid]
            else :
                meansOut[tid] = 0.0
    '''                      
    meansOut = normalizeMeans(meansOut)
    
    sys.stderr.write("after normalization toatl means = "+str(sum(meansOut))+"\n") ##################
   
    return meansOut

    
def expectedLogLikelihood_(means, counts, te_features) :
    sampProbs = [0.0]*len(means)
    
    for tid in range(len(counts)) : 
        relativeAbundance = means[tid]
    
    if len(means) > 0 :
        return logLikelihood( sampProbs) / len(means)
    else :
        return 0.0
  

def dotProd_(u, v) :
    dot = 0.0
    dot = sum(map(operator.mul,u,v))
    
    return dot
        

def EMestimate(te_features,multi_reads,multi_counts,numItr,estimatedReadLength,readMappability,annot_TE_position): #multi_counts =Ft
    #E-step 
    
    sys.stderr.write("multi-reads = %s " % (str(len(multi_reads))))
    means0 = []
    ############################################
    for tid in range(len(multi_counts)) :
        tlen = te_features.getLength(tid)
        if tlen < 0 :
            sys.stderr.write("Error in optimization: the TE does not exist!\n")
            raise
        effectiveLength = tlen  - estimatedReadLength + 1
        #if effectiveLength < 0 :
        #    effectiveLength = 1
        if effectiveLength > 0 :
            means0.append(1.0 * readMappability[tid]/effectiveLength)
        else :
            #sys.stderr.write("effective length is less that read length mean0.\n")
            means0.append(0.0)
    
    # relative abundance    
    means0 = normalizeMeans(means0)
    ###############################################
    '''
    for tid in range(len(multi_counts)) :
        means0.append(1.0 * readMappability[tid])

    means0 = normalizeMeans(means0)
    '''
    ####################
    sys.stderr.write("after normalization toatl means0 = "+str(sum(means0))+"\n")
    '''
    /**
         * Defaults for these values taken from the R implementation of
         * [SQUAREM](http://cran.r-project.org/web/packages/SQUAREM/index.html).
         */
    '''
    minStep0 = 1.0
    minStep = 1.0
    maxStep0 = 1.0
    maxStep = 1.0
    mStep = 4.0
    nonMonotonicity = 1.0
    negLogLikelihoodOld = float("inf")
    negLogLikelihoodNew = float("inf")
    cur_iter = 0
    t_size = len(multi_counts)
    r = [0]*t_size
    r2 = [0] * t_size
    v = [0]*t_size
    meansPrime = [0.0] * t_size
    outerIteration = 1
    while cur_iter < numItr :
        cur_iter += 1
        sys.stderr.write("SQUAREM iteraton [" + str(cur_iter) + "]\n")
        sys.stderr.write("1/3\n")
        try :
            means1 = EMUpdate(means0, te_features,multi_reads,estimatedReadLength,readMappability,annot_TE_position)
        except :
            sys.stderr.write("Error in EMupdate\n")
            raise
        
        if negLogLikelihoodOld != float("inf") :
            negLogLikelihoodOld = -expectedLogLikelihood_(means0)
        
        sys.stderr.write("2/3\n")
        try:
            means2 = EMUpdate(means1, te_features,multi_reads,estimatedReadLength,readMappability,annot_TE_position)
        except :
            sys.stderr.write("Error in EMupdate\n")
            raise
        
        
        for tid in range(len(means0)) :
            r[tid] = means1[tid] - means0[tid]
            r2[tid] = means2[tid] - means1[tid]
            v[tid] = (means2[tid] - means1[tid]) - r[tid]
   
        rNorm = math.sqrt(dotProd_(r,r))
        r2Norm = math.sqrt(dotProd_(r2,r2)) 
        vNorm = math.sqrt(dotProd_(v,v))
        rr = dotProd_(r,v)
        rvNorm = math.sqrt(abs(rr))

        if vNorm == 0 :
            means0 = means1
            sys.stderr.write("at iteration " + str(cur_iter) + " vNorm == 0 \n")
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
            sys.stderr.write("alpha = " + str(alphaS) + ".\n ")
            sys.stderr.write("Performing a stabilization step.\n")
            try :
                meansPrime = EMUpdate(meansPrime, te_features,multi_reads,estimatedReadLength,readMappability,annot_TE_position)
            except :
                sys.stderr.write("Error in EMupdate\n")
                raise
    
                if alphaS == maxStep :
                    maxStep = max(maxStep0, maxStep/mStep)
                    alphaS = 1.0
                
        sys.stderr.write("alpha = " + str(alphaS) + ", ")
        
        if alphaS == maxStep :
            maxStep = mStep * maxStep
        
        if minStep < 0 and alphaS == minStep :
            minStep = mStep * minStep
        
        meansPrime, means0 = means0, meansPrime
        
        if not math.isnan(negLogLikelihoodNew) :
            negLogLikelihoodOld = negLogLikelihoodNew
    
    if cur_iter >= numItr :
        sys.stderr.write("not converge.....\n")
    else :
        sys.stderr.write("converge at iteration " + str(cur_iter)+"\n")    
    new_multi_counts = computeAbundances(te_features,means0,multi_reads,estimatedReadLength)
    return new_multi_counts

def computeAbundances(te_features,meansIn,multi_reads,estimatedReadLength): #estimatedReadLength is a constant,estimate from read length
    #M-step 1st equation   
    size = len(meansIn)
    multi_counts = [0] * size
    
    sys.stderr.write("num of multi reads = "+str(len(multi_reads))+"\n")
    ############################
    for kid in range(len(multi_reads)) :
        TE_transcripts = multi_reads[kid]
        totalMass = 0.0
        for tid in TE_transcripts :
              totalMass += meansIn[tid]
        
        if totalMass > 0.0 :
                  norm = 1.0 / totalMass
        else :
                  norm = 0.0
                  
        for tid in TE_transcripts :
              multi_counts[tid] += meansIn[tid] * norm
    #####################################
    '''
    for kid in range(len(multi_reads)) : #every read
        TE_transcripts = multi_reads[kid] #a set of TE copy overlap with every read
        totalMass = 0.0
        tlen = te_features.getLength(kid)
        if tlen <0 :
            sys.stderr.write("Error in optimization: the TE does not exist!\n")
            raise
            
        estimatedReadLengthP = ReadLengthDistribution (estimatedReadLength)
        effectiveLength = (tlen - estimatedReadLength + 1)*estimatedReadLengthP
        
        for tid in TE_transcripts :
              totalMass += meansIn[tid]*effectiveLength
        
        if totalMass > 0.0 :
                  norm = 1.0 / totalMass
        else :
                  norm = 0.0
                        
        for tid in TE_transcripts :
              multi_counts[tid] += meansIn[tid] * effectiveLength * norm
    '''
    sys.stderr.write("total multi counts = "+ str(sum(multi_counts))+"\n")
    return multi_counts
