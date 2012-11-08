# Association testing
from fisher import pvalue
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import FloatVector,StrVector
import math
from itertools import chain
from collections import Counter
import numpy as np
import scipy as sp

def logistic_regression(phenotype_list=[],genotype_list=[]):
	"""

	"""
	print "Made it into function"
	# Grabbing R env object
	r = robjects.r

	stats = importr("stats")
	print "Loaded R instance"
	# Converting genos and phenos to R vectors
	phenotypes = FloatVector(phenotype_list)
	genotypes = StrVector(genotype_list)
	
	# Grabbing alleles to return later
	alleles = sorted(set(genotypes))
	n_alleles = len(set(genotypes))
	r.cat("Hello")
	print "About to load env vars"
	# Model fitting
	robjects.globalenv["phenotypes"] = phenotypes
	robjects.globalenv["genotypes"] = genotypes

	
	sums = r.sum(phenotypes)
	print sums

	print "About to fit model"
	lm = stats.glm("phenotypes ~ genotypes - 1",family = "binomial")
	

	print "Model fitted"
	# Compiling dict to return association_dict[allele] = [p-value,odds ratio]
	association_dict = {}
	for i in range(n_alleles):
		association_dict[alleles[i]] = [r.summary(lm).rx2('coefficients')[i + (3 * n_alleles)], math.exp(r.summary(lm).rx2('coefficients')[i])]

	return association_dict

def allelic_association(phenotype_list=[],genotype_list=[]):
	"""
		##Fisher P-Values
		#http://en.wikipedia.org/wiki/Fisher's_exact_test

		#			CaseA	ControlA
		#With V		a		b	
		#Without	c		d
	"""
	case_alleles = []
	control_alleles = []
	for i in range(len(phenotype_list)):
		loc_alleles = genotype_list[i].split(",")
		if phenotype_list[i] == 1:
			for a in loc_alleles:
				case_alleles.append(a)
		else:
			for a in loc_alleles:
				control_alleles.append(a)

	allele = list(set(chain(case_alleles,control_alleles)))
	
	# Implementing slow scipy chi-square if we have more than two allles
	#if len(allele) > 2:
	#	table = np.zeros(shape)
	
	case_counts = Counter(case_alleles)
	control_counts = Counter(control_alleles)
	p = pvalue(case_counts[allele[0]], control_counts[allele[0]], case_counts[allele[1]], control_counts[allele[1]]).two_tail	
	print p

if __name__=='__main__':
	# Genos and Phenos
	phenos = [1,1,1,1,1,0,0,0,0,0,0,0,0]
	genos = ["131,131","131,131","131,131","131,131","131,131","131,132","131,132","132,132","131,132","132,132","132,132","131,132","132,132"]
	
	# Function testing
	print "Testing function"
	assoc = logistic_regression(phenos,genos)
	for k in assoc.keys():
		print "Allele",k,"\t","P-value",assoc[k][0],"\t","Odds Ratio",assoc[k][1]

	# Allelic Association testing
	allelic_association(phenos,genos)
