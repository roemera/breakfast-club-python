# Association testing
from fisher import pvalue
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import FloatVector,StrVector
import math
from itertools import chain
from collections import Counter
import numpy as np
from scipy import stats


r = robjects.r
r_stats = importr("stats")

def logistic_regression(phenotype_list=[],genotype_list=[]):
	# Converting genos and phenos to R vectors
	phenotypes = FloatVector(phenotype_list)
	genotypes = StrVector(genotype_list)
	
	# Grabbing alleles to return later
	alleles = sorted(set(genotypes))
	n_alleles = len(set(genotypes))

	# Model fitting
	robjects.globalenv["phenotypes"] = phenotypes
	robjects.globalenv["genotypes"] = genotypes

	lm = r_stats.glm("phenotypes ~ genotypes - 1",family = "binomial")	

	# Compiling dict to return association_dict[allele] = [p-value,odds ratio]
	association_dict = {}
	for i in range(n_alleles):
		association_dict[alleles[i]] = [r.summary(lm).rx2('coefficients')[i + (3 * n_alleles)], math.exp(r.summary(lm).rx2('coefficients')[i])]

	return association_dict



def allelic_association(phenotype_list=[],genotype_list=[]):
	"""
		
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

	# Getting set of alleles and their counts
	allele = list(set(chain(case_alleles,control_alleles)))
	case_counts = Counter(case_alleles)
	control_counts = Counter(control_alleles)

	# Implementing slow scipy chi-square if we have more than two allles
	if len(allele) > 2:
		table = np.zeros(shape=(2,len(allele)))
		for i in range(len(allele)):
			table[0,i] = case_counts[allele[i]]
			table[1,i] = control_counts[allele[i]]

		chi2, p, dof, ex = stats.chi2_contingency(table)
		return p


	p = pvalue(case_counts[allele[0]], control_counts[allele[0]], case_counts[allele[1]], control_counts[allele[1]]).two_tail	
	return p

if __name__=='__main__':
	# Genos and Phenos
	phenos = [1,1,1,1,1,0,0,0,0,0,0,0,0]
	genos = ["131,131","131,131","131,131","131,131","131,131","131,132","131,132","132,132","131,133","132,132","132,132","131,132","132,132"]
	
	# Function testing
	print "Testing function"
	#assoc = logistic_regression(phenos,genos)
	#for k in assoc.keys():
	#	print "Allele",k,"\t","P-value",assoc[k][0],"\t","Odds Ratio",assoc[k][1]

	# Allelic Association testing
	p_value = allelic_association(phenos,genos)
	print "Allelic Association pvalue was",p_value