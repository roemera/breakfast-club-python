# Association testing
from fisher import pvalue
import rpy2.robjects as robjects
from rpy2.robjects import FloatVector,StrVector
import math

def logistic_regression(phenotype_list=[],genotype_list=[]):
	"""

	"""
	# Grabbing R env object
	r = robjects.r

	# Converting genos and phenos to R vectors
	phenotypes = FloatVector(phenotype_list)
	genotypes = StrVector(genotype_list)
	
	# Grabbing alleles to return later
	alleles = sorted(set(genotypes))
	n_alleles = len(set(genotypes))

	# Model fitting
	robjects.globalenv["phenotypes"] = phenotypes
	robjects.globalenv["genotypes"] = genotypes
	lm = r.glm("phenotypes ~ genotypes - 1",family = "binomial")
	
	# Compiling dict to return association_dict[allele] = [p-value,odds ratio]
	association_dict = {}
	for i in range(n_alleles):
		association_dict[alleles[i]] = [r.summary(lm).rx2('coefficients')[i + (3 * n_alleles)], math.exp(r.summary(lm).rx2('coefficients')[i])]

	return association_dict

"""
def allelic_association(case_alleles=[],control_alleles=[]):
	a=0
	b=
	alleles = list(set(unlist([case_alleles,control_alleles])))
	p = pvalue(a, b, c, d).two_tail	
	return p
"""

if __name__=='__main__':
	# Genos and Phenos
	phenos = [1,1,1,1,1,0,0,0,0,0,0,0,0]
	genos = ["131,131","131,131","131,131","131,131","131,131","131,132","131,132","132,132","131,132","132,132","132,132","131,132","132,132"]
	
	# Function testing
	print "Testing function"
	assoc = logistic_regression(phenos,genos)
	for k in assoc.keys():
		print "Allele",k,"\t","P-value",assoc[k][0],"\t","Odds Ratio",assoc[k][1]

	"""
	# Allelic Association testing
	cases = []
	controls = []
	for i in range(len(phenos)):
		alleles = genos[i].split(",")
		if phenos[i] == 1:
			for a in alleles:
				cases.append(a)
		else:
			for a in alleles:
				controls.append(a)

	p_val = allelic_association()
	"""