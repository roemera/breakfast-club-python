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

class AssociationTesting:
	"""
		Class for association of MSATs/INDELs:
	"""
	def __init__(self):
		self.r = robjects.r

	def single_maker_logistic_association(self, phenotype_list=[], genotype_list=[]):
		"""
			Computes single marker logistic regressionassociation for 
			lists of phenotypes and genotypes of equal length 
		"""
		# Make sure phenotype and genotype lists are same size
		if len(phenotype_list) != len(genotype_list):
			return None

		# Converting genos and phenos to R vectors
		phenotypes = FloatVector(phenotype_list)
		genotypes = StrVector(genotype_list)

		# Grabbing alleles to return later
		alleles = sorted(set(genotypes))
		n_alleles = len(set(genotypes))

		# Model fitting
		robjects.globalenv["phenotypes"] = phenotypes
		robjects.globalenv["genotypes"] = genotypes

		lm = self.r.glm("phenotypes ~ genotypes",family = "binomial")	
		p_value = self.r.pchisq(lm.rx2('null.deviance') - lm.rx2('deviance'), lm.rx2('df.null') - lm.rx2('df.residual'), lower.tail = FALSE)

		# Compiling dict to return association_dict[allele] = [p-value,odds ratio]
		#association_dict = {}
		#for i in range(n_alleles):
		#	# Don't want to report results of intercept field
		#	if i > 0:
		#		association_dict[alleles[i]] = [self.r.summary(lm).rx2('coefficients')[i + (3 * n_alleles)], math.exp(self.r.summary(lm).rx2('coefficients')[i])]

		return p_value

	def single_maker_allelic_association(self, phenotype_list=[], genotype_list=[]):
		"""
			Computes single marker logistic regressionassociation for 
			lists of phenotypes and genotypes of equal length 
		"""
		# Make sure phenotype and genotype lists are same size
		if len(phenotype_list) != len(genotype_list):
			return None

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

		# Running fast Fisher's algoritm if 
		p = pvalue(case_counts[allele[0]], control_counts[allele[0]], case_counts[allele[1]], control_counts[allele[1]]).two_tail	
		return p


if __name__=='__main__':
	# Runs some tests if run as standalone

	# Genos and Phenos
	phenos = [1,1,1,1,1,0,0,0,0,0,0,0,0]
	genos = ["131,131","131,131","131,131","131,131","131,131","131,132","131,132","132,132","131,133","132,132","132,132","131,132","132,132"]
	
	# Initialize analysis class
	tester = AssociationTesting()

	# Allelic Association testing
	p_value = tester.single_maker_allelic_association(phenos, genos)
	print "Allelic Association pvalue was",p_value

	# Logistic regression testing
	p_values = tester.single_maker_logistic_association(phenos, genos)
	print "Logistic Association pvalue was",p_values

