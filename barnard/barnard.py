from rpy2.robjects.packages import importr
from rpy2.robjects import r

class Barnard:
	def __init__(self):
		self.barnard = importr("Barnard")
	def test(self, a, b, c, d):
		print self.barnard.barnardw_test(a,b,c,d, dp = 0.001)
		p_values = self.barnard.barnardw_test(a,b,c,d, dp = 0.001)[5]
		one_sided, two_sided = p_values
		return round(one_sided,4)
