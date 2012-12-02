# Python lib for interacting with BLAT server
import os,sys
import math

"""
	BLAT server info:
		gfServer start localhost 8123 canFam3.2bit &
		gfServer start localhost 8124 hg19.2bit &
		gfServer start localhost 8125 mm10.2bit &
		gfServer -trans -mask start localhost 8126 canFam3.2bit &
		gfServer -trans -mask start localhost 8127 hg19.2bit &
		gfServer -trans -mask start localhost 8128 mm10.2bit &
"""

class BlatClient:
	def __init__(self, host_name, port_number, is_protein_bool):
		self.host = host_name
		self.port = port_number
		self.is_protein = is_protein_bool

	def blat_sequence(self, nuc_string=""):
		"""
			Take a nucleotide sequence. Returns a list of BlatAlignment objects. 
		"""	
		blat_alignments = []

		to_blat = ">blatter\n%s" % nuc_string
		print to_blat
		print "echo '%s' | gfClient %s %i /blat_server/ stdin stdout" % (to_blat ,self.host, self.port)
		blat_results = os.popen("echo '%s' | gfClient %s %i /blat_server/ stdin stdout" % (to_blat ,self.host, self.port)).read().split("\n")
		for result in blat_results:
			spl = result.split("\t")
			if spl[0].isdigit():			
				algn = BlatAlignment(result,self.is_protein)
				blat_alignments.append(algn)

		return blat_alignments

class BlatAlignment:
	"""
		PSL format line and returns a BlatAlignment object
	"""
	def __init__(self, psl_line, is_protein_bool):
		matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts = psl_line.split("\t")
		
		self.start = int(tStart.strip())
		self.stop = int(tEnd.strip())
		self.chromosome = tName.strip()
		self.strand = strand.strip()
		self.base_matches = int(matches.strip())
		self.mismatches = int(misMatches.strip())
		self.repeat_matches = int(repMatches.strip())
		self.query_size = int(qSize.strip())
		self.pct_match = float(self.base_matches) / float(self.query_size)
		self.span = self.stop - self.start
		self.num_gaps = int(tNumInsert.strip())
		self.total_gap_size = int(tBaseInsert.strip()) 
		self.n_masked_bases = int(nCount.strip())
		self.query_start = int(qStart.strip())
		self.query_stop = int(qEnd.strip())
		self.query_gap_count = int(qNumInsert.strip())
		self.query_gap_bases = int(qBaseInsert.strip())
		self.block_count = int(blockCount.strip())
		self.block_sizes = [int(x.strip()) for x in blockSizes.strip().split(",") if x != ""]
		self.query_block_starts = [int(x.strip()) for x in qStarts.strip().split(",") if x != ""]
		self.block_starts = [int(x.strip()) for x in tStarts.strip().split(",") if x != ""]
		self.is_protein = is_protein_bool

	def get_pct_ident(self):
		"""
			Calculates UCSCs percent identity score
		"""
		if self.is_protein:
			sizeMul = 3
		else:
			sizeMul = 1

		# Calculate alignment size
		qAliSize = sizeMul * (self.query_stop - self.query_start)
		tAliSize = self.stop - self.start
		aliSize = min(qAliSize, tAliSize)

		# Break if something fucked up
		if aliSize <= 0:
			return 0

		# Get difference in alignment size
		sizeDif = qAliSize - tAliSize
		if sizeDif < 0:
			sizeDif = -sizeDif
		
		insertFactor = self.query_gap_count

		total = (sizeMul * (self.base_matches + self.repeat_matches + self.mismatches))
		if total != 0:
			milliBad = (1000 * (self.mismatches * sizeMul + insertFactor + round(3 * math.log10(1 + sizeDif)))) / total
		
		pct_ident = 100.0 - milliBad * 0.1
		
		return pct_ident


