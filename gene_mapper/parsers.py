from csv import DictReader

class GenotypeReportGenotype:
	def __init__(self, report_line_dict):
		#Use csv dict as init to fill genotype vars
		sample_alias,date,well_pos,_ = report_line_dict["Sample File"].rsplit("_",3)
		#Remove underscore from sample name if presnet
		self.sample_alias = sample_alias.replace("_", " ")
		#Marker name
		self.marker_name = report_line_dict["Marker"].replace("_", " ")
		#Row & Col
		self.row = well_pos[0]; self.col = int(well_pos[1:])
		#Allels
		self.alleles = []

		allele_1 = report_line_dict["Allele 1"]
		allele_2 = report_line_dict["Allele 2"]

		if allele_1.strip() != "?" and allele_1.strip() != "":
			self.alleles.append(int(allele_1))

			if allele_2.strip() == "" or allele_2.strip() == "?":
				self.alleles.append(int(allele_1))

			elif allele_2.strip() != "" and allele_2.strip() != "?":
				self.alleles.append(int(allele_2))


class GenotypeReportReader:
	def __init__(self, file_handle):
		#Use csv dict reader to read in report
		csv_reader = DictReader(file_handle, delimiter = "\t")
		self.genotypes = []
		for report_line_dict in list(csv_reader):
			if not report_line_dict["Sample File"].lower().startswith("blank"):
				self.genotypes.append(GenotypeReportGenotype(report_line_dict))

	def __len__(self):
		return len(self.genotypes)

	def __getitem__(self, key):
		if type(key) == int:
			return self.genotypes[key]
		elif type(key) == str:
			return self.get_genotype_by_marker_name(key)

	def get_genotype_by_marker_name(self, key):
		genotypes = []
		for genotype in self.genotypes:
			if genotype.marker_name == key:
				genotypes.append(genotype)
		return genotypes

	def get_marker_names(self):
		marker_names = []
		for genotype in self.genotypes:
			marker_names.append(genotype.marker_name)
		return list(set(marker_names))

	def get_marker_alleles(self, marker_name):
		alleles = []
		for genotype in self[marker_name]:
			if len(genotype.alleles) > 0:
				alleles.extend(",".join(genotype.alleles)
		return list(set(alleles))

	def get_genotype(self, marker_name = "", sample_alias = ""):
		for genotype in self.genotypes:
			if genotype.marker_name == marker_name and genotype.sample_alias == sample_alias:
				return genotype
		return None

