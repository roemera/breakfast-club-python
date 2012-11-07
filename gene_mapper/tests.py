import unittest
from parsers import GenotypeReportReader

class GenotypeReportReaderTests(unittest.TestCase):
    def setUp(self):
        self.reader = GenotypeReportReader(open("data/TEAL_by_DOBB_WOBB_25_1_1_CG Lineup Report.txt"))
    def test(self):
        #Right number of genotypes
        self.failUnless(len(self.reader) == 672 - 14)
        #Check first genotype
        first_gtype = self.reader[0]
        self.failIf(first_gtype.sample_alias == "BLANK")
        self.failUnless(first_gtype.sample_alias == "BS CK193")
        self.failUnless(first_gtype.marker_name == "INDEL 5785")
        self.failUnless(first_gtype.row == "A")
        self.failUnless(first_gtype.col == 1)
        self.failUnless(first_gtype.alleles == [138, 138])
        #Check other genotypes
        self.failUnless(self.reader[3].alleles == [134, 138])
        self.failUnless(len(self.reader["INDEL 5785"]) == 94)

def main():
    unittest.main()

if __name__ == '__main__':
    main()