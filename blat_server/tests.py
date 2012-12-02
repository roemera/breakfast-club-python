from blat_server import *

def main():
    blatter = BlatClient('localhost',8123,False)
    dna = "AATCGAATCGAATCAAATCAAATCAAATCGAATAGGAAAGAATCGAATTGAATCAAGTAGAATCAAATCAAAAAGAATCGAATCGAGTTGAAAAAATTGAATCGAATTTAATCGAAAAGAATTGAATCAAAAATAATCGAATCGAATCGAAAGGAATCGAACCAAATCAAATAGAAAAGAATCGATTCGAATCGAATCAAAATAATCGAATTGAGTCGAATAAAAAAATCGAATCAAACGAATCGAATCGAAAGAATCGAATCAAAAAGAATTGAAACGAATAGAATCAAGTCGAATCGAAAAGAATCGAATGTAATCGAATGGAAAAGAATCGAATCAAAGAGAATTGA"
    aligns = blatter.blat_sequence(dna)
    for a in aligns:
        print a.chromosome,a.start,a.stop,a.strand
        print "Base Matches",a.base_matches
        print "Base mismatches",a.mismatches
        print "Percent Match",a.pct_match

        print "Percent Identitty",a.get_pct_ident()


if __name__ == '__main__':
    main()