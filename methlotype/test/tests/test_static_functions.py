import unittest
import pandas as pd
#import TestCase

import sys
# append the path of the
# parent directory
sys.path.append("../..")

# from methlotype.analyze_methylBam import *
from methlotype.utils.bam_functions import *


class Test(unittest.TestCase):
    def test_rev_comp(self):
        print('\ntest_rev_comp')
        nucs = ['A', 'C', 'G', 'T']
        rev_nucs = ['T', 'G', 'C', 'A']
        for i, base in enumerate(nucs):
            self.assertEqual(rev_comp(base), rev_nucs[i])

    def test_identify_methyl_tags(self):
        print('\ntest_identify_methyl_tags')
        bamfile = pysam.AlignmentFile('tiny.2.bam', "rb")
        for read in bamfile.fetch('chr20', 1000360, 1000860):
            mm, ml = identify_methyl_tags(read)
            self.assertEqual(mm, 'Mm')
            self.assertEqual(ml, 'Ml')

    def test_identify_methyl_tags_uppercase(self):
        print('\ntest_identify_methyl_tags_uppercase')
        bamfile = pysam.AlignmentFile('tiny.bam', "rb")
        for read in bamfile.fetch('chr1', 33, 50):
            mm, ml = identify_methyl_tags(read)
            self.assertEqual(mm, 'MM')
            self.assertEqual(ml, 'ML')

    def test_get_hp_tag(self):
        print('\ntest_get_hp_tag')
        bamfile = pysam.AlignmentFile('tiny.2.bam', "rb")
        for read in bamfile.fetch('chr20', 1000360, 1000860):
            self.assertEqual(get_hp_tag(read), 1)


    def test_get_cigar_ref_pos(self):
        """test cigar string ref position parser.
        need sample read with cigar string, make one like Jean does :)"""
        print('\ntest_get_cigar_ref_pos')
        ref_positions = [0, 0, 0, 0, 0, 32, 33, 34, 35, 36, 38, 39, 40, 41, 42,
                         43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
                         57, 58, 58, 58, 58, 58, 59, 60]
        bamfile = pysam.AlignmentFile('tiny.bam', "rb")
        for read in bamfile.fetch('chr1', 33, 50):
            print(read.query_name, read.is_reverse)
            print(get_cigar_ref_pos(read))
            print(read.get_reference_positions(full_length=True))
            self.assertEqual(get_cigar_ref_pos(read), ref_positions)

    def test_get_cigar_ref_pos_rev(self):
        """test cigar string ref position parser."""
        print('\ntest_get_cigar_ref_pos rev read alignment')
        ref_positions = [33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45,
                         46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                         58, 59, 60, 61, 62, 63, 64, 64, 64, 64, 64, 65, 66]
        bamfile = pysam.AlignmentFile('tiny.rev.bam', "rb") # tiny.bam
        for read in bamfile.fetch('chr1', 33, 50):
            print(read.query_name, read.is_reverse)
            print(list(read.query_sequence))
            print(get_cigar_ref_pos(read))
            print(read.get_reference_positions(full_length=True))
            self.assertEqual(get_cigar_ref_pos(read), ref_positions)




if __name__ == '__main__':


    unittest.main()