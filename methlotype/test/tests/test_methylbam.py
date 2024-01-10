import unittest
import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal
#import TestCase

import sys
# append the path of the
# parent directory
sys.path.append("../..")

from methlotype.analyze_methylBam import *
from methlotype.utils.methylBam import *

class Test(unittest.TestCase):
    # def test_MethylBam_class(self):
    def methylBam(self):
        print('test_MethylBam_class')

        # run methyl analyze per read using url
        mb = MethylBam('tiny.3.bam', 255 * .05, 0, 'test', ['chr20', '1000360', '1000860'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()

        truth_mod_base_in_read = {'e58e0b04-626e-46d6-80b8-fb5c2fd682fc':
                                      {'chr20': {1000438: [1000359, 1001305, 1, 253, 106, '+'],
                                                 1000813: [1000359, 1001305, 3, 150, 477, '+'],
                                                 1000923: [1000359, 1001305, 4, 131, 587, '+'],
                                                 1000987: [1000359, 1001305, 6, 127, 650, '+'],
                                                 1001185: [1000359, 1001305, 9, 62, 850, '+'],
                                                 1001206: [1000359, 1001305, 10, 255, 871, '+'],
                                                 1001264: [1000359, 1001305, 11, 255, 929, '+'],
                                                 1001285: [1000359, 1001305, 12, 198, 950, '+']}},
                                  'e58e0b04-626e-46d6-80b8-fb5c2fd682fx':
                                      {'chr20': {1000438: [1000359, 1001305, 1, 253, 106, '+'],
                                                 1000813: [1000359, 1001305, 3, 150, 477, '+'],
                                                 1000923: [1000359, 1001305, 4, 131, 587, '+'],
                                                 1000987: [1000359, 1001305, 6, 127, 650, '+'],
                                                 1001185: [1000359, 1001305, 9, 62, 850, '+'],
                                                 1001206: [1000359, 1001305, 10, 255, 871, '+'],
                                                 1001264: [1000359, 1001305, 11, 255, 929, '+'],
                                                 1001285: [1000359, 1001305, 12, 198, 950, '+']}},
                                  '83fb086e-7e17-4041-ac70-f795c1fd286a':
                                      {'chr20': {1000812: [1000559, 1001130, 3, 39, 274, '+'],
                                                 1000949: [1000559, 1001130, 4, 16, 410, '+']}}}

        truth_rdf = pd.DataFrame(1, index=['e58e0b04-626e-46d6-80b8-fb5c2fd682fc',
                                           'e58e0b04-626e-46d6-80b8-fb5c2fd682fx',
                                           '83fb086e-7e17-4041-ac70-f795c1fd286a'],
                                 columns=[1000438, 1000812, 1000813])
        truth_rdf.at['e58e0b04-626e-46d6-80b8-fb5c2fd682fc', 1000812] = 0
        truth_rdf.at['e58e0b04-626e-46d6-80b8-fb5c2fd682fx', 1000812] = 0
        truth_rdf.at['83fb086e-7e17-4041-ac70-f795c1fd286a', 1000438] = -1
        truth_rdf.at['83fb086e-7e17-4041-ac70-f795c1fd286a', 1000813] = 0

        truth_read_strands = {'e58e0b04-626e-46d6-80b8-fb5c2fd682fc': ['+', 1, 1000359, 1001305],
                              'e58e0b04-626e-46d6-80b8-fb5c2fd682fx': ['+', 1, 1000359, 1001305],
                              '83fb086e-7e17-4041-ac70-f795c1fd286a': ['+', 0, 1000559, 1001130]}

        self.assertEqual(mod_base_in_read, truth_mod_base_in_read, 'mod_base_in_read failed, tiny.3.bam')
        self.assertEqual(truth_rdf.to_dict(), rdf.to_dict(), 'rdf failed, tiny.3.bam')
        self.assertEqual(read_strands, truth_read_strands, 'read_strands failed, tiny.3.bam')

    def test_methylBam2(self):
        print('test_MethylBam_class2')

        # run methyl analyze per read using url
        mb = MethylBam('tiny.merged.2.bam', 255 * .05, 0, 'test', ['chr1', '33', '70'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()

        truth_mod_base_in_read = {'read_28833_29006_6954':
                                      {'chr1': {34: [32, 60, 0, 231, 6, '+'],
                                                43: [32, 60, 1, 200, 14, '+'],
                                                45: [32, 60, 2, 146, 16, '+'],
                                                60: [32, 60, 3, 220, 35, '+']}},
                                  'read_28833_29006_6954_2':
                                      {'chr1': {34: [32, 60, 0, 231, 6, '+'],
                                                43: [32, 60, 1, 200, 14, '+'],
                                                45: [32, 60, 2, 146, 16, '+'],
                                                60: [32, 60, 3, 220, 35, '+']}},
                                  'read_28833_29006_6954_r': {
                                      'chr1': {60: [33, 61, 0, 231, 36, '-'],
                                               45: [33, 61, 1, 200, 17, '-'],
                                               43: [33, 61, 2, 146, 15, '-'],
                                               34: [33, 61, 3, 220, 7, '-']}},
                                  'read_28833_29006_6945':
                                      {'chr1': {63: [33, 66, 0, 231, 30, '-'],
                                                55: [33, 66, 1, 200, 22, '-'],
                                                53: [33, 66, 2, 146, 20, '-'],
                                                33: [33, 66, 3, 220, 1, '-']}}}
        truth_rdf = pd.DataFrame(1, index=['read_28833_29006_6945',
                                           'read_28833_29006_6954_2',
                                           'read_28833_29006_6954_r',
                                           'read_28833_29006_6954'],
                                 columns=[33, 34, 43, 45, 53, 55, 60, 63])
        for read in ['read_28833_29006_6954', 'read_28833_29006_6954_2', 'read_28833_29006_6954_r']:
            for pos in [33, 53, 55]:
                truth_rdf.at[read, pos] = 0
            truth_rdf.at[read, 63] = -1
        for pos in [34, 43, 45, 60]:
            truth_rdf.at['read_28833_29006_6945', pos] = 0


        truth_read_strands = {'read_28833_29006_6954': ['+', 0, 32, 60],
                              'read_28833_29006_6954_2': ['+', 0, 32, 60],
                              'read_28833_29006_6954_r': ['-', 0, 33, 61],
                              'read_28833_29006_6945': ['-', 1, 33, 66]}
        # print('modbase in read', mod_base_in_read)
        print('rdf', rdf)
        # print(truth_rdf)
        # print('read strands', read_strands)
        self.assertEqual(mod_base_in_read, truth_mod_base_in_read, 'mod_base_in_read failed, tiny.merged.2.bam')
        self.assertEqual(truth_rdf.to_dict(), rdf.to_dict(), 'rdf failed, tiny.merged.2.bam')
        self.assertEqual(read_strands, truth_read_strands, 'read_strands failed, tiny.merged.2.bam')


    # def test_overlap_function(self):
         # print('test_overlap_function.tiny3')

         # mb = MethylBam('tiny.3.bam', 255 * .05, 0, 'test', ['chr20', '1000360', '1000860'])
         # mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
         # read_scores = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'test')
         # print('tiny3: mod_base', mod_base_in_read)
         # print('rdf',  rdf)
         # print('read strands', read_strands)
         # print('read scores', read_scores)

    def test_overlap_function2(self):
        print('test_overlap_function.merged')

        mb = MethylBam('tiny.merged.2.bam', 255 * .05, 0, 'tiny.merged.2.bam', ['chr1', '33', '70'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'tiny.merged.2.bam')

        # print('rdf', rdf)
        # print('read scores', read_scores.to_dict())

        # make truth data
        read_hap_labels = ['read_28833_29006_6954_+_0', 'read_28833_29006_6954_2_+_0',
                           'read_28833_29006_6954_r_-_0 ', 'read_28833_29006_6945_-_1']
        match1 = 1.423593519226589 #1.7251552301180073
        match2 = 9.31276870800169  #5.145182434739834
        mismatch = 0 #-11.101234942578339
        truth_read_scores = pd.DataFrame(np.nan, index=read_hap_labels, columns=read_hap_labels)
        truth_read_scores.at[read_hap_labels[0], read_hap_labels[0]] = match1
        truth_read_scores.at[read_hap_labels[0], read_hap_labels[1]] = match1
        truth_read_scores.at[read_hap_labels[0], read_hap_labels[2]] = match1
        truth_read_scores.at[read_hap_labels[0], read_hap_labels[3]] = mismatch
        truth_read_scores.at[read_hap_labels[1], read_hap_labels[1]] = match1
        truth_read_scores.at[read_hap_labels[1], read_hap_labels[2]] = match1
        truth_read_scores.at[read_hap_labels[1], read_hap_labels[3]] = mismatch
        truth_read_scores.at[read_hap_labels[2], read_hap_labels[2]] = match1
        truth_read_scores.at[read_hap_labels[2], read_hap_labels[3]] = mismatch
        truth_read_scores.at[read_hap_labels[3], read_hap_labels[3]] = match2

        print('truth', truth_read_scores)
        message = '3 bam read scores not identical, expect to be ran with replacement'
        self.assertTrue( truth_read_scores.equals(read_scores), message)



    def overlap_function3(self):
        print('test_overlap_function3: 6 bams 2sets of matching reads, 1 different')

        mb = MethylBam('tiny.merged.6.bam', 255 * .05, 0, 'tiny.merged.6.bam', ['chr1', '33', '70'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'tiny.merged.6.bam')

        truth_read_scores = {'read_28833_29006_6954_+_0':
                                 {'read_28833_29006_6954_+_0': 5.015067828019353,
                                  'read_28833_29006_6954_2_+_0': np.nan,
                                  'read_28833_29006_6945_-_1': np.nan,
                                  'read_28833_29006_6946_-_1': np.nan,
                                  'read_28833_29006_6946_2_-_1': np.nan},
                             'read_28833_29006_6954_2_+_0':
                                 {'read_28833_29006_6954_+_0': 5.015067828019353,
                                  'read_28833_29006_6954_2_+_0': 5.015067828019353,
                                  'read_28833_29006_6945_-_1': np.nan,
                                  'read_28833_29006_6946_-_1': np.nan,
                                  'read_28833_29006_6946_2_-_1': np.nan},
                             'read_28833_29006_6945_-_1':
                                 {'read_28833_29006_6954_+_0': -7.4552398885714455,
                                  'read_28833_29006_6954_2_+_0': -7.4552398885714455,
                                  'read_28833_29006_6945_-_1': 9.418602572271237,
                                  'read_28833_29006_6946_-_1': np.nan,
                                  'read_28833_29006_6946_2_-_1': np.nan},
                             'read_28833_29006_6946_-_1':
                                 {'read_28833_29006_6954_+_0': -10.917387670126178,
                                  'read_28833_29006_6954_2_+_0': -10.917387670126178,
                                  'read_28833_29006_6945_-_1': -12.42912907328786,
                                  'read_28833_29006_6946_-_1': 6.4321581405801025,
                                  'read_28833_29006_6946_2_-_1': np.nan},
                             'read_28833_29006_6946_2_-_1':
                                 {'read_28833_29006_6954_+_0': -10.917387670126178,
                                  'read_28833_29006_6954_2_+_0': -10.917387670126178,
                                  'read_28833_29006_6945_-_1': -12.42912907328786,
                                  'read_28833_29006_6946_-_1': 6.4321581405801025,
                                  'read_28833_29006_6946_2_-_1': 6.4321581405801025}}

        truth_read_scores_df = pd.DataFrame.from_dict(truth_read_scores)
        print('truth', truth_read_scores_df)
        message = '6 bam read scores not identical, expect to be ran with replacement; tiny.merged.6.bam'
        self.assertTrue( truth_read_scores_df.equals(read_scores), message+'_assertTrue' )



if __name__ == '__main__':


    unittest.main()