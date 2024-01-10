import unittest
import time
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
    def test_methylGraph(self):
        mb = MethylBam('tiny.merged.2.bam', 255 * .05, 0, 'tiny.merged.2.bam', ['chr1', '33', '70'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()

        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'tiny.merged.2.bam')
        print(edges)
        graphing.make_graph(edges, read_strands, 'tiny.merged.2')


    def test_methylGraph6(self):
        mb = MethylBam('tiny.merged.6.bam', 255 * .05, 0, 'tiny.merged.6.bam', ['chr1', '33', '70'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'tiny.merged.6.bam')
        print(edges)
        graphing.make_graph(edges, read_strands, 'tiny.merged.6')

    def methylGraph_tiny3(self):
        mb = MethylBam('tiny.3.bam', 255 * .05, 0, 'tiny.3', ['chr20', '1000360', '1001560'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'tiny.3')
        print(edges)
        graphing.make_graph(edges, read_strands, 'tiny.3')

    def methylGraph_tiny3(self):
        mb = MethylBam('tiny.4.bam', 255 * .05, 0, 'tiny.4', ['chr20', '1000360', '1001560'])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'tiny.4')
        print(edges)
        graphing.make_graph(edges, read_strands, 'tiny.4')

    def methylGraph_snrpn(self):
        st = time.time()
        hg02 = 'http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam'
        mb = MethylBam(hg02, 255 * .95, 0, 'hg02.snrpn', ['chr15', 24955941, 24955999])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        mb_time = time.time()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'hg02.snrpn')
        read_scores_time = time.time()
        # print(edges)
        graphing.make_graph(edges, read_strands, 'hg02.snrpn')
        graph_time = time.time()
        print('mb_time:', mb_time - st, '\neverlap_edgewts time:', read_scores_time - mb_time,
              'graph_time', graph_time - read_scores_time, '\ntotal time:', graph_time - st)

    def methylGraph_snrpn_2k(self):
        st = time.time()
        hg02 = 'http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam'
        mb = MethylBam(hg02, 255 * .95, 0, 'hg02.snrpn_2k.', ['chr15', 24954682, 24956837])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        mb_time = time.time()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'hg02.snrpn_2k')
        read_scores_time = time.time()
        # print(edges)
        graphing.make_graph(edges, read_strands, 'hg02.snrpn_2k')
        graph_time = time.time()
        print('_2k\n mb_time:', mb_time - st, '\neverlap_edgewts time:', read_scores_time - mb_time,
              'graph_time', graph_time - read_scores_time, '\ntotal time:', graph_time - st)

    # chr15:24955689-24955918
    def methylGraph_snrpn_2c(self):
        st = time.time()
        hg02 = 'http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam'
        mb = MethylBam(hg02, 255 * .95, 0, 'hg02.snrpn_2c.', ['chr15', 24955689, 24955918])
        mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()
        mb_time = time.time()
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, 'hg02.snrpn_2c')
        read_scores_time = time.time()
        # print(edges)
        graphing.make_graph(edges, read_strands, 'hg02.snrpn_2c')
        graph_time = time.time()
        print('_2k\n mb_time:', mb_time - st, '\neverlap_edgewts time:', read_scores_time - mb_time,
              'graph_time', graph_time - read_scores_time, '\ntotal time:', graph_time - st)

if __name__ == '__main__':


    unittest.main()



