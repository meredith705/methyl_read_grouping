#!/usr/bin/env python3
import argparse
import os
import pysam
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd
import itertools
from sknetwork.data import from_edge_list
from functools import reduce
from scipy.spatial.distance import hamming

# from utils import bam_functions as bf
# from utils import graph_utils as graphing
# from utils import methylBam as MethylBam

# static bam functions
import sys

# append the path of the current directory
sys.path.append(".")
# from .bam_functions import *

from methlotype.utils import methylBam as mB
from methlotype.utils import bam_functions as bf
from methlotype.utils import graph_utils as graphing


def get_cov_at_each_pos(df, sample_name):
    """ store the frequency of each methylation character at each position for easy look up """

    # positional_frequencies = pd.DataFrame(np.nan, index=[-1, 0, 1], columns=list(df.columns))
    read_cov = pd.DataFrame(np.nan, index=[-1, 0, 1], columns=list(df.columns))
    # get a count for each character at each position
    for pos in df.columns:
        vdf = df.loc[:, pos].value_counts()
        # added: divide by # of total values to make it a fq
        for c in vdf.index:
            # positional_frequencies.at[c, pos] = vdf[c]#/vdf.sum()
            read_cov.at[c, pos] = vdf[c]

    # print(positional_frequencies)
    print('read cov:', read_cov)

    read_cov.T.plot.bar().get_figure().savefig(sample_name + '.fqsT.bar.png')
    plt.clf()

    return read_cov


def make_overlap_edge_weights(reads_dict, rdf, read_strands, sample_name):
    ''' find reads that overlap methyl positions and score the overlap '''
    # consider dataframe of positions, rows are reads, where alignments are already lined up
    # calculate fq of C, M at each position in overlap between the two reads

    scores = []
    edges = []
    read_hap_labels = ['_'.join([key, str(items[0]), str(items[1])]) for key, items in read_strands.items()]
    read_scores = pd.DataFrame(np.nan, index=read_hap_labels, columns=read_hap_labels)

    # Conditional Probability Distribution
    cpd = pd.DataFrame({'C': [0.975, 0.025],
                        'M': [0.025, 0.975]},
                       index=[0, 1])
    # print('cpd', cpd)

    # dynamic proramming speeds up frequency look up at each position
    read_cov = get_cov_at_each_pos(rdf, sample_name)

    # iterate over all pairs of reads ( with replacement? )
    for a, b in itertools.combinations_with_replacement(reads_dict, 2):
    # for a, b in itertools.combinations(reads_dict, 2):
        if bf.get_chrom(reads_dict[a]) == bf.get_chrom(reads_dict[b]):
            a_b_positional_values = []
            a_read_alignment = sorted(list(reads_dict[a][bf.get_chrom(reads_dict[a])].keys()))
            b_read_alignment = sorted(list(reads_dict[b][bf.get_chrom(reads_dict[b])].keys()))
            # print('a', a, a_read_alignment)
            # print('b', b, b_read_alignment)

            # if b overlaps a; calculate probabilities
            if (a_read_alignment[0] <= b_read_alignment[0]) and (a_read_alignment[-1] >= b_read_alignment[0]):

                overlap_start = bf.find_overlap_index(a_read_alignment, b_read_alignment)
                overlap_end = bf.find_overlap_end_index(a_read_alignment, b_read_alignment)
                # print('overlap start, end: ', overlap_start, overlap_end)
                # make sure the end of the overlap doesn't go past the end of b
                if overlap_end:
                    if (overlap_end - overlap_start) > len(b_read_alignment):
                        overlap_end = overlap_start + len(b_read_alignment) - 1

                # isolate the dataframe to only where these 2 reads overlap
                read_pair_df = rdf.loc[[a, b], a_read_alignment[overlap_start]:a_read_alignment[overlap_end]]
                # print('read_pair_df', read_pair_df)
                # iterate through all positions that overlap
                for i, aPos in enumerate(read_pair_df.columns):
                    if aPos in rdf.columns:
                        # store the methyl character for each read at the current position
                        a_char = rdf.loc[a, aPos]
                        b_char = rdf.loc[b, aPos]

                        # read coverage: total number of reads with 1 or 0 at this position
                        read_cov_at_pos = read_cov.loc[[0, 1], aPos].sum()
                        # print(i,aPos, 'read_cov_at_pos', read_cov_at_pos)

                        # random overlap ( denominator )
                        # frequency of the symbol at position i in read a and read b.
                        # read coverage / total number of reads at that position
                        a_char_freq = read_cov.loc[a_char, aPos] / read_cov_at_pos # qxi
                        b_char_freq = read_cov.loc[b_char, aPos] / read_cov_at_pos # qyi
                        # print('a_char_fq', a_char_freq, 'b', b_char_freq)

                        # cdp probabilities
                        # If there are 0(cannonical C) or 1(methyl C) characters calculate their frequency
                        # for c in cpd.columns:
                        freq_c_pos_i = 0        #qci = 0
                        freq_methyl_c_pos_i = 0 #qmi = 0
                        # if there is a value there, calculate the frequency of 0(C) and 1(5mC) at this position
                        if not np.isnan(read_cov.loc[0, aPos]):
                            freq_c_pos_i = read_cov.loc[0, aPos] / read_cov_at_pos  # rdf.shape[0]
                        if not np.isnan(read_cov.loc[1, aPos]):
                            freq_methyl_c_pos_i = read_cov.loc[1, aPos] / read_cov_at_pos  # rdf.shape[0]
                        # print('freq_c_pos_i', freq_c_pos_i, 'freq_methyl_c_pos_i', freq_methyl_c_pos_i )
                        # multiply the frequency by the cpd value for that character for both reads
                        if a_char != -1 and b_char != -1:
                            # product of
                            prob_xi_given_chari = freq_c_pos_i * cpd.loc[a_char, 'C'] * cpd.loc[b_char, 'C']
                            prob_yi_given_chari = freq_methyl_c_pos_i * cpd.loc[a_char, 'M'] * cpd.loc[b_char, 'M']
                            # then add them together
                            prob_overlap_c_methyl_c = prob_xi_given_chari + prob_yi_given_chari
                            prob_overlap_by_chance = a_char_freq * b_char_freq
                            # store the log of the ratio
                            a_b_positional_values.append(np.log(prob_overlap_c_methyl_c / prob_overlap_by_chance))
                            # print('prob_xi_given_chari',  prob_xi_given_chari,   'prob_yi_given_chari', prob_yi_given_chari)
                            # print('prob_overlap_c_methyl_c', prob_overlap_c_methyl_c,
                            #       'prob_overlap_by_chance', prob_overlap_by_chance)
                            # print(i, aPos, 'log', np.log(prob_overlap_c_methyl_c / prob_overlap_by_chance))
                    else:
                        print('not in df?', aPos)

        score = sum(a_b_positional_values)

        a_label = "_".join([a, str(read_strands[a][0]), str(read_strands[a][1])])
        b_label = "_".join([b, str(read_strands[b][0]), str(read_strands[b][1])])
        # print(a_label, b_label, 'score:', score, '\n\n')
        read_scores.at[a_label, b_label] = score
        scores.append(score)

        # append edge to list of edges, if it is a matching score
        if score > 0:# and a_label != b_label:
            edges.append( (a_label, b_label, score) )
        plt.scatter(a, score)

    plt.savefig(sample_name + "score.scatter.png")
    plt.clf()

    print(sample_name, read_scores)
    # print('read score mode\n', read_scores.mode(axis=1))
    read_scores.to_csv('read_scores_dataframe.' + sample_name + '.csv')
    # make_quick_hm(read_scores, sample_name + '.scores.region')
    # make_pca(read_scores.fillna(0), sample_name + '.region.scores')
    # ax = read_scores.plot.hist(legend=False, colormap='Blues').get_figure().savefig(sample_name + '.scores.hist.png')
    num_bins = int(max(scores) - min(scores) / 2)
    plt.hist(scores, bins = num_bins)
    plt.title(sample_name+" Scores Histogram")
    plt.savefig(sample_name + ".plt.scores_hist.png")
    plt.clf()

    return read_scores, edges




def main(bamPath, min_prob, min_reads, sample_name, region, scoreGraph):

    # methyl probability filter
    ml_min_prob = min_prob
    ml_cutoff = ml_min_prob * 255

    # run methyl analyze per read using url
    methylbam = mB.MethylBam(bamPath, ml_cutoff, min_reads, sample_name, region)
    mod_base_in_read, rdf, read_strands = methylbam.get_methyl_per_read()

    # plot_modbase_locations(mod_base_in_read)

    # make_overlap_edges(mod_base_in_read, 10,read_strands)
    if scoreGraph:
        print('making edge weights')
        read_scores, edges = make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, sample_name)
        print('graphing')
        graphing.make_graph( edges, read_strands, sample_name)


if __name__ == "__main__":
    '''
	Python script to cluster reads in a bam by their methyl tag provile

	Usage: 
	python3.9 analyze_methylBam.py -b http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam
	-s sample_name -r chr20 63664503 63678780
	
	# 'chr15', 24953000, 24958133):#'chr20',63664503,63678780):#'chr2',51026934,51027176):#
	    # neuron loc chr20:63664503-63664830 chr20:63664503-63678780  'chr20',63664503,63678780):
	    chr20:58815000-58895000
        # neuron loc #1: chr2:51026934-51027176
        # SNRPN chr15:24953000-24958133 'chr15', 24953000, 24958133):     'chr15',24955941,24955999):#    'chr15', 24954860,24956812):#
	python3.9 ../../../methlotype/analyze_methylBam.py -b http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam -s HG002.snrpn -r chr15 24953000 24958133
	
    https://public.gi.ucsc.edu/~memeredith/methyl/beds/methylation_R9/HG002_R9_R10/R9_R10_HG002/HG002_R9_phased.mod.hg38.bam

    '''

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-b",
        required=True,
        type=str,
        help="Input BAM, can be file path or http URL "
    )

    parser.add_argument(
        "-s",
        required=True,
        type=str,
        help="sample name"
    )

    parser.add_argument(
        "-r",
        nargs=3,
        required=True,
        type=str,
        help="space separated list of region components: chromosome start end"
    )

    parser.add_argument(
        "-p",
        required=False,
        type=float,
        default=0.5,
        help="Minimum probability cutoff for the likelyhood that a position is methylated (ml tag) "
    )

    parser.add_argument(
        "-c",
        required=False,
        type=int,
        default=0,
        help="Minimum read coverage of a methyl position to be included in analysis "
    )

    parser.add_argument(
        "-g","--score_graph",
        required=False,
        type=int,
        default=1,
        help="Boolean to run the scoring and graphing "
    )

    args = parser.parse_args()

    print('score_graph',args.score_graph)

    main(bamPath=args.b, min_prob=args.p, min_reads=args.c, sample_name=args.s, region=args.r,
         scoreGraph=args.score_graph)
