#!/usr/bin/env python3
import argparse
import os
import pysam
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from multiprocessing import Pool
import re
import numpy as np
import pandas as pd
import itertools
# from sknetwork.data import from_edge_list
# from functools import reduce
# from scipy.spatial.distance import hamming

from utils import bam_functions as bf
from utils import graph_utils as graphing
from utils import methylBam as mB
"""
# static bam functions
import sys

# append the path of the current directory
sys.path.append(".")
# from .bam_functions import *

from methlotype.utils import methylBam as mB
from methlotype.utils import bam_functions as bf
from methlotype.utils import graph_utils as graphing
"""

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
    # print('read cov:', read_cov)

    # plot barchart of methylation / canonical C in region
    # read_cov.T.plot.bar().get_figure().savefig(sample_name + '.fqsT.bar.png')
    # plt.clf()

    return read_cov


def aggregated_methylation_per_region(input_arg_list):
    #bamPath, min_prob, min_reads,min_cpgs, sample_name, region, scoreGraph, heatmap, slop):

    # methyl probability filter
    ml_min_prob = input_arg_list[1] #min_prob
    ml_cutoff = ml_min_prob * 255
    min_reads=int(input_arg_list[2])
    min_cpgs = int(input_arg_list[3])
    region = input_arg_list[5].split("_")
    sample_name=input_arg_list[4]
    slop=input_arg_list[-1]

    # run methyl analyze per read using url
    methylbam = mB.MethylBam(input_arg_list)
        #bamPath, ml_cutoff, min_reads, min_cpgs, sample_name, region, heatmaps=heatmap, slop=slop)
    mod_base_in_read, rdf, read_strands = methylbam.get_methyl_per_read()

    # print('\n\n',region,'shape',rdf.shape)
    # frequency look up at each position
    read_cov = get_cov_at_each_pos(rdf, sample_name)
    # if any column has less than min_reads remove it
    # remove columns with less than min_reads reads supporting that methyl position
    read_cov = read_cov.loc[:, read_cov.sum(axis=0) > min_reads]
    mean_region_read_cov = read_cov.sum(axis=0).mean()
    if np.isnan(mean_region_read_cov):
        mean_region_read_cov = rdf.shape[0]
    # print('mean_region_read_cov',mean_region_read_cov)
    # also remove positions that only have 1 methylated read
    # TODO think if this is the right way to measue this in a low coverage areas
    # Could be a ratio of 1/ mean_region_read_cov
    read_cov = read_cov.loc[:, read_cov.loc[1,:]>1]
    # print('\nread_cov', read_cov.shape,'\n',read_cov)

    # make a heatmap after filtering, but it no longer looks readable.
    # bf.make_quick_hm(read_cov, sample_name+'filtered')

    if rdf.shape[1]==0 or read_cov.shape[1]==0:
        agg_pct_mean = agg_region_mean =region_std = "-1"
        read_depth = str(mean_region_read_cov)
        num_cpgs="0"
        pass_min_cpgs="False"
        numModReads = "-1"
        numCanonReads = "-1"
        return region[0], str(int(region[1])-slop), str(int(region[2])+slop), agg_pct_mean, agg_region_mean, region_std, read_depth, num_cpgs, pass_min_cpgs, numModReads, numCanonReads


    else:    
        
        # regional meth percentage as methyated measurments / 
        row_sums = read_cov.sum(axis=1)
        # print('row sums\n', row_sums)

        # regional meth percentage as mean of methyl percent at each cpg + strand position
        pct_meth = read_cov.T
        pct_meth['pct_methyl'] = pct_meth[1]/(pct_meth[0]+pct_meth[1]) * 100
        agg_pct_mean = str(round((pct_meth['pct_methyl'].sum())/(pct_meth.shape[0]) ,4))
        # print('pct meth\n', pct_meth,'\n')

        
        numModReads = str(int(row_sums[1]))
        numCanonReads = str(int(row_sums[0]))
        agg_region_mean = str(round((row_sums[1] *100) / (row_sums[0]+row_sums[1]) ,4))
        region_std = str(round(rdf.std().mean() ,4))
        read_depth = str(round(len(mod_base_in_read.keys()) ,4))
        num_cpgs = str(round(read_cov.shape[1] ,4))
        pass_min_cpgs = str(read_cov.shape[1]>=min_cpgs)

        # print("\nchr","start","end","agg_pct_mean", "agg_region_mean", "region_std", "read_depth", "num_cpgs", "pass_min_cpgs", "numModReads", "numCanonReads", sep="\t")
        # print(region[0], str(int(region[1])-slop), str(int(region[2])+slop),agg_pct_mean, agg_region_mean, region_std, 
        #     read_depth, num_cpgs, pass_min_cpgs, numModReads, numCanonReads, sep="\t")


        return region[0], str(int(region[1])-slop), str(int(region[2])+slop), agg_pct_mean, agg_region_mean, region_std, read_depth, num_cpgs, pass_min_cpgs, numModReads, numCanonReads



def main(bamPath, bedPath, region_name, min_prob, min_reads, sample_name, 
            region, scoreGraph, heatmap, slop, min_cpgs):

    

    # store arguments to be passed to agg_meth_per_region in Pool multithreading 
    task_args_pre = [bamPath, min_prob, min_reads, min_cpgs, sample_name]
    task_args_suff = [heatmap, slop]

    # store the regions as a list of tasks to be parallelized 
    regions=[]

    with open(bedPath, "r") as bedfile:
        for line in bedfile:
            bed_entry = line.strip().split("\t")
            # skip any header lines
            if (bed_entry[0] == "chr") | (line[0] == "#"):
                continue
            if len(bed_entry)>3:
                bed_entry=bed_entry[:3]
            regions.append(tuple(task_args_pre+["_".join(bed_entry)]+task_args_suff) ) 

    # set the number of processes run in advance?
    # num_processes = 1
    print('regions\n',regions[:10])

    # Create a Pool of processes
    with Pool() as pool:
        # Map the task function to the Pool of processes
        results = pool.map(aggregated_methylation_per_region, regions)

    print('results\n',results)
    # for region in regions:
    #     aggregated_methylation_per_region(bamPath, min_prob, min_reads, sample_name, region, scoreGraph, heatmap, slop, min_cpgs)

    outFile = sample_name+"_"+region_name+".bed"
    with open(outFile,'w') as out_bed:
        header = ["#chr","start","end","agg_pct_mean", "agg_region_mean", "region_std", "read_depth", "num_cpgs", "pass_min_cpgs", "numModReads", "numCanonReads",str(sample_name),"\n"]
        out_bed.write("\t".join(header))
        for entry in results:
            out_bed.write("\t".join(entry))
            out_bed.write('\n')

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
	python3.9 methlotype/regionalMethylation.py -b https://public.gi.ucsc.edu/~memeredith/methyl/beds/methylation_R9/HG002_R9_R10/R9_R10_HG002/HG002_R9_phased.mod.hg38.bam \
    -s HG002_R9 -d Hs_EPDnew_006_hg38.bed --heatmap 0 -e EPDnew06_promoters -p 10
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
        required=False,
        type=str,
        help="space separated list of region components: chromosome start end"
    )

    parser.add_argument(
        "-d",
        required=False,
        type=str,
        help="Input BED file of Regions to aggregate methylation over"
    )

    parser.add_argument(
        "-e", "--region_name",
        required=False,
        type=str,
        help="Name of regions being aggregated"
    )

    parser.add_argument(
        "-p",
        required=False,
        type=float,
        default=0.5,
        help="Minimum probability cutoff for the likelyhood that a position is methylated (ml tag) "
    )

    parser.add_argument(
        "-n","--min_cpgs",
        required=False,
        type=float,
        default=5,
        help="Minimum number of cpg calls to be called as a region"
    )

    parser.add_argument(
        "-c",
        required=False,
        type=int,
        default=0,
        help="Minimum read coverage of a methyl position to be included in analysis "
    )

    parser.add_argument(
        "-x","--extra_bps",
        required=False,
        type=int,
        default=0,
        help="Number of extra base pairs to add to the start and end position of an input region"
    )

    parser.add_argument(
        "-g","--score_graph",
        required=False,
        type=int,
        default=0,
        help="1 to run the scoring and graphing, default 0: to not run scoring and graphing"
    )

    parser.add_argument(
        "-m","--heatmap",
        required=False,
        type=int,
        default=0,
        help="1 to plot heatmap of read methylation in the region, default 0: to not plot "
    )

    args = parser.parse_args()

    main(bamPath=args.b, bedPath=args.d, region_name=args.region_name, min_prob=args.p, min_reads=args.c, sample_name=args.s, region=args.r,
         scoreGraph=args.score_graph, heatmap=args.heatmap, slop=args.extra_bps, min_cpgs=args.min_cpgs)
