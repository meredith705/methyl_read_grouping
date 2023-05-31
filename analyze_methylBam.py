#!/usr/bin/env python3
import argparse
import os
import pysam
from Bio.Seq import Seq
import re
import numpy as np
import pandas as pd
import itertools
from functools import reduce
from scipy.spatial.distance import hamming

# heatmap 
import seaborn as sn
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# louvain
from scipy import sparse
from sknetwork.clustering import Louvain, get_modularity
from sknetwork.data import from_edge_list
from sknetwork.embedding import Spring
from sknetwork.linalg import normalize
from sknetwork.utils import get_membership
from sknetwork.visualization import svg_graph, svg_bigraph
from IPython.display import SVG


# some static functions for the MethylBam class
def identify_methyl_tags(read):
    """ the methyl tags can be either Mm/Ml or MM/ML for R10 or R9 """
    if read.has_tag('Mm') and read.has_tag('Ml'):
        mm_tag = 'Mm'
        ml_tag = 'Ml'
    elif read.has_tag('MM') and read.has_tag('ML'):
        mm_tag = 'MM'
        ml_tag = 'ML'
    else:
        print("Read", read.query_name, " does not have either Mm/Ml or MM/ML tags.")
        return -1, -1
    return mm_tag, ml_tag


def rev_comp(base):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return complement.get(base, base)


def get_hp_tag(read):
    if read.has_tag('HP'):
        return read.get_tag('HP')
    else:
        return 0


def get_cigar_ref_pos(read):
    """ parse the cigar string to obtain methylated cytosine positions
        in refernce coordinate space """
    pos = read.reference_start
    refcigarpositions = []

    # for anti-sense strand read alignments the cigar tuples need to be read from back to front,
    # they are created with respect to the reference.
    if read.is_reverse:
        cigartups = read.cigartuples  # [::-1]
        pos += 1
    else:
        cigartups = read.cigartuples

    for i, c in enumerate(cigartups):
        # Read consuming CIGAR operations
        # 4 : S : soft clip; set empty position holders because Mm tag includes these positions; can be H_S_query_S_H at ends
        # 5 : H : hard clip; clipped sequences NOT present in SEQ don't count these positions;
        # 			only present as the first and/or last operation.
        # not sure if H clipped sequences are included in the Mm tag; if so add more ref position placeholders.
        if c[0] == 4:
            if i < 2:
                refcigarpositions[0:0] = [0] * c[1]
            else:
                refcigarpositions.extend([0] * c[1])
        # 1 : INS; add place holder ref coordinate for each read-only base consumed by the insertion
        if c[0] in [1]:
            refcigarpositions.extend([pos] * c[1])

        # 0 : M :  match/mismatch add ref positions to the array for every M base
        # 7 : = :  add ref positions to the array for every = base
        # 8 : X :  mismatch add ref position as if it was =
        if c[0] in [0, 7, 8]:
            refcigarpositions.extend([i for i in np.arange(pos, pos + c[1])])
            pos += c[1]

        # 2 : D :  DEL;skip the ref position
        # 3 : N : REF_SKIP;
        if c[0] in [2, 3]:
            pos += c[1]

    return refcigarpositions


class MethylBam:
    """ Class for storing and parsing a BAM file that contains methyl tags from ONT sequencing """
    def __init__(self,bamPath, ml_cutoff, min_reads, sample_name, region):
        self.bamPath = bamPath
        self.ml_cutoff = ml_cutoff
        self.min_reads = min_reads
        self.sample_name = sample_name
        # this won't always be run on a region
        self.region = region
        self.mod_base_in_read = {}
        self.modbase = None
        self.mm = None


    def get_methyl_per_read(self):
        """ iterate over all the read mapping to a specified region using fetch().
        Each iteration returns a AlignedSegment object which represents a single read along with its fields and optional tags:
        """
        # open alignment file
        samfile = pysam.AlignmentFile(self.bamPath, "rb")

        logfile = open('log.txt', 'w')

        read_set = set()
        read_strands = {}
        chrom_modbase_pos = {}

        # Fetch reads from region: chr20:1000000-1020000
        mmTag, mlTag = '', ''
        if len(self.region) == 3:
            chromosome = self.region[0]
            start = int(self.region[1])
            end = int(self.region[2])
        else:
            return -1
        # func to store strand as +/- instead of boolean
        strand = lambda x: '-' if x else '+'

        # TODO This won't always be per region....
        # TODO maybe filter out short read alignments, of some cutoff
        for read in samfile.fetch(chromosome, start, end):
            # store some easy access meta data
            read_set.add(read.query_name)
            # read_strands[read.query_name]= [strand(read.is_reverse), self.get_hp_tag(read)]
            bfile = open('beds/' + str(read.query_name) + '.mods.bed', 'w')
            logfile.write(' '.join([ '\n', str(strand(read.is_reverse)), str(read.mapq), str(read.reference_name),
                                     str(read.reference_start), str(read.reference_end), str(read.query_name) ]) )

            # get mm and ml tag information per read
            # the methyl tags can be either Mm/Ml or MM/ML for R10 or R9
            if len(mmTag) == 0:
                mm_tag, ml_tag = identify_methyl_tags(read)
                if ml_tag == -1:
                    return -1

            # mm tag needs to be split and skip bases list isolted:
            self.parse_mm_tag(read, mm_tag)
            ml = read.get_tag(ml_tag)

            # store the seq in the proper orientation relative to the Ml tag
            seq = read.seq  # seq = identify_seq_orientation(read)

            # store the index position of every base of interest in the read
            mod_base_i = [b.start() for b in re.finditer(self.modbase.upper(), str(seq).upper())]

            # get reference positions that account for INDELs/Clips by parsing the cigar string
            refcigarpositions = get_cigar_ref_pos(read)

            # store ref start and end positions per read for easy access
            read_strands[read.query_name] = [strand(read.is_reverse), get_hp_tag(read),
                                             refcigarpositions[read.query_alignment_start],
                                             refcigarpositions[read.query_alignment_end - 1]]

            # for a REVERSE(-) alignment move through the alignment from back to front
            # counting rev complement C's (aka G) in increments reported by mm tag
            if read.is_reverse:
                # an index for where the next modified base is
                mod_base_i_index = 0
                for i, i_skip in enumerate(self.mm):
                    # increment the read position by the number of modbase indices to skip
                    # the base after skipped ones is the modified base, so add 1
                    mod_base_i_index += int(i_skip) + 1

                    # only count modified bases within the read.query_alignment_start and end positions
                    if read.query_alignment_start < mod_base_i[-mod_base_i_index] < read.query_alignment_end:
                        # sanity check that the base at this position is the modified base in the mm tag
                        self.check_modbase_rev(seq, mod_base_i, mod_base_i_index)
                        # if the ml probability of being modified is > cutoff,
                        if ml[i] >= self.ml_cutoff:
                            # the reference position of the mod base
                            # store this as pos - 1 for reverse strands, so the mod base is on the C of the + strand
                            pos = refcigarpositions[mod_base_i[-mod_base_i_index]] - 1

                            # write to bed:
                            bfile.write(str(read.reference_name) + '\t' + str(pos - 1) + '\t' + str(pos) + '\t' + str(
                                ml[i]) + '\t' + '-' + '\n')

                            self.fill_dictionary(read, chrom_modbase_pos, pos,
                                            [refcigarpositions[read.query_alignment_start],
                                             refcigarpositions[read.query_alignment_end - 1], i, ml[i],
                                             mod_base_i[-mod_base_i_index], '-'], logfile)

            # for a FORWARD(+) alignment move through the alignment by moving in the increments reported by the mm tag
            else:
                read_modbase_pos = 0
                for i, i_skip in enumerate(self.mm):
                    # it is the modbase after the # to skip
                    if i == 0:
                        read_modbase_pos = int(i_skip)
                    else:
                        # increment the read position by the number of modbases indicies to skip
                        read_modbase_pos = read_modbase_pos + int(i_skip) + 1

                    # only count modified bases after the read.query_alignment_start position
                    if read.query_alignment_start < mod_base_i[read_modbase_pos] < read.query_alignment_end:
                        self.check_modbase(seq, mod_base_i, read_modbase_pos)

                        # if the ml probability of being modified is > cutoff,
                        # store the modified base location as it 0-based position in the read sequence string
                        if ml[i] >= self.ml_cutoff:
                            pos = refcigarpositions[mod_base_i[read_modbase_pos]] + 1

                            # write to bed
                            bfile.write(str(read.reference_name) + '\t' + str(pos - 1) + '\t' + str(pos) + '\t' + str(
                                ml[i]) + '\t' + '+' + '\n')

                            self.fill_dictionary(read, chrom_modbase_pos, pos,
                                            [refcigarpositions[read.query_alignment_start],
                                             refcigarpositions[read.query_alignment_end - 1], i, ml[i],
                                             mod_base_i[read_modbase_pos], '+'], logfile)

            bfile.close()

        # for each chromosome make a new methylation df
        for chrom, position_dt in chrom_modbase_pos.items():
            print('chrom',chrom)
            # make an empty dataframe of the correct size; set this to all -1
            df = pd.DataFrame(-1, index=read_set, columns=sorted(list(position_dt.keys())))

            # using query start and end account for reads not overlapping a position to be counted
            # separate from unmethylated positions in an overlapping read
            # position not covered: -1   unmethylated position: 0   metylated position: 1
            for read, values in read_strands.items():
                df.at[read, values[2]:values[3]] = 0
            # store methylated positions per read
            for pos in sorted(list(position_dt.keys())):
                # for each methylated position, add 1 to the dataframe entry
                for read in list(position_dt[pos].keys()):
                    # df.at[row,col]= new val
                    df.at[read, pos] += 1

        # remove columns with less than min_reads reads supporting that methyl position
        # df = df.loc[:, df.sum(axis=0) > self.min_reads]

        # make_overlap_links(df, read_set)
        # run_louvain(df)

        make_quick_hm(df, self.sample_name)
        make_pca(df, self.sample_name)

        # select only the region, instead of the flanks with low coverage artifact of pysam fetch
        df_region_only = df.loc[:, start:end]
        make_quick_hm(df_region_only, self.sample_name + '.reads.region')
        make_pca(df_region_only, self.sample_name + '.region')

        samfile.close()

        logfile.close()

        # write_read_mods_to_bed(chrom_modbase_pos)

        return self.mod_base_in_read, df_region_only, read_strands

    def parse_mm_tag(self, read, mm_tag):
        """ mm tag needs to be split and skip bases list isolated: Mm C+m?,1,1,6,5;
            ml tag is stored as an array from the get_tag method, mm is not """
        mm_arr = read.get_tag(mm_tag)
        # TODO: Consider storing the +, in the event it is ever a - instead
        self.modbase = mm_arr[0]
        # replace the mod base with its complement if the alignment is reverse
        if read.is_reverse:
            self.modbase = rev_comp(self.modbase)
        # semicolon delimiter at the end of the list separating the tags is removed
        self.mm = mm_arr[:-1].split(',')[1::]

    def check_modbase_rev(self, seq, mod_base_i, mod_base_i_index):
        # sanity check that that base is the modbase from the Mm Tag
        if seq[mod_base_i[-mod_base_i_index] - 1] != rev_comp(self.modbase).upper():
            print('not cpg', seq[mod_base_i[-mod_base_i_index] - 1])
        if seq[mod_base_i[-mod_base_i_index]].upper() != self.modbase.upper():
            print('base is not modbase', seq[:-mod_base_i_index].upper())
            return -1

    def check_modbase(self, seq, mod_base_i, read_modbase_pos):
        # sanity check that that base is the modbase from the Mm Tag
        if seq[mod_base_i[read_modbase_pos] + 1].upper() != rev_comp(self.modbase).upper():
            print('not cpg', seq[mod_base_i[read_modbase_pos] + 1])
        if seq[mod_base_i[read_modbase_pos]].upper() != self.modbase.upper():
            print('base is not modbase', seq[:read_modbase_pos].upper())
            return -1

    def fill_dictionary(self, read, chrom_modbase_pos, pos, info_list, logfile):
        """ Fill the 2 dictionaries I'm currenlty using, hopefully make into just one soon """
        # check that the chromosome is a key already
        # chrom_modbase_pos: {ref_chrom : ref position : readName :
        # [read.query_alignment_start, read.query_alignment_end, position in mm tag, ml value, position in read, strand] }
        if read.reference_name not in chrom_modbase_pos.keys():
            chrom_modbase_pos[read.reference_name] = {pos: {read.query_name: info_list}}

        # check if position is in the dictionary
        if pos not in chrom_modbase_pos[read.reference_name].keys():
            chrom_modbase_pos[read.reference_name][pos] = {read.query_name: info_list}
            logfile.write(
                '\t' + str(info_list[2]) + ' -strand pos ' + '\t' + str(pos) + ' ' + str(read.query_name) + '\n')

        elif pos in chrom_modbase_pos[read.reference_name].keys():
            chrom_modbase_pos[read.reference_name][pos][read.query_name] = info_list
            logfile.write(
                '\t' + str(info_list[2]) + ' - existing strand pos ' + '\t' + str(pos) + ' ' + str(read.query_name) + '\n')

        # dict #2 don't need this soon
        # check for query name in dict keys
        if read.query_name not in self.mod_base_in_read.keys():
            # create a dictionary entry with the read as key
            # read mod base dict : {readName : ref chrom : ref position : [position in mm tag, ml value, position in read, strand]}
            self.mod_base_in_read[read.query_name] = {read.reference_name: {pos: info_list}}

        if pos not in self.mod_base_in_read[read.query_name][read.reference_name].keys():
            self.mod_base_in_read[read.query_name][read.reference_name][pos] = info_list


def make_quick_hm(data, title):
    ''' sort and heatmap'''

    data.sort_values(by=list(data), axis=0, inplace=True, ascending=False)
    # print('post sort\n',data)
    data.to_csv('methyl_dataframe.' + title + '.csv')
    hm = sn.heatmap(data=data, yticklabels=1, xticklabels=1)

    plt.title(title + ' Heatmap')
    # plt.tick_params( labelleft=False)
    # # this won't be right, use the values ....
    # first_tick = 0
    # last_tick = 280
    # print([data.columns[0], data.columns[-1]],first_tick, last_tick)
    # plt.xticks( [first_tick, last_tick], [data.columns[0], data.columns[-1]]) #, rotation = 'horizontal' )
    plt.subplots_adjust(bottom=0.2)
    # plt.tight_layout()
    plt.savefig(title + ".hm.png")
    plt.clf()


def make_pca(data, title):
    fig = plt.figure(1, figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d", elev=-150, azim=110)
    pca = PCA(n_components=3)
    pca.fit(data)
    X_reduced = pca.transform(data)
    ax.scatter(
        X_reduced[:, 0],
        X_reduced[:, 1],
        X_reduced[:, 2],
        cmap=plt.cm.Set1,
        edgecolor="k",
        s=40,
    )

    ax.set_title(title)
    ax.set_xlabel(str(round(list(pca.explained_variance_ratio_)[0], 4)) + "%")
    ax.xaxis.set_ticklabels([])
    ax.set_ylabel(str(round(list(pca.explained_variance_ratio_)[1], 4)) + "%")
    ax.yaxis.set_ticklabels([])
    ax.set_zlabel(str(round(list(pca.explained_variance_ratio_)[2], 4)) + "%")
    ax.zaxis.set_ticklabels([])

    plt.savefig(title + ".pca3D.png")
    plt.clf()


def run_louvain(df):
    ''' run louvain clustering on my dataframe  '''
    # print('louvain:')
    sdf = df.astype(pd.SparseDtype('int', 0))
    cols = list(df)
    rows = list(df.index)

    # print('sparse',sdf.head())
    biadj = sdf.sparse.to_coo()
    biadj = sparse.csr_matrix(biadj)

    louvain = Louvain()
    louvain.fit(biadj)
    # labels_row = louvain.labels_row_
    # labels_col = louvain.lables_col_

    image = svg_bigraph(biadj, rows, cols, filename='biadj_graph.svg')


def get_chrom(read):
    if len(list(read.keys())) == 1:
        return list(read.keys())[0]
    else:
        print(read.keys())


def find_overlap_index(a, b):
    ''' find the index where b begins to overlap a '''
    for i, ai in enumerate(a):
        if ai >= b[0]:
            return i


def find_overlap_end_index(a, b):
    ''' find the index where b ends overlapping a '''
    end_index = None
    for i, ai in enumerate(a):
        if ai >= b[-1]:
            end_index = i - 1
    if not end_index:
        end_index = len(a) - 1
    return end_index


def get_freq_at_each_pos(df, sample_name):
    """ store the frequency of each methylation character at each position for easy look up """

    positional_frequencies = pd.DataFrame(np.nan, index=[-1, 0, 1], columns=list(df.columns))
    # get a count for each character at each position
    for pos in df.columns:
        vdf = df.loc[:, pos].value_counts()
        for c in vdf.index:
            positional_frequencies.at[c, pos] = vdf[c]

    print(positional_frequencies)

    positional_frequencies.T.plot.bar().get_figure().savefig(sample_name + '.fqsT.bar.png')
    plt.clf()

    return positional_frequencies


def make_overlap_edge_weights(reads_dict, rdf, read_strands, sample_name):
    ''' find reads that overlap methyl positions and score the overlap '''
    # consider dataframe of positions, rows are reads, where alignments are already lined up
    # calculate fq of C, M at each position in overlap between the two reads

    scores = []
    reads = list(reads_dict.keys())
    read_hap_labels = ['_'.join([key, str(items[0]), str(items[1])]) for key, items in read_strands.items()]
    read_scores = pd.DataFrame(np.nan, index=read_hap_labels, columns=read_hap_labels)

    # Conditional Probability Distribution
    cpd = pd.DataFrame({'C': [0.975, 0.025],
                        'M': [0.025, 0.975]},
                       index=[0, 1])
    print('cpd', cpd)

    # dynamic proramming speeds up frequency look up at each position
    positional_frequencies = get_freq_at_each_pos(rdf, sample_name)

    # iterate over all pairs of reads
    for a, b in itertools.combinations(reads_dict, 2):
        if get_chrom(reads_dict[a]) == get_chrom(reads_dict[b]):
            a_b_positional_values = []
            a_read_alignment = sorted(list(reads_dict[a][get_chrom(reads_dict[a])].keys()))
            b_read_alignment = sorted(list(reads_dict[b][get_chrom(reads_dict[b])].keys()))

            # if b overlaps a; calculate probabilities
            if (a_read_alignment[0] < b_read_alignment[0]) and (a_read_alignment[-1] > b_read_alignment[0]):

                overlap_start = find_overlap_index(a_read_alignment, b_read_alignment)
                overlap_end = find_overlap_end_index(a_read_alignment, b_read_alignment)
                # make sure the end of the overlap doesn't go past the end of b
                if overlap_end:
                    if (overlap_end - overlap_start) > len(b_read_alignment):
                        overlap_end = overlap_start + len(b_read_alignment) - 1

                # isolate the dataframe to only where these 2 reads overlap
                read_pair_df = rdf.loc[[a, b], a_read_alignment[overlap_start]:a_read_alignment[overlap_end]]

                # iterate through all positions that overlap
                for i, aPos in enumerate(read_pair_df.columns):
                    if aPos in rdf.columns:
                        # store the methyl character for each read at the current position
                        a_char = rdf.loc[a, aPos]
                        b_char = rdf.loc[b, aPos]

                        # read coverage: total number of reads with 1 or 0 at this position
                        read_cov_at_pos = positional_frequencies.loc[[0, 1], aPos].sum()

                        # random overlap ( denominator )
                        # frequency of the symbol at position i in read a and read b.
                        qxi = positional_frequencies.loc[a_char, aPos] / read_cov_at_pos
                        qyi = positional_frequencies.loc[b_char, aPos] / read_cov_at_pos

                        # cdp probabilities
                        # If there are 0(cannonical C) or 1(methylated C) characters at this position calculate their frequency
                        # for c in cpd.columns:
                        qci = 0
                        qmi = 0
                        # if there is a value there, calculate the frequency of 0(C) and 1(5mC) at this position
                        if not np.isnan(positional_frequencies.loc[0, aPos]):
                            qci = positional_frequencies.loc[0, aPos] / read_cov_at_pos  # rdf.shape[0]
                        if not np.isnan(positional_frequencies.loc[1, aPos]):
                            qmi = positional_frequencies.loc[1, aPos] / read_cov_at_pos  # rdf.shape[0]

                        # multiply the frequency by the cpd value for that character for both reads
                        if a_char != -1 and b_char != -1:
                            pxi_giv_c = qci * cpd.loc[a_char, 'C'] * cpd.loc[b_char, 'C']
                            pxi_giv_m = qmi * cpd.loc[a_char, 'M'] * cpd.loc[b_char, 'M']
                            # then add them together
                            # print('\n', aPos, 'P(xy|M)', pxi_giv_c, '+', pxi_giv_m, '=', pxi_giv_c + pxi_giv_m )
                            # print('\tP(xy|R)', qxi * qyi)
                            pxi_yi_sum = pxi_giv_c + pxi_giv_m
                            px_y_r_prod = qxi * qyi
                            # store the log of the ratio
                            a_b_positional_values.append(np.log(pxi_yi_sum / px_y_r_prod))

                    else:
                        print('not in df?', aPos)

        score = sum(a_b_positional_values)

        a_label = "_".join([a, str(read_strands[a][0]), str(read_strands[a][1])])
        b_label = "_".join([b, str(read_strands[b][0]), str(read_strands[b][1])])
        print(a_label, b_label, 'score:', score)
        read_scores.at[a_label, b_label] = score
        scores.append(score)
        plt.scatter(a, score)

    plt.savefig(sample_name + "score.scatter.png")
    plt.clf()

    print(read_scores)
    make_quick_hm(read_scores, sample_name + '.scores.region')
    make_pca(read_scores.fillna(0), sample_name + '.region.scores')
    ax = read_scores.plot.hist(legend=False, colormap='Blues').get_figure().savefig(sample_name + '.scores.hist.png')
    plt.hist(scores)
    plt.title(sample_name+" Scores Histogram")
    plt.savefig(sample_name + ".plt.scores_hist.png")
    plt.clf()


def graph(edges, read_strands):
    ''' make the graph from edge list '''
    graph = from_edge_list(edges)
    adjacency = graph.adjacency
    print(adjacency.shape)
    # print(graph.names)

    algorithm = Spring()
    position = algorithm.fit_transform(adjacency)
    image = svg_graph(adjacency, position, names=graph.names, filename='spring_adj_graph.svg')

    algo = Louvain()
    labels = algo.fit_predict(adjacency)
    print(labels)

    label_list = []
    for i, label in enumerate(labels):
        label_list.append((label, read_strands[graph.names[i]][0], read_strands[graph.names[i]][1], graph.names[i]))

    label_list.sort(key=lambda x: x[0])
    for lab in label_list:
        print(lab)

    # print(position)
    image = svg_graph(adjacency, position, names=graph.names, labels=labels, filename='louvain_adj_graph.svg')


def main(bamPath, min_prob, min_reads, sample_name, region):

    # methyl probability filter
    ml_min_prob = min_prob
    ml_cutoff = ml_min_prob * 255

    # run methyl analyze per read using url
    mb = MethylBam(bamPath, ml_cutoff, min_reads, sample_name, region)
    mod_base_in_read, rdf, read_strands = mb.get_methyl_per_read()

    # plot_modbase_locations(mod_base_in_read)

    # make_overlap_edges(mod_base_in_read, 10,read_strands)

    make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, sample_name)


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
        help="sample name "
    )

    parser.add_argument(
        "-r",
        nargs=3,
        required=True,
        type=str,
        help="space separated list of region conponents: chromosome start end"
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

    args = parser.parse_args()
    main(bamPath=args.b, min_prob=args.p, min_reads=args.c, sample_name=args.s, region=args.r)


# todo: figure out which scores are going to which reads, and somehow tag the reads with scores that
#        match the groups I am expecting.