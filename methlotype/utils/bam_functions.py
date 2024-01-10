#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
import itertools

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
        if c[0] in [4,5]:
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


def get_chrom(read):
    if len(list(read.keys())) == 1:
        return list(read.keys())[0]
    else:
        print(read.keys())


def find_overlap_index(a, b):
    """ find the index where b begins to overlap a """
    for i, ai in enumerate(a):
        if ai >= b[0]:
            return i


def find_overlap_end_index(a, b):
    """ find the index where b ends overlapping a """
    end_index = None
    if a[-1] == b[-1]:
        end_index = len(a) - 1
    else:
        for i, ai in enumerate(a):
            # print('ai', ai, b[-1])
            if ai >= b[-1]:
                end_index = i - 1
                # print('overlap end', end_index)
    if not end_index:
        end_index = len(a) - 1
        # print('overlap end changed', end_index)
    return end_index


# some graphing functions
def make_quick_hm(data, title):
    ''' sort and heatmap'''

    # this is the SettingWithCopyWarning:
    # A value is trying to be set on a copy of a slice from a DataFrame
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


