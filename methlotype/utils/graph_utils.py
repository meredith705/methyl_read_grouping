#!/usr/bin/env python3

import numpy as np
import pandas as pd

# louvain
from scipy import sparse
from sknetwork.clustering import Louvain, get_modularity
from sknetwork.data import from_edge_list, from_adjacency_list, from_graphml, from_csv
from sknetwork.visualization import svg_graph, svg_bigraph
from sknetwork.embedding import Spring
from sknetwork.linalg import normalize
from sknetwork.utils import get_membership

from IPython.display import SVG

def make_graph(edges, read_strands, sample):
    ''' make the graph from edge list '''

    graph = from_edge_list(edges)#, directed=True, weighted=False)
    adjacency = graph.adjacency
    print('adjacency.shape', adjacency.shape)

    algorithm = Spring()
    position = algorithm.fit_transform(adjacency)
    image = svg_graph(adjacency, position, names=graph.names, filename=sample+'_spring_adj_graph.svg')

    algo = Louvain()
    labels = algo.fit_predict(adjacency)
    print('labels', labels)

    label_list = []
    node_names = []
    for i, label in enumerate(labels):
        # label_list.append((label, read_strands[graph.names[i]][0], read_strands[graph.names[i]][1], graph.names[i]))
        label_list.append( '_'.join([str(label), graph.names[i][-3:]]) )

    # label_list.sort(key=lambda x: x[0])
    # for i,lab in enumerate(label_list):
    #     # print(lab, graph.names[i])
    #     print(lab)

    # print(position)
    # print(label_list)
    image = svg_graph(adjacency, position, names=label_list, labels=labels, filename=sample+'_louvain_adj_graph.svg')

    # bipartite:
    graphb = from_edge_list(edges, bipartite=True)
    biadj = graphb.biadjacency
    names = graphb.names
    names_col = graphb.names_col
    bimage = svg_bigraph(biadj, names_row=names, names_col=names_col, filename=sample+'_biadj_graph.svg')


def run_louvain(df, sample):
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

    image = svg_bigraph(biadj, rows, cols, filename=sample+'biadj_graph.svg')