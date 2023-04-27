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


def identify_metyl_tags(read):
	""" the methyl tags can be either Mm/Ml or MM/ML for R10 or R9"""
	if read.has_tag('Mm') and read.has_tag('Ml'):
		mm_tag='Mm'
		ml_tag='Ml'
	elif read.has_tag('MM') and read.has_tag('ML'):
		mm_tag='MM'
		ml_tag='ML'
	else:
		print("Read", read.query_name, " does not have either Mm/Ml or MM/ML tags.")
		return -1, -1

	return mm_tag, ml_tag


def rev_comp(base):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return complement.get(base,base)


def parse_mm_tag(read, mm_tag):
	""" mm tag needs to be split and skip bases list isolated: Mm C+m?,1,1,6,5;
		ml tag is stored as an array from the get_tag method, mm is not """
	mm_arr = read.get_tag(mm_tag)
	# TODO: Consider storing the +, in the event it is ever a - instead
	modbase = mm_arr[0]
	# replace the mod base with it's complement if the alignment is reverse
	if read.is_reverse:
		modbase = rev_comp(modbase)
	# semicolon delimiter at the end of the list separating the tags is removed
	mm = mm_arr[:-1].split(',')[1::]

	return modbase, mm


def get_hp_tag(read):
	if read.has_tag('HP'):
		return read.get_tag('HP')
	else:
		return 0


def check_modbase_rev(seq,modbase,mod_base_i,mod_base_i_index):
	# sanity check that that base is the modbase from the Mm Tag
	if seq[ mod_base_i[-mod_base_i_index]-1] != rev_comp(modbase).upper():
		print('not cpg',seq[ mod_base_i[-mod_base_i_index]-1])
	if seq[ mod_base_i[-mod_base_i_index] ].upper() != modbase.upper():
		print( 'base is not modbase',seq[:-mod_base_i_index].upper() )
		return -1


def check_modbase(seq,modbase,mod_base_i,read_modbase_pos):
	# sanity check that that base is the modbase from the Mm Tag
	if seq[ mod_base_i[read_modbase_pos]+1 ].upper() != rev_comp(modbase).upper():
		print('not cpg',seq[ mod_base_i[read_modbase_pos]+1 ])
	if seq[ mod_base_i[read_modbase_pos] ].upper() != modbase.upper():
		print( 'base is not modbase',seq[:read_modbase_pos].upper() )
		return -1


def get_cigar_ref_pos(read):
	""" parse the cigar string to obtain methylated cytosine positions
		in refernce coordinates for each read """


	pos=read.reference_start#+1
	refcigarpositions = []#[read.reference_start+1]

	# for anti-sense strand read alignments the cigar tuples need to be read from back to front,
	# they are created with respect to the reference.
	if read.is_reverse:
		cigartups = read.cigartuples #[::-1]
		pos+=1
	else:
		cigartups = read.cigartuples

	for i,c in enumerate(cigartups):
		# Read comsuming CIGAR operations
		# 4 : S : soft clip; set empty position holders because Mm tag includes these positions; can be H_S_query_S_H at ends
		# 5 : H : hard clip; clipped sequences NOT present in SEQ don't count these positions;
		# 			only present as the first and/or last operation.
		# not sure if H clipped sequences are included in the Mm tag; if so add more ref position placeholders.
		if c[0] == 4:
			if i<2:
				refcigarpositions[0:0] = [0]*c[1]
			else:
				refcigarpositions.extend([0]*c[1])
		# 1 : INS; add place holder ref coordinate for each read only base consumed by the insertion
		if c[0] in [1]:
			refcigarpositions.extend([pos]*c[1])

		# 0 : M :  match/mismatch add ref positions to the array for every M base
		# 7 : = :  add ref positions to the array for every = base
		# 8 : X :  mismatch add ref position as if it was =
		if c[0] in [0,7,8]:
			# if read.query_name == 'ea1a4055-aa66-4319-8cec-ea002f7c8ba4' and len(refcigarpositions) < 100:
			# 	print('M:',c,pos,refcigarpositions)
			refcigarpositions.extend( [i for i in np.arange(pos,pos+c[1]) ] )
			pos+=c[1]

		# 2 : D :  DEL;skip the ref position
		# 3 : N : REF_SKIP;
		if c[0] in [2,3]:
			pos+=c[1]

	return refcigarpositions


def fill_dictionary(read, chrom_modbase_pos, mod_base_in_read, pos, info_list, logfile):
	""" Fill the 2 dictionaries I'm currenlty using, hopefully make into just one soon """
	# check that the chromosome is a key already
	# chrom_modbase_pos: {ref_chrom : ref position : readName :
	# [read.query_alignment_start, read.query_alignment_end, position in mm tag, ml value, position in read, strand] }
	if read.reference_name not in chrom_modbase_pos.keys():
		chrom_modbase_pos[read.reference_name] = {	pos: {read.query_name: info_list}}

	# check if position is in the dictionary
	if pos not in chrom_modbase_pos[read.reference_name].keys():
		chrom_modbase_pos[read.reference_name][pos] = {read.query_name: info_list}
		logfile.write(
			'\t' + str(info_list[2]) + ' -strand pos ' + '\t' + str(pos) + ' ' + str(read.query_name) + '\n')

	elif pos in chrom_modbase_pos[read.reference_name].keys():
		chrom_modbase_pos[read.reference_name][pos][read.query_name] = info_list
		logfile.write('\t' + str(info_list[2]) + ' - existing strand pos ' + '\t' + str(pos) + ' ' + str(read.query_name) + '\n')

	# dict #2 don't need this soon
	# check for query name in dict keys
	if read.query_name not in mod_base_in_read.keys():
		# create a dictionary entry with the read as key
		# read mod base dict : {readName : ref chrom : ref position : [position in mm tag, ml value, position in read, strand]}
		mod_base_in_read[read.query_name] = {read.reference_name: {pos: info_list}}

	if pos not in mod_base_in_read[read.query_name][read.reference_name].keys():
		mod_base_in_read[read.query_name][read.reference_name][pos] = info_list


def get_methyl_per_read(filePath, ml_cutoff, min_reads, mod_base_in_read, sample_name):
	""" iterate over all the read mapping to a specified region using fetch().
	Each iteration returns a AlignedSegment object which represents a single read along with its fields and optional tags:
	"""
	# open alignment file
	samfile = pysam.AlignmentFile(filePath, "rb")

	logfile = open('log.txt', 'w')

	read_set = set()
	read_strands = {}
	chrom_modbase_pos = {}

	# Fetch reads from region: chr20:1000000-1020000
	# neuron loc chr20:63664503-63664830 chr20:63664503-63678780  'chr20',63664503,63678780):
	# neuron loc #1: chr2:51026934-51027176
	# SNRPN chr15:24953000-24958133 'chr15', 24953000, 24958133):     'chr15',24955941,24955999):#    'chr15', 24954860,24956812):#

	mmTag, mlTag = '', ''

	# TODO This won't always be per region....
	# TODO maybe filter out short read alignments, of some cutoff
	for read in samfile.fetch('chr15', 24953000, 24958133):#'chr20',63664503,63678780):#'chr2',51026934,51027176):#

		strand = lambda x: '-' if x else '+'
		# store some easy access meta data
		read_set.add(read.query_name)
		# read_strands[read.query_name]= [strand(read.is_reverse), get_hp_tag(read)]

		bfile = open('beds/'+str(read.query_name)+'.mods.bed','w')
		logfile.write('\n'+str(strand(read.is_reverse))+' '+ str(read.mapq)+' '+ str(read.reference_name)+' '+ str(read.reference_start)+' '+str(read.reference_end)+' '+str(read.query_name) )

		# get mm and ml tag information per read
		# the methyl tags can be either Mm/Ml or MM/ML for R10 or R9
		if len(mmTag) == 0:
			mm_tag, ml_tag = identify_metyl_tags(read)
			if ml_tag == -1:
				return -1

		# mm tag needs to be split and skip bases list isolted:
		modbase, mm = parse_mm_tag(read, mm_tag)
		ml = read.get_tag(ml_tag)

		# store the seq in the proper orientation relative to the Ml tag
		seq = read.seq  #seq = identify_seq_orientation(read)

		# store the index position of every modified base in the read
		mod_base_i = [b.start() for b in re.finditer(modbase.upper(), str(seq).upper())]

		# get reference positions that account for INDELs/Clips by parsing the cigar string
		refcigarpositions = get_cigar_ref_pos(read)

		# store ref start and end positions per read for easy access
		read_strands[read.query_name] = [strand(read.is_reverse), get_hp_tag(read), refcigarpositions[read.query_alignment_start], refcigarpositions[read.query_alignment_end-1]]

		# for a REVERSE(-) alignment move through the alignment from back to front
		# counting rev complement C's (aka G) in increments reported by mm tag
		if read.is_reverse:

			# an index for where the next modified base is
			mod_base_i_index = 0
			for i, i_skip in enumerate(mm):
				# increment the read position by the number of modbases indicies to skip
				# the base after the # to skip is the modified base, so add 1
				mod_base_i_index  += int(i_skip) + 1

				# only count modified bases within the read.query_alignment_start and end positions
				if read.query_alignment_start < mod_base_i[-mod_base_i_index] < read.query_alignment_end:
					# sanity check that the base at this position is the modified base in the mm tag
					check_modbase_rev(seq, modbase, mod_base_i, mod_base_i_index)
					# if the ml probability of being modified is > cutoff,
					if ml[i] >= ml_cutoff:

						# the reference position of the mod base
						# store this as pos - 1 for reverse strands, so the mod base is on the C of the + strand
						pos = refcigarpositions[mod_base_i[-mod_base_i_index]] - 1

						# write to bed:
						bfile.write(str(read.reference_name)+'\t'+str(pos-1)+'\t'+str(pos)+'\t'+str(ml[i])+'\t'+'-'+'\n' )

						# this is outside of the list index

						fill_dictionary(read, chrom_modbase_pos, mod_base_in_read, pos, [refcigarpositions[read.query_alignment_start], refcigarpositions[read.query_alignment_end-1], i, ml[i], mod_base_i[-mod_base_i_index], '-'], logfile)

		# for a FORWARD(+) alignment move through the alignment by moving in the increments reported by the mm tag
		else:
			read_modbase_pos = 0
			for i, i_skip in enumerate(mm):
				# it is the modbase after the # to skip
				if i==0:
					read_modbase_pos=int(i_skip)
				else:
					# increment the read position by the number of modbases indicies to skip
					# ( the number of mod bases to skip to land on the base that is observed as modified )
					read_modbase_pos = read_modbase_pos + int(i_skip) + 1

				# only count modified bases after the read.query_alignment_start position
				if read.query_alignment_start < mod_base_i[read_modbase_pos] < read.query_alignment_end:

					# sanity check that that base is the modbase from the Mm Tag
					check_modbase(seq,modbase,mod_base_i,read_modbase_pos)

					# if the ml probability of being modified is > cutoff,
					# store the modified base location as it 0-based position in the read sequence string
					if ml[i] >= ml_cutoff:

						pos = refcigarpositions[mod_base_i[read_modbase_pos]]+1

						# write to bed
						bfile.write(str(read.reference_name)+'\t'+str(pos-1)+'\t'+str(pos)+'\t'+str(ml[i])+'\t'+'+'+'\n' )

						fill_dictionary(read, chrom_modbase_pos, mod_base_in_read, pos,
										[refcigarpositions[read.query_alignment_start], refcigarpositions[read.query_alignment_end-1], i, ml[i], mod_base_i[read_modbase_pos], '+'], logfile)

		bfile.close()

	# for each chromosome make a new methylation df
	for chrom,position_dt in chrom_modbase_pos.items():
		# make an empty dataframe of the correct size; set this to all -1
		df = pd.DataFrame(-1, index=read_set, columns = sorted(list(position_dt.keys())))

		# using query start and end account for reads not overlapping a position to be counted
		# separate from unmethylatd positions in an overlapping read
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
	# df = df.loc[:, df.sum(axis=0) > min_reads]

	# make_overlap_links(df, read_set)
	# run_louvain(df)

	make_quick_hm(df, sample_name)
	# make_pca(df, sample_name)

	samfile.close()

	logfile.close()

	# write_read_mods_to_bed(chrom_modbase_pos)

	return mod_base_in_read, df, read_strands


def make_quick_hm(data,title):
	''' sort and heatmap'''

	data.sort_values(by=list(data), axis=0, inplace=True, ascending=False)
	# print('post sort\n',data)
	data.to_csv('methyl_dataframe.csv')
	hm = sn.heatmap(data = data, yticklabels=1, xticklabels=1)

	plt.title(title+' Heatmap')
	# plt.tight_layout()
	plt.savefig(title+".hm.png")
	plt.clf()
	# try clustermap
	# cm = sn.clustermap(data=data, annot=True, fmt='d')
	# cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
	# plt.title(title+' Clustermap')
	# # plt.tight_layout()
	# plt.savefig(title + ".cm.png")
	# plt.clf()


def make_pca(data,title):
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
	ax.set_xlabel(str(round(list(pca.explained_variance_ratio_)[0],4))+"%")
	ax.xaxis.set_ticklabels([])
	ax.set_ylabel(str(round(list(pca.explained_variance_ratio_)[1],4))+"%")
	ax.yaxis.set_ticklabels([])
	ax.set_zlabel(str(round(list(pca.explained_variance_ratio_)[2],4))+"%")
	ax.zaxis.set_ticklabels([])

	plt.savefig(title+".pca3D.png")
	plt.clf()


def run_louvain(df):
	''' run louvain clustering on my dataframe  '''
	# print('louvain:')
	sdf = df.astype(pd.SparseDtype('int',0))
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
	# SVG(image)


def get_chrom(read):
	if len(list(read.keys())) == 1:
		return list(read.keys())[0]
	else:
		print(read.keys())


def exact_overlap_k(pos1,pos2,k):
	''' check if k suffix/prefix methylated positions overlap'''
	return pos1[-k:]==pos2[:k]


def hamming_overlap_k(pos1,pos2,k):
	''' check if k suffix/prefix methylated positions overlap'''
	return hamming(pos1[-k:],pos2[:k])


def find_overlap_index(a,b):
	''' find the index where b begins to overlap a '''
	for i,ai in enumerate(a):
		if ai>=b[0]:
			return i


def find_overlap_end_index(a,b):
	''' find the index where b ends overlapping a '''
	end_index = None
	for i,ai in enumerate(a):
		if ai>=b[-1]:
			end_index = i-1
	if not end_index:
		end_index = len(a)-1
	return end_index


def make_overlap_edges(reads_df,k, read_strands):
	''' find reads that overlap methyl positions '''
	# TODO:
	# consider dataframe of positions, rows are reads, where alignments are already lined up

	## Edge Parameters
	match_count_min = 100
	match_ratio_min = 0.70

	match_percentages = open('matchInfo.tsv','w')
	links = []
	for a,b in itertools.combinations(reads_df,2):
		# print('a',reads_df[a],'b',reads_df[b])

		if get_chrom(reads_df[a])==get_chrom(reads_df[b]):
			a_pos = sorted(list(reads_df[a][get_chrom(reads_df[a])].keys()))
			b_pos = sorted(list(reads_df[b][get_chrom(reads_df[b])].keys()))
			# print('a',a_pos[:10],reads_df[a][get_chrom(reads_df[a])][a_pos[0]])
			# print('b',b_pos[:10],reads_df[b][get_chrom(reads_df[b])][b_pos[0]])
			# check if reads are aligned in the same orientation
			# if reads_df[a][get_chrom(reads_df[a])][a_pos[0]][-1] == reads_df[b][get_chrom(reads_df[b])][b_pos[0]][-1]:
				# print(a,b)
			# if b overlaps a check for suffix of a matching prefix of b
			if (a_pos[0] < b_pos[0]) and (a_pos[-1] > b_pos[0]):

				overlap_start = find_overlap_index(a_pos,b_pos)
				overlap_end   = find_overlap_end_index(a_pos,b_pos)

				if overlap_end:
					if (overlap_end - overlap_start) > len(b_pos):
						overlap_end = overlap_start+len(b_pos)-1

				# print('b over a: overlap_start',overlap_start,'-',overlap_end,len(a_pos),'length of overlap to end of a match',len(a_pos[overlap_start:overlap_end]),'b',len(b_pos))
				# if exact_overlap_k(a_pos[overlap_start::],b_pos,len(a_pos[overlap_start::])):
				# print('exact:',exact_overlap_k(a_pos[overlap_start:overlap_end],b_pos[:len(a_pos[overlap_start:overlap_end])],len(a_pos[overlap_start::])-5) )

				# print('ham:',hamming_overlap_k(a_pos[overlap_start:overlap_end],b_pos[:len(a_pos[overlap_start:overlap_end])],len(a_pos[overlap_start:overlap_end])-5) )

				# count each position that matches between these reads
				match_count = 0
				for i, aa in enumerate(a_pos[overlap_start:overlap_end]):
					if aa in b_pos[:len(a_pos[overlap_start:overlap_end])] :# or aa+1 in b_pos or aa-1 in b_pos:
						# print(aa, b_pos.index(aa)) #,b_pos.index(aa+1),b_pos.index(aa-1))
						match_count+=1
					elif aa+1 in b_pos[:len(a_pos[overlap_start:overlap_end])]:
						# print(aa,'+1',b_pos.index(aa+1))
						match_count+=1
					elif aa-1 in b_pos[:len(a_pos[overlap_start:overlap_end])]:
						# print(aa,'-1',b_pos.index(aa-1))
						match_count+=1
				if len(a_pos[overlap_start:overlap_end])>100:
					# print(a,b,'match_count:',match_count, match_count/len(a_pos[overlap_start:overlap_end]))
					match_percentages.write(('\t').join([a,b,str(match_count), str(round(match_count/len(a_pos[overlap_start:overlap_end]),4)),'\n']) )

					if match_count > match_count_min and match_count/len(a_pos[overlap_start:overlap_end]):
						links.append((a,b))


			elif (b_pos[0] < a_pos[0]) and (b_pos[-1] > a_pos[0]):
				overlap_start = find_overlap_index(b_pos,a_pos)
				overlap_end   = find_overlap_end_index(b_pos,a_pos)
				# print('a over b: overlap_start',overlap_start,len(b_pos[overlap_start:overlap_end]),a,b)
				# if exact_overlap_k(a_pos[overlap_start::],b_pos,len(a_pos[overlap_start:overlap_end])):
				# print(exact_overlap_k(b_pos[overlap_start::],a_pos,len(b_pos[overlap_start:overlap_end])))
				# print(hamming_overlap_k(b_pos[overlap_start:overlap_end],a_pos,len(b_pos[overlap_start:overlap_end])))
				# count each position that matches between these reads
				match_count = 0
				for i, aa in enumerate(a_pos[overlap_start:overlap_end]):
					if aa in b_pos[:len(a_pos[overlap_start:overlap_end])] :# or aa+1 in b_pos or aa-1 in b_pos:
						# print(aa, b_pos.index(aa)) #,b_pos.index(aa+1),b_pos.index(aa-1))
						match_count+=1
					elif aa+1 in b_pos[:len(a_pos[overlap_start:overlap_end])]:
						# print(aa,'+1',b_pos.index(aa+1))
						match_count+=1
					elif aa-1 in b_pos[:len(a_pos[overlap_start:overlap_end])]:
						# print(aa,'-1',b_pos.index(aa-1))
						match_count+=1
				if len(a_pos[overlap_start:overlap_end])>100:
					# print(a,b,'match_count:',match_count, match_count/len(a_pos[overlap_start:overlap_end]))
					match_percentages.write(('\t').join([a,b,str(match_count), str(round(match_count/len(a_pos[overlap_start:overlap_end]),4)),'\n']) )
					if match_count > match_count_min and match_count/len(a_pos[overlap_start:overlap_end]):
						links.append((a,b))
				# print(b_pos[overlap_start::],'\n',a_pos)
						# break
	match_percentages.close()

	# print(links)
	graph(links, read_strands)


def get_freq_at_each_pos(df):
	""" store the frequency of each methylation character at each position for easy look up """

	positional_frequencies = pd.DataFrame(np.nan, index=[-1,0,1], columns=list(df.columns) )

	for pos in df.columns:
		vdf = df.loc[:, pos].value_counts()
		for c in vdf.index:
			positional_frequencies.at[c, pos] = vdf[c]

	print(positional_frequencies)
	return positional_frequencies



def make_overlap_edge_weights(reads_dict, rdf, read_strands, sample_name):
	''' find reads that overlap methyl positions '''
	# TODO: maybe limit the frequency calc to the fetch columns instead of the whole read coverage df
	# consider dataframe of positions, rows are reads, where alignments are already lined up\
	# calculate fq of C, M at each position in overlap between the two reads
		# should be done using pandas in the rdf column

	# match_percentages = open('matchInfo.tsv','w')
	links = []
	read_scores = pd.DataFrame(np.nan, index=list(reads_dict.keys()), columns = list(reads_dict.keys() ) )
	read_score_list =[]

	# Conditional Probability Distribution
	cpd = pd.DataFrame({'C': [0.95, 0.05],
						'M': [0.05, 0.95]},
						index=[0, 1])
	print('cpd',cpd)

	positional_frequencies = get_freq_at_each_pos(rdf)

	# iterate over all pairs of reads
	for a, b in itertools.combinations(reads_dict,2):

		if get_chrom(reads_dict[a])==get_chrom(reads_dict[b]):
			a_b_positional_values = []
			a_pos = sorted(list(reads_dict[a][get_chrom(reads_dict[a])].keys()))
			b_pos = sorted(list(reads_dict[b][get_chrom(reads_dict[b])].keys()))

			# if b overlaps a calculate probabilities
			if (a_pos[0] < b_pos[0]) and (a_pos[-1] > b_pos[0]):

				overlap_start = find_overlap_index(a_pos,b_pos)
				overlap_end   = find_overlap_end_index(a_pos,b_pos)

				if overlap_end:
					if (overlap_end - overlap_start) > len(b_pos):
						overlap_end = overlap_start+len(b_pos)-1

				# calc fq of M, !M at each position in read overlap
				# print('\na & b row', rdf.loc[[a, b], (a_pos[overlap_start]):(a_pos[overlap_end])])
				read_pair_df = rdf.loc[[a, b], a_pos[overlap_start]:a_pos[overlap_end]]

				frequency_product_list = []
				for i, aPos in enumerate(read_pair_df.columns):

					if aPos in rdf.columns:
						# print('a', aPos, [rdf.loc[a, aPos]], 'count:', rdf.loc[:, aPos].value_counts()[rdf.loc[a, aPos]])
						# print('b', [rdf.loc[b, aPos]], 'count', rdf.loc[:, aPos].value_counts()[rdf.loc[b, aPos]])

						# count only 0,1 as denominator of frequency calculation
						# dynamic proramming can help speed up this look up at each position for each read
						# df_values = rdf.loc[:, aPos].value_counts()
						# values = list(df_values.index)

						a_char = rdf.loc[a, aPos]
						b_char = rdf.loc[b, aPos]

						# print(df_values)
						# print('avalue', df_values[[rdf.loc[a, aPos]]] )
						# print('bvalue', df_values[[rdf.loc[b, aPos]]] )

						# count number of reads with 1 or 0 at this position
						read_cov_at_pos = positional_frequencies.loc[:, aPos].sum()			#rdf.loc[:, aPos].value_counts()[values].sum()


						# random overlap ( denominator )
						# frequency of the symbol at position i in read a and read b.           #Wrt the # of reads that cover that position
						# qxi = df_values[a_char] / rdf.shape[0]  # read_cov_at_pos
						# qyi = df_values[b_char] / rdf.shape[0]  # read_cov_at_pos

						qxi = positional_frequencies.loc[a_char, aPos] / rdf.shape[0]  # read_cov_at_pos
						qyi = positional_frequencies.loc[b_char, aPos] / rdf.shape[0]  # read_cov_at_pos

						# cdp probabilities
						# for c in cpd.columns:
						qci = 0
						qmi = 0
						if not np.isnan(positional_frequencies.loc[0, aPos]):
							qci = positional_frequencies.loc[0, aPos] / rdf.shape[0] #read_cov_at_pos
						if not np.isnan(positional_frequencies.loc[1, aPos]):
							qmi = positional_frequencies.loc[1, aPos] / rdf.shape[0] #read_cov_at_pos

						# print(rdf.loc[a, aPos], cpd.loc[rdf.loc[a, aPos], 'C'], rdf.loc[b, aPos], cpd.loc[rdf.loc[b, aPos], 'C'] ) #cpd[rdf.loc[a, aPos], 'C'], cpd[rdf.loc[b, aPos], 'C'])
						# print(rdf.loc[a, aPos], cpd[rdf.loc[a, aPos], 'M'], cpd[rdf.loc[b, aPos], 'M'])

						if a_char != -1 and b_char != -1:
							pxi_giv_c = qci * cpd.loc[a_char, 'C'] * cpd.loc[b_char, 'C']
							pxi_giv_m = qmi * cpd.loc[a_char, 'M'] * cpd.loc[b_char, 'M']
							# then add them together
							# print('\n', aPos, 'P(xy|M)', pxi_giv_c, '+', pxi_giv_m, '=', pxi_giv_c + pxi_giv_m )
							# print('\tP(xy|R)', qxi * qyi)

							pxi_yi_sum = pxi_giv_c + pxi_giv_m
							px_y_r_prod = qxi * qyi
							# print('\t',pxi_yi_sum/ px_y_r_prod,  np.log( pxi_yi_sum/ px_y_r_prod))
							a_b_positional_values.append( np.log( pxi_yi_sum/ px_y_r_prod) )
						else:
							# some problem here where -1 is in the overlap... lok at finding overlap start and overlap end
							print(aPos, a_char, b_char)
						# frequency_product_list.append( qxi * qyi )

					else:
						print('not in df?', aPos)

					# if i >4:
					# 	break
				# print(np.prod(frequency_product_list))
				# print(reduce((lambda x,y: x*y), frequency_product_list))

				# break
		print(a, b, 'score:', sum(a_b_positional_values))
		read_scores.at[a, b] = sum(a_b_positional_values)

	print(read_scores)
	make_quick_hm(read_scores, sample_name+'.scores')
	ax = read_scores.plot.hist(legend=False).get_figure().savefig(sample_name+'.scores.hist.png')

	# match_percentages.close()

	# print(links)
	# graph(links, read_strands)


def graph(edges, read_strands):
	''' make the graph from edge list '''
	graph = from_edge_list(edges)
	adjacency = graph.adjacency
	print(adjacency.shape)
	# print(graph.names)

	algorithm = Spring()
	position = algorithm.fit_transform(adjacency)
	image = svg_graph(adjacency,position,names=graph.names, filename='spring_adj_graph.svg')

	algo = Louvain()
	labels = algo.fit_predict(adjacency)
	print(labels)

	label_list = []
	for i,label in enumerate(labels):
		label_list.append((label,read_strands[graph.names[i]][0],read_strands[graph.names[i]][1], graph.names[i]))

	label_list.sort(key= lambda x: x[0])
	for lab in label_list:
		print( lab)

	# print(position)
	image = svg_graph(adjacency,position,names=graph.names,labels=labels, filename='louvain_adj_graph.svg')


def plot_modbase_locations(mod_base_in_read):
	''' plot the reads '''

	for read,chrom in mod_base_in_read.items():
		print(read)
		c=0
		for chrom,pos in chrom.items():
			print(chrom, list(pos.keys())[:5])

	return


def overlap_matrix(mod_base_in_read):
	''' make a matrix of every cpg site and fill rows with positions that each read is methylated
		reads x cpg sites
	 '''

	ref_cpgs = open('http://public.gi.ucsc.edu/~memeredith/methyl/beds/cpg_hg38.bed')


def main(bamPath, min_prob, min_reads, sample_name):

	ml_min_prob = min_prob
	ml_cutoff = ml_min_prob * 255

	# dictionary to store reads and modified bases in those reads
	mod_base_in_read = {}

	# run methyl analyze per read using url
	mod_base_in_read, rdf, read_strands = get_methyl_per_read(bamPath,ml_cutoff,min_reads, mod_base_in_read,sample_name)
	# mod_base_in_read = get_methyl_per_read("http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam",ml_cutoff, mod_base_in_read)
	# mod_base_in_read = get_methyl_per_read("HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.chr20_run2.haplotagged.bam",ml_cutoff, mod_base_in_read)

	# plot_modbase_locations(mod_base_in_read)

	# make_overlap_edges(mod_base_in_read, 10,read_strands)

	make_overlap_edge_weights(mod_base_in_read, rdf, read_strands, sample_name)






if __name__ == "__main__":
	'''
	Python script to cluster reads in a bam by their methyl tag provile

	Usage: python3.9 analyze_methylBam.py -b http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam
	-s sample_name 
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
		"-p",
		required=False,
		type=float,
		default=0.5,
		help="Minimum probability cutoff for the likelyhood that a position is methylated (ml tag) "
	)

	parser.add_argument(
		"-r",
		required=False,
		type=int,
		default=0,
		help="Minimum number of reads to support a methyl position to be included in analysis "
	)

	args = parser.parse_args()
	main(bamPath=args.b, min_prob=args.p, min_reads=args.r, sample_name=args.s)

	# parser.add_argument(
	# 	"-r",
	# 	required=True,
	# 	type=str,
	# 	help="Input lifted yak bed, sorted by the 4th column of shasta_component_id_pos_endPos"
	# )

	# parser.add_argument(
	# 	"-o",
	# 	required=True,
	# 	type=str,
	# 	help="Output directory"
	# )



	# main(yak_path=args.y, lifted_path=args.r, output_directory=args.o)







#################### try pileup? #########################
''' iterating over each base of a specified region using the pileup() method. Each iteration returns a PileupColumn which represents all the reads in the SAM file that map to a single base in the reference sequence. The list of reads are represented as PileupRead objects in the PileupColumn.pileups property:
	Not easy to get methyl tags from the pile up 
'''
# column_count = 0
# for pileupcolumn in samfile.pileup('chr20', 1000000, 1000020):
# 	column_count+=1
# 	print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n), dir(pileupcolumn))
# 	for pileupread in pileupcolumn.pileups:
# 		if not pileupread.is_del and not pileupread.is_refskip:
# 				# query position is None if is_del or is_refskip is set.
# 				print('\tbase in read %s = %s' % (pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))
# 		print(dir(pileupread))
# 	if column_count > 0: 
# 		break

# samfile.close()

'''
[__... 
'aend', 'alen', 'aligned_pairs', 
'bin', 'blocks', 
'cigar', 'cigarstring', 'cigartuples', 'compare', 
'flag', 'from_dict', 'fromstring', 

'get_aligned_pairs', 'get_blocks', 'get_cigar_stats', 'get_forward_qualities', 'get_forward_sequence', 'get_overlap', 
	'get_reference_positions', 'get_reference_sequence', 'get_tag', 'get_tags', 

'has_tag', 'header', 
'infer_query_length', 'infer_read_length', 'inferred_length', 'is_duplicate', 'is_forward', 'is_mapped', 'is_paired', 
	'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary', 'is_unmapped', 'isize', 

'mapping_quality', 'mapq', 'mate_is_forward', 'mate_is_mapped', 'mate_is_reverse', 'mate_is_unmapped', 'modified_bases', 'modified_bases_forward', 'mpos', 'mrnm', 

'next_reference_id', 'next_reference_name', 'next_reference_start', 
'opt', 'overlap', 
'pnext', 'pos', 'positions', 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end', 'query_alignment_length', 'query_alignment_qualities', 'query_alignment_sequence', 'query_alignment_start', 'query_length', 'query_name', 'query_qualities', 'query_sequence', 
'reference_end', 'reference_id', 'reference_length', 'reference_name', 'reference_start', 'rlen', 'rname', 'rnext', 
'seq', 'setTag', 'set_tag', 'set_tags', 
'tags', 'template_length', 'tid', 'tlen', 'to_dict', 'to_string', 'tostring']
'''


