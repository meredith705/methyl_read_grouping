#!/usr/bin/env python3
import argparse
import os
import pysam
from Bio.Seq import Seq
import re
import numpy as np
import pandas as pd
import itertools
from scipy.spatial.distance import hamming

# heatmap 
import seaborn as sn
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# louvain
from scipy import sparse
from sknetwork.clustering import Louvain, get_modularity
from sknetwork.linalg import normalize
from sknetwork.utils import get_membership
from sknetwork.visualization import svg_graph, svg_bigraph
from IPython.display import SVG


def identify_metyl_tags(read):
	''' the methyl tags can be either Mm/Ml or MM/ML for R10 or R9'''
	if read.has_tag('Mm') and read.has_tag('Ml'):
		mmTag='Mm'
		mlTag='Ml'
	elif read.has_tag('MM') and read.has_tag('ML'):
		mmTag='MM'
		mlTag='ML'
	else:
		print("Read", read.query_name, " does not have either Mm/Ml or MM/ML tags.")
		return -1, -1

	return mmTag, mlTag

def rev_comp(base):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return complement.get(base,base)

def parse_mm_tag(read):
	'''mm tag needs to be split and skip bases list isolted: Mm C+m?,1,1,6,5;
		ml tag is stored as an array from the get_tag method, mm is not'''
	mm_arr = read.get_tag('Mm')
	# TODO: Consider storing the +, in the event it is ever a - instead
	modbase = mm_arr[0]
	# replace the mod base with it's complement if the alignment is reverse
	if read.is_reverse:
		modbase = rev_comp(modbase)
	# semi-colon delimiter at the end of the list separating the tags is removed 
	mm = mm_arr[:-1].split(',')[1::]

	return modbase, mm

def check_modbase_rev(seq,modbase,mod_base_i,mod_base_i_index):
	# sanity check that that base is the modbase from the Mm Tag
	if seq[ mod_base_i[-mod_base_i_index]-1] != rev_comp(modbase).upper():
		print('not cpg',seq[ mod_base_i[-mod_base_i_index]-1])
	if seq[ mod_base_i[-mod_base_i_index] ].upper() != modbase.upper():
		print('base is not modbase',seq[:-mod_base_i_index].upper(), 'at',i,'. i_skip:', i_skip)
		return -1
def check_modbase(seq,modbase,mod_base_i,read_modbase_pos):
	# sanity check that that base is the modbase from the Mm Tag
	if seq[ mod_base_i[read_modbase_pos]+1 ].upper() != rev_comp(modbase).upper():
		print('not cpg',seq[ mod_base_i[read_modbase_pos]+1 ])
	if seq[ mod_base_i[read_modbase_pos] ].upper() != modbase.upper():
		print('base is not modbase',seq[:read_modbase_pos].upper(), 'at',i,'. i_skip:', i_skip)
		return -1

def get_cigar_ref_pos(read):
	''' parse the cigar string to obtain methylated cytosine positions
		in refernce coordinates for each read '''


	pos=read.reference_start#+1
	readcigarpositions = []#[read.reference_start+1]

	# for anti-sense strand read alignments the cigar tuples need to be read from back to front, 
	# they are created with respect to the reference. 
	if read.is_reverse:
		cigartups = read.cigartuples#[::-1]
		pos+=1
	else:
		cigartups = read.cigartuples

	for i,c in enumerate(cigartups):
		# Read comsuming CIGAR operations
		# 4 : S : soft clip; set empty position holders because Mm tag includes these positions; can be H_S_query_S_H at ends 
		# 5 : H : hard clip; clipped sequences NOT present in SEQ don't count these positions; only present as the first and/or last operation.
				# not sure if H clipped sequences are included in the Mm tag; if they are I'll need to add more ref position placeholders.
		if c[0] == 4:
			if i<2: 
				readcigarpositions[0:0] = [0]*c[1]
			else:
				readcigarpositions.extend([0]*c[1])
		# 1 : INS; add place holder ref coordinate for each read only base consumed by the insertion
		if c[0] in [1]:
			readcigarpositions.extend([pos]*c[1])

		# 0 : M :  match/mismatch add ref positions to the array for every M base
		# 7 : = :  add ref positions to the array for every = base
		# 8 : X :  mismatch add ref position as if it was =
		if c[0] in [0,7,8]:
			# if read.query_name == 'ea1a4055-aa66-4319-8cec-ea002f7c8ba4' and len(readcigarpositions) < 100:
			# 	print('M:',c,pos,readcigarpositions)
			readcigarpositions.extend( [i for i in np.arange(pos,pos+c[1]) ] )
			pos+=c[1] 

		# 2 : D :  DEL;skip the ref position
		# 3 : N : REF_SKIP; 
		if c[0] in [2,3]:
			pos+=c[1]

	return readcigarpositions

def write_read_mods_to_bed(chrom_modbase_pos):
	''' write all read mods to a bed file per read'''
	# chrom_modbase_pos: {ref_chrom : ref position : readName : [positon in mm tag, ml value, position in read, strand] }

	# for chrom,position_dt in chrom_modbase_pos.items():
	# 	# make an empty dataframe of the correct size. Maybe the reads can be stored in a set instead of dict
	# 	df = pd.DataFrame(0, index=read_set, columns = sorted(list(position_dt.keys())) )
	# 	# print(df)
	# 	for pos in sorted(list(position_dt.keys())):
	# 		print(pos,len(position_dt[pos].keys()),position_dt[pos].keys() ) 
	# 		for read in list(position_dt[pos].keys()):
	# 			# df.at[row,col]= new val
	# 			df.at[read,pos] += 1
	# 	print('filled\n',df)


	for read,chroms in mod_base_in_read.items():

		bfile = open('beds/'+str(read)+'.mods.bed','w')

		for chrom , refPos in chroms.items():
			for pos, modinfo in refPos.items():
				# adjust the position based on the strand
				if modinfo[3] == '+':
					bfile.write(str(chrom)+'\t'+str(pos)+'\t'+str(pos+1)+'\t'+str(modinfo[1])+'\t'+str(modinfo[3])+'\n' )
				elif modinfo[3] == '-':
					bfile.write(str(chrom)+'\t'+str(pos-1)+'\t'+str(pos)+'\t'+str(modinfo[1])+'\t'+str(modinfo[3])+'\n' )
		bfile.close()

def get_methyl_per_read(filePath, ml_cutoff, min_reads,mod_base_in_read,sample_name):
	''' iterate over all of the read mapping to a specified region using fetch(). 
	Each iteration returns a AlignedSegment object which represents a single read along with its fields and optional tags:
	'''
	# open alignment file
	samfile = pysam.AlignmentFile(filePath, "rb")

	logfile = open('log.txt','w')

	read_set = set()
	chrom_modbase_pos = {}

	# Fetch reads from region: chr20:1000000-1020000
	# neuron loc chr20:63664503-63664830 chr20:63664503-63678780  'chr20',63664503,63678780):
	# neuron loc #1: chr2:51026934-51027176
	# SNRPN chr15:24953000-24958133 'chr15', 24953000, 24958133):
	read_count = 0
	mmTag, mlTag = '', ''

	# TODO This won't always be per region....
	for read in samfile.fetch('chr20',63664503,63678780):
		strand = lambda x : '-' if x else '+'
		# print('\n',strand(read.is_reverse), read.mapq, read.reference_name, read.reference_start,read.reference_end, read.query_name)
		read_set.add(read.query_name)
		bfile = open('beds/'+str(read.query_name)+'.mods.bed','w')
		logfile.write('\n'+str(strand(read.is_reverse))+' '+ str(read.mapq)+' '+ str(read.reference_name)+' '+ str(read.reference_start)+' '+str(read.reference_end)+' '+str(read.query_name) )
		read_count+=1
		# get mm and ml tag information per read
		# the methyl tags can be either Mm/Ml or MM/ML for R10 or R9
		if len(mmTag)==0:
			mmTag, mlTag = identify_metyl_tags(read)
			if mlTag==-1:
				return -1

		# mm tag needs to be split and skip bases list isolted: 
		modbase, mm = parse_mm_tag(read)
		ml = read.get_tag('Ml')

		# store the seq in the proper orientation relative to the Ml tag
		seq = read.seq  #seq = identify_seq_orientation(read)

		# store the index position of every modified base in the read
		mod_base_i = [b.start() for b in re.finditer(modbase.upper(),str(seq).upper())]
	
		# get reference positions that account for INDELs/Clips by parsing the cigar string
		readcigarpositions = get_cigar_ref_pos(read)

		# for a REVERSE(-) alignment move through the alignment from back to front counting rev complement C's (aka G) in increments reported by mm tag
		if read.is_reverse:
			mod_base_i_index = 0 
			for i, i_skip in enumerate(mm):
				# increment the read position by the number of modbases indicies to skip
				# the base after the # to skip is the modified base, so add 1
				mod_base_i_index  += int(i_skip) + 1

				# only count modified bases after the read.query_alignment_start position
				if mod_base_i[-mod_base_i_index] > read.query_alignment_start and mod_base_i[-mod_base_i_index] < read.query_alignment_end:
					check_modbase_rev(seq,modbase,mod_base_i,mod_base_i_index)
					# if the ml probability of being modified is > cutoff,
					if ml[i] >= ml_cutoff:
						# write to bed: 
						bfile.write(str(read.reference_name)+'\t'+str(readcigarpositions[mod_base_i[-mod_base_i_index]]-1)+'\t'+str(readcigarpositions[mod_base_i[-mod_base_i_index]])+'\t'+str(ml[i])+'\t'+'-'+'\n' )

						pos = readcigarpositions[mod_base_i[-mod_base_i_index]]
						# check that the chromomsome is a key already
						# chrom_modbase_pos: {ref_chrom : ref position : readName : [positon in mm tag, ml value, position in read, strand] }
						if read.reference_name not in chrom_modbase_pos.keys():
							chrom_modbase_pos[read.reference_name]={pos : {read.query_name : [ i,ml[i],mod_base_i[-mod_base_i_index],'-' ]} }

						# check if position is in the ditionary 
						if pos not in chrom_modbase_pos[read.reference_name].keys():
							chrom_modbase_pos[read.reference_name][pos] = {read.query_name : [ i,ml[i],mod_base_i[-mod_base_i_index],'-' ]}
							logfile.write('\t'+str(i)+' '+str(i_skip)+' -strand pos '+'\t'+str(pos)+' '+str(read.query_name)+'\n')

						elif readcigarpositions[mod_base_i[-mod_base_i_index]] in chrom_modbase_pos[read.reference_name].keys():

							chrom_modbase_pos[read.reference_name][pos][read.query_name] = [ i,ml[i],mod_base_i[-mod_base_i_index],'-' ]
							logfile.write('\t'+str(i)+' '+str(i_skip)+' - existing strand pos '+'\t'+str(pos)+' '+str(read.query_name)+'\n')

				# if i>5:
				# 	break

		# for a FORWARD(+) alignment move through the alignment by moving in the increments reported by the mm tag
		else:
			read_modbase_pos = 0
			for i, i_skip in enumerate(mm):
				# increment the read position by the number of modbases indicies to skip 
				# ( the number of mod bases to skip to land on the base that is observed as modified )
				read_modbase_pos = read_modbase_pos + int(i_skip)
				
				# it is the modbase after the # to skip
				if i>0:
					read_modbase_pos+=1

				# only count modified bases after the read.query_alignment_start position
				if mod_base_i[read_modbase_pos] > read.query_alignment_start and mod_base_i[read_modbase_pos] < read.query_alignment_end:

					# sanity check that that base is the modbase from the Mm Tag
					check_modbase(seq,modbase,mod_base_i,read_modbase_pos)

					# if the ml probability of being modified is > cutoff, 
					# store the modified base location as it 0-based position in the read sequence string
					if ml[i] >= ml_cutoff:
						# write to bed 
						bfile.write(str(read.reference_name)+'\t'+str(readcigarpositions[mod_base_i[read_modbase_pos]])+'\t'+str(readcigarpositions[mod_base_i[read_modbase_pos]]+1)+'\t'+str(ml[i])+'\t'+'+'+'\n' )

						pos = readcigarpositions[mod_base_i[read_modbase_pos]]+1
						# check that the chromomsome is a key already
						# chrom_modbase_pos: {ref_chrom : ref position : readName : [positon in mm tag, ml value, position in read, strand] }
						if read.reference_name not in chrom_modbase_pos.keys():
							chrom_modbase_pos[read.reference_name]={pos : {read.query_name : [ i,ml[i],mod_base_i[read_modbase_pos],'+' ]} }
						# check if position is in the ditionary 
						if pos not in chrom_modbase_pos[read.reference_name].keys():
							# print('\t',i,i_skip,'+ strand pos:',readcigarpositions[mod_base_i[read_modbase_pos]])
							logfile.write('\t'+str(i)+' '+str(i_skip)+' + strand pos '+str(pos)+' '+str(read.query_name)+'\n' )
							chrom_modbase_pos[read.reference_name][pos] = {read.query_name : [ i,ml[i],mod_base_i[read_modbase_pos],'+' ]}

						elif pos in chrom_modbase_pos[read.reference_name].keys():
							chrom_modbase_pos[read.reference_name][pos][read.query_name] = [ i,ml[i],mod_base_i[read_modbase_pos],'+' ]
							# print('\t',i,i_skip,'+ strand existing pos:',readcigarpositions[mod_base_i[read_modbase_pos]])
							logfile.write('\t'+str(i)+' '+str(i_skip)+' + strand pos '+str(pos)+' '+str(read.query_name)+'\n' )				

		bfile.close()

	for chrom,position_dt in chrom_modbase_pos.items():
		# make an empty dataframe of the correct size. Maybe the reads can be stored in a set instead of dict
		df = pd.DataFrame(0, index=read_set, columns = sorted(list(position_dt.keys())) )
		# print(df)
		for pos in sorted(list(position_dt.keys())):
			# print(pos,len(position_dt[pos].keys()),position_dt[pos].keys() ) 
			for read in list(position_dt[pos].keys()):
				# df.at[row,col]= new val
				if position_dt[pos][read][3] == '-':
					if pos-1 in position_dt.keys():
						if df.at[read,pos-1] > 0:
							# print('\n',read, pos, '\n\t',position_dt[pos],'\n\t',position_dt[pos-1])
							logfile.write('minus strand with more than 1 methyl call at position: '+str(pos)+'\n'+str(read) +'\t'+str(pos) +'\n\t'+str(position_dt[pos])+'\n\t'+str(position_dt[pos-1]))
							df.at[read,pos]+=1
						else:
							df.at[read,pos-1]+=1
					else:
						df.at[read,pos] += 1
				else:
					df.at[read,pos] += 1

	# remove columns with less than 3 reads supporting that methyl position
	df = df.loc[:, df.sum(axis=0) > min_reads]


	# run_louvain(df)
	make_quick_hm(df,sample_name)
	make_pca(df,sample_name)

	print(df.T.corr())
	hm = sn.heatmap(data = df.T.corr())
	plt.title(sample_name)
	plt.savefig(sample_name+".corr.hm.png")
	plt.clf()

	samfile.close()

	logfile.close()

	# write_read_mods_to_bed(chrom_modbase_pos)

	return chrom_modbase_pos

def make_quick_hm(data,title):
	''' sort and heatmap'''

	data.sort_values(by=list(data),axis=0, inplace=True,ascending=False)
	print('post sort\n',data)
	data.to_csv('methyl_dataframe.csv')
	hm = sn.heatmap(data = data)
	plt.title(title)
	plt.savefig(title+".hm.png")
	plt.clf()


def make_pca(data,title):
	fig = plt.figure(1, figsize=(8, 6))
	ax = fig.add_subplot(111, projection="3d", elev=-150, azim=110)
	X_reduced = PCA(n_components=3).fit_transform(data)
	ax.scatter(
	X_reduced[:, 0],
	X_reduced[:, 1],
	X_reduced[:, 2],
	cmap=plt.cm.Set1,
	edgecolor="k",
	s=40,
	)

	ax.set_title(title)
	ax.set_xlabel("1st eigenvector")
	ax.xaxis.set_ticklabels([])
	ax.set_ylabel("2nd eigenvector")
	ax.yaxis.set_ticklabels([])
	ax.set_zlabel("3rd eigenvector")
	ax.zaxis.set_ticklabels([])

	plt.savefig(title+".pca3D.png")
	plt.clf()



def run_louvain(df):
	''' run louvain clustering on my dataframe  '''
	print('louvain:')
	sdf = df.astype(pd.SparseDtype('int',0))
	cols = list(df)
	rows = list(df.index)

	print('sparse',sdf.head())
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

def make_overlap_edges(reads_df,k):
	''' find reads that overlap methyl positions '''
	# TODO:
	# consider dataframe of positions, rows are reads, where alignments are already lined up

	links = []
	for a,b in itertools.combinations(reads_df,2):
		# print('a',reads_df[a],'b',reads_df[b])
		
		if get_chrom(reads_df[a])==get_chrom(reads_df[b]):
			a_pos = sorted(list(reads_df[a][get_chrom(reads_df[a])].keys()))
			b_pos = sorted(list(reads_df[b][get_chrom(reads_df[b])].keys()))
			# print('a',a_pos[:10],reads_df[a][get_chrom(reads_df[a])][a_pos[0]])
			# print('b',b_pos[:10],reads_df[b][get_chrom(reads_df[b])][b_pos[0]])
			# check if reads are aligned in the same orientation
			if reads_df[a][get_chrom(reads_df[a])][a_pos[0]][-1] == reads_df[b][get_chrom(reads_df[b])][b_pos[0]][-1]:
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
					print(a,b,'match_count:',match_count, match_count/len(a_pos[overlap_start:overlap_end]))

				elif (b_pos[0] < a_pos[0]) and (b_pos[-1] > a_pos[0]):
					overlap_start = find_overlap_index(b_pos,a_pos)
					everlap_end   = find_overlap_end_index(b_pos,a_pos)
					print('a over b: overlap_start',overlap_start,len(b_pos[overlap_start:overlap_end]),a,b)
					# if exact_overlap_k(a_pos[overlap_start::],b_pos,len(a_pos[overlap_start:overlap_end])):
					print(exact_overlap_k(b_pos[overlap_start::],a_pos,len(b_pos[overlap_start:overlap_end])))
					print(hamming_overlap_k(b_pos[overlap_start:overlap_end],a_pos,len(b_pos[overlap_start:overlap_end])))
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
					print(a,b,'match_count:',match_count, match_count/len(a_pos[overlap_start:overlap_end]))
					# print(b_pos[overlap_start::],'\n',a_pos)
							# break


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



def main(bamPath,min_prob,min_reads,sample_name):

	# set ml prob cutoff
	# if not (min_prob):
	# 	min_prob=0.5
	ml_min_prob = min_prob
	ml_cutoff = ml_min_prob * 255

	# dictionary to store reads and modified bases in those reads
	# mod_base_in_read['read_name'] = {base_pos_in_seq : [mm_index, ml_prob]}
	mod_base_in_read = {}
	
	# run methyl analyze per read using url
	mod_base_in_read = get_methyl_per_read(bamPath,ml_cutoff,min_reads, mod_base_in_read,sample_name)
	# mod_base_in_read = get_methyl_per_read("http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam",ml_cutoff, mod_base_in_read)
	# mod_base_in_read = get_methyl_per_read("HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.chr20_run2.haplotagged.bam",ml_cutoff, mod_base_in_read)

	# plot_modbase_locations(mod_base_in_read)

	make_overlap_edges(mod_base_in_read, 10)






if __name__ == "__main__":
	'''
	Python script to cluster reads in a bam by their methyl tag provile

	Usage: python3.9 analyze_methylBam.py -b http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam
	-b bam 
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
		default=3,
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


