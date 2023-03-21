#!/usr/bin/env python3
import argparse
import os
import pysam
from Bio.Seq import Seq
import re
import numpy as np


#TODO: 
'''	test the function
	get reference position instead of read position
	get mod probabilities for (-) read'''

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

def identify_seq_orientation(read):
	'''
		Mm tags are created as the read comes off the sequencer, always starting at the 5' end.
		So the orientation of the tags and the seq fies may not be the same if the seq field 
		was reverse complemented to align to the reference sequence. 

	'''
	if read.is_forward:

		# ? should this be read.query_sequence, query_alignment_sequence, query_alignment_start ? 
		# seq = Seq(read.seq)
		seq = read.seq
		print('\n +', read.reference_name, read.reference_start, read.query_name,'seq',seq[:100] )
		print('\t',read.cigarstring[:100]) 
		print('\t',read.get_tag('Mm')[:100])
	#TODO
	elif read.is_reverse:
		seq = read.seq #get_forward_sequence() #seq #Seq(read.seq).reverse_complement()
		print('\n -',read.reference_name, read.reference_start, read.query_name,'\nread.seq',read.seq[:100],'\n seq',seq[:100],'...',seq[-120:] )
		print('\t',read.cigarstring)
		print('\t',read.get_tag('Mm')[:100]) 
	else: 
		print('whats the problem?', read)
	return seq

def parse_mm_tag(read):
	'''mm tag needs to be split and skip bases list isolted: Mm C+m?,1,1,6,5;
		ml tag is stored as an array from the get_tag method, mm is not	'''
	mm_arr = read.get_tag('Mm')
	modbase = mm_arr[0]
	# semi-colon delimiter at the end of the list separating the tags is removed 
	mm = mm_arr.split(',')[1:-1]
	return modbase, mm

def get_cigar_ref_pos(read):
	''' parse the cigar string to obtain methylated cytosine positions
		in refernce coordinates for each read '''

	base_counter=0
	pos=read.reference_start+1
	readcigarpositions = [read.reference_start+1]

	# for anti-sense strand read alignments the cigar tuples need to be read from back to front, 
	# they are created with respect to the reference. 
	if read.is_reverse:
		cigartups = read.cigartuples[::-1]
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
				base_counter+=c[1]
			else:
				readcigarpositions.extend([0]*c[1])
		# 1 : INS; add place holder ref coordinate for each read only base consumed by the insertion
		if c[0] in [1]:
			readcigarpositions.extend([pos]*c[1])
			base_counter+=c[1]


		# 0 : M :  match/mismatch add ref positions to the array for every M base
		# 7 : = :  add ref positions to the array for every = base
		# 8 : X :  mismatch add ref position as if it was =
		if c[0] in [0,7,8]:
			readcigarpositions.extend( [i for i in np.arange(pos,pos+c[1]) ] )
			pos+=c[1] 
			base_counter+=c[1]

		# 2 : D :  DEL;skip the ref position
		# 3 : N : REF_SKIP; 
		if c[0] in [2,3]:
			pos+=c[1]

	return readcigarpositions

def write_read_mods_to_bed(mod_base_in_read):
	''' write all read mods to a bed file per read'''
	for read,chroms in mod_base_in_read.items():

		bfile = open(str(read)+'.mods.bed','w')

		for chrom , refPos in chroms.items():
			for pos, modinfo in refPos.items():
				bfile.write(str(chrom)+'\t'+str(pos)+'\t'+str(pos+1)+'\t'+str(modinfo[1])+'\n' )
		bfile.close()

def get_methyl_per_read(filePath,ml_cutoff, mod_base_in_read):
	''' iterate over all of the read mapping to a specified region using fetch(). 
	Each iteration returns a AlignedSegment object which represents a single read along with its fields and optional tags:
	'''
	# open alignment file
	samfile = pysam.AlignmentFile(filePath, "rb")

	# Fetch reads from region: chr20:1000000-1020000
	# neuron loc chr20:63664503-63664830
	read_count = 0
	mmTag, mlTag = '', ''

	for read in samfile.fetch('chr20', 1000000, 1000020):
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

		#TODO: reverse strand seqs are messed up either in seq or in CIGAR string
		# store the seq in the proper orientation relative to the Ml tag
		seq = identify_seq_orientation(read)

		# store the index position of every modified base in the read
		mod_base_i = [b.start() for b in re.finditer(modbase.upper(),str(seq).upper())]

		# get reference positions that account for INDELs/Clips by parsing the cigar string
		readcigarpositions = get_cigar_ref_pos(read)

		if read.is_reverse:
			for i,s in enumerate(seq[:250]):
				print(s,readcigarpositions[i])

		# move through the read by moving in the increments reported by the mm tag
		read_modbase_pos = 0
		for i, i_skip in enumerate(mm):
			# increment the read position by the number of modbases indicies to skip 
			# ( the number of C's to skip to land on the C that is observed as modified )
			read_modbase_pos = read_modbase_pos + int(i_skip)
			
			# it is the modbase after the # to skip
			if i>0:
				read_modbase_pos+=1

			# only count alignments after the read.query_alignment_start position
			if mod_base_i[read_modbase_pos] > read.query_alignment_start and mod_base_i[read_modbase_pos] < read.query_alignment_end:
				##############################
				if i < 20:
					print(i,i_skip,'read_modbase_pos',read_modbase_pos, 'modbase_i',mod_base_i[read_modbase_pos],
						seq[ mod_base_i[read_modbase_pos]-10:mod_base_i[read_modbase_pos]+10 ].upper(),
						read.positions[ mod_base_i[read_modbase_pos]-read.query_alignment_start ], readcigarpositions[mod_base_i[read_modbase_pos]] )
					print('\t base is modbase?:',seq[ mod_base_i[read_modbase_pos] ].upper() == modbase.upper() )
				##############################
				# sanity check that that base is the modbase from the Mm Tag
				if seq[ mod_base_i[read_modbase_pos] ].upper() != modbase.upper():
					print('base is not modbase',seq[:read_modbase_pos].upper(), 'at',i,'. i_skip:', i_skip)
					break

				# if the ml probability of being modified is > cutoff, 
				# store the modified base location as it 0-based position in the read sequence string
				# read mod base dict : {readName : ref chrom : ref position : [positon in mm tag, ml value, position in read]}
				# TODO Store the strand
				if ml[i] >= ml_cutoff:
					if read.query_name not in mod_base_in_read.keys():
						mod_base_in_read[read.query_name]={read.reference_name : {readcigarpositions[mod_base_i[read_modbase_pos]] : [ i,ml[i],mod_base_i[read_modbase_pos] ]} }
					if readcigarpositions[mod_base_i[read_modbase_pos]] not in mod_base_in_read[read.query_name][read.reference_name].keys():
						mod_base_in_read[read.query_name][read.reference_name][readcigarpositions[mod_base_i[read_modbase_pos]]] = [ i,ml[i],mod_base_i[read_modbase_pos] ]
			else:
				if i < 20:
					print('outside of clip:', mod_base_i[read_modbase_pos],i,i_skip,'read_modbase_pos',read_modbase_pos)


		# if read_count>2:
		# 	break


	samfile.close()

	write_read_mods_to_bed(mod_base_in_read)

	return mod_base_in_read

def main():

	# set ml prob cutoff
	ml_min_prob = 0.5
	ml_cutoff = ml_min_prob * 255

	# dictionary to store reads and modified bases in those reads
	# mod_base_in_read['read_name'] = {base_pos_in_seq : [mm_index, ml_prob]}
	mod_base_in_read = {}
	
	# run methyl analyze per read using url
	mod_base_in_read = get_methyl_per_read("http://public.gi.ucsc.edu/~memeredith/methyl/HG002_card/HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.haplotagged.bam",ml_cutoff, mod_base_in_read)
	# mod_base_in_read = get_methyl_per_read("HG002_ONT_card_2_GRCh38_PEPPER_Margin_DeepVariant.chr20_run2.haplotagged.bam",ml_cutoff, mod_base_in_read)

if __name__ == "__main__":
	'''
	Python script to cluster reads in a bam by their methyl tag provile

	Usage: python3.9 analyze_methylBam.py 
	-y hap1.yak.switch-error.mat.sorted.bed -r hap2.yak.switch-error.txt.mat.grch38.sorted.bed -o merged_lifted_files/
	'''
	main()
	
	# parser = argparse.ArgumentParser()

	# parser.add_argument(
	# 	"-y",
	# 	required=True,
	# 	type=str,
	# 	help="Input yak altered bed, sorted by the component ID "
	# )

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

	# args = parser.parse_args()

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





















