#!/usr/bin/env python3
import pysam
import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd
import itertools
from methlotype.utils import bam_functions as bf
from methlotype.utils import graph_utils as graphing

# import bam_functions as bf
# import graph_utils as graphing

class MethylBam:
    """ Class for storing and parsing a BAM file that contains methyl tags from ONT sequencing """
    def __init__(self, bamPath, ml_cutoff, min_reads, sample_name, region):
        self.bamPath = bamPath
        self.ml_cutoff = ml_cutoff
        self.min_reads = min_reads
        self.sample_name = sample_name
        # this won't always be run on a region
        self.region = region
        self.mod_base_in_read = {}
        self.modbase = None
        self.mm = None

    def strand(self, x):
        """store strand as +/- instead of boolean"""
        return '-' if x else '+'

    def get_methyl_per_read(self):
        """ iterate over all the read mapping to a specified region using fetch().
        Each iteration returns a AlignedSegment a single read object with its fields and optional tags:
        """
        # open alignment file
        samfile = pysam.AlignmentFile(self.bamPath, "rb")
        logfile = open('log.txt', 'w')

        read_set = set()
        read_strands = {}
        chrom_modbase_pos = {}

        # Fetch reads from region: chr startPos endPos
        mmTag, mlTag = '', ''
        if len(self.region) == 3:
            chromosome = self.region[0]
            start = int(self.region[1])
            end = int(self.region[2])
        else:
            return -1

        # TODO This won't always be per region....
        # TODO maybe filter out short read alignments, of some cutoff
        for read in samfile.fetch(chromosome, start, end):
            # store some easy access metadata
            read_set.add(read.query_name)
            bfile = open('beds/' + str(read.query_name) + '.mods.bed', 'w')
            logfile.write(' '.join([ '\n', str(self.strand(read.is_reverse)), str(read.mapq), str(read.reference_name),
                                     str(read.reference_start), str(read.reference_end), str(read.query_name) ]) )

            # get mm and ml tag information per read
            # the methyl tags can be either Mm/Ml or MM/ML for R10 or R9
            if len(mmTag) == 0:
                mm_tag, ml_tag = bf.identify_methyl_tags(read)
                if ml_tag == -1:
                    return -1
            # mm tag needs to be split and skip bases list isolated:
            self.parse_mm_tag(read, mm_tag)
            ml = read.get_tag(ml_tag)

            # store the seq in the proper orientation relative to the Ml tag
            seq = read.seq  # seq = identify_seq_orientation(read)

            # store the index position of every base of interest in the read
            mod_base_i = [b.start() for b in re.finditer(self.modbase.upper(), str(seq).upper())]

            # get reference positions that account for INDELs/Clips by parsing the cigar string
            refcigarpositions = bf.get_cigar_ref_pos(read)

            # store ref start and end positions per read for easy access
            read_strands[read.query_name] = [self.strand(read.is_reverse), bf.get_hp_tag(read),
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
                        if int(ml[i]) >= self.ml_cutoff:
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
                        if int(ml[i]) >= self.ml_cutoff:
                            pos = refcigarpositions[mod_base_i[read_modbase_pos]] + 1  #CHANGEd

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
            print('chrom', chrom)
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

        bf.make_quick_hm(df, self.sample_name)
        # make_pca(df, self.sample_name)

        # select only the region, instead of the flanks with low coverage artifact of pysam fetch
        df_region_only = df.loc[:, start:end]
        bf.make_quick_hm(df_region_only, self.sample_name + '.reads.region')
        # make_pca(df_region_only, self.sample_name + '.region')

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
            self.modbase = bf.rev_comp(self.modbase)
        # semicolon delimiter at the end of the list separating the tags is removed
        self.mm = mm_arr[:-1].split(',')[1::]

    def check_modbase_rev(self, seq, mod_base_i, mod_base_i_index):
        # sanity check that that base is the modbase from the Mm Tag
        if seq[mod_base_i[-mod_base_i_index] - 1] != bf.rev_comp(self.modbase).upper():
            print(seq, '\nnot cpg', seq[mod_base_i[-mod_base_i_index] - 1], mod_base_i[-mod_base_i_index], mod_base_i)
        if seq[mod_base_i[-mod_base_i_index]].upper() != self.modbase.upper():
            print('base is not modbase', seq[:-mod_base_i_index].upper())
            return -1

    def check_modbase(self, seq, mod_base_i, read_modbase_pos):
        # sanity check that that base is the modbase from the Mm Tag
        if seq[mod_base_i[read_modbase_pos] + 1].upper() != bf.rev_comp(self.modbase).upper():
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

