# Original Author: Jean Monlong
# 2023-06-27 Melissa Meredith
import numpy
import random
import pysam
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# array of nucleotides
_nuc = numpy.array(["A", "T", "C", "G"])


def randSeq(length):
    ''' fill an array with length number of nucleotides by selecting a
        random integer between 0-2 to select a base from the _nuc array '''
    seqArray = _nuc[[random.randint(0, 3) for i in range(length)]]
    return(MutableSeq("".join(seqArray)))


# parameters
N = 100000
snp_rate = .01
recomb_rate = .001
N_ont_reads = 100
ont_read_size = 11000

# "ancestral" sequence
anc = randSeq(N)
# write as "reference" chromosome 1
SeqIO.write([SeqRecord(anc, id='chr1')], "ref.fa", "fasta")

# 2 sets of parent haplotypes with some SNPs
parents = [[], []]
for parent in range(2):
    for haplotype in range(2):
        # MutableSeq allows for editing the ancestral seq in place
        par = MutableSeq(str(anc))
        # for each base in the parent sequence add SNPs based on snp_rate
        for pos in range(len(par)):
            if random.random() < snp_rate:
                new_base = _nuc[random.randint(0, 3)]
                # keep selecting new bases until it is different than the anc base
                while par[pos] == new_base:
                    new_base = _nuc[random.randint(0, 3)]
                par[pos] = new_base
        # store that seq as a parent haplotype
        parents[parent].append(par)

# recombine parents into child haplotype pair
haps = ['', '']
for parent in range(2):
    cur_par = 0
    for pos in range(len(parents[parent][0])):
        haps[parent] += (parents[parent][cur_par][pos])
        # switch parents based on recombination_rate per base
        if random.random() < recomb_rate:
            cur_par = 1 - cur_par

# simulate ONT reads
ont_f = open('test_ont.fastq', 'wt')
readid = 0
for rr in range(N_ont_reads):
    hap = random.randint(0, 1)
    pos = random.randint(0, N - ont_read_size)
    read = haps[hap][pos:(pos+ont_read_size)]
    ont_f.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*ont_read_size))
    readid += 1
ont_f.close()

# write out a bam that includes mm & ml tags
header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'}] }
                   # {'LN': 1584, 'SN': 'chr2'}

#read_28833_29006_6945
with pysam.AlignmentFile('tiny.rev.bam', "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6945"
    a.query_sequence="CGGCCAAGACCAAGATATACGCGTAGCTACGTAAGCT"
    #                "CG-CCAAGACCAAGATATACG-G-AGCTACG-AAGCT"-
    a.flag = 16
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((7, 10), (2, 1), (7, 10), (8, 5), (7, 5), (1, 4),(7, 3))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<<<")
    a.tags = (("NM", 1),
              ("MM", "C+m?,1,1,0,3;"),
              ("ML", (231, 200, 146, 220)),
              ("HP", 1))
    outf.write(a)

with pysam.AlignmentFile('tiny.rev.2.bam', "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6946"
    a.query_sequence="ACGCGCAAGCGAAAGATATACGCGTCGATACGTACGATCGCTTACG"
                    #"AC-C-CAAGC-AAAGATATACG-G-CGATACG-ACGATCGCTTACG"
                    #  * *    s*    s      * * *     *   s  *      s
    a.flag = 16
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((7, 10), (2, 1), (7, 13), (8, 5), (7, 5), (1, 4), (7, 9))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<<<")
    a.tags = (("NM", 1),
              ("MM", "C+m?,1,1,0,0,0,1,1,0;"),
              ("ML", (231, 200, 146, 220, 231, 200, 146, 220)),
              ("HP", 1))
    outf.write(a)

with pysam.AlignmentFile('tiny.rev.3.bam', "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6946_2"
    a.query_sequence="ACGCGCAAGCGAAAGATATACGCGTCGATACGTACGATCGCTTACG"
                    #"AC-C-CAAGC-AAAGATATACG-G-CGATACG-ACGATCGCTTACG"
                    #  * *    s*    s      * * *     *   s  *      s
    a.flag = 16
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((7, 10), (2, 1), (7, 13), (8, 5), (7, 5), (1, 4), (7, 9))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<<<")
    a.tags = (("NM", 1),
              ("MM", "C+m?,1,1,0,0,0,1,1,0;"),
              ("ML", (231, 200, 146, 220, 231, 200, 146, 220)),
              ("HP", 1))
    outf.write(a)

#read_28833_29006_6954
with pysam.AlignmentFile('tiny.bam', "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6954"
    a.query_sequence="AGCTTACGTAGCTACGCGTATATCTTGGTCTTGGCCG"
    #.query_sequence="AGCTT-CGTAGCT-C-CGTATATCTTGGTCTTGG-CG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((4, 5), (7, 5), (2, 1), (7, 10), (8, 5), (7, 5), (1, 4),(7, 3))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<<<")
    a.tags = (("NM", 1),
              ("MM", "C+m?,1,1,0,3;"),
              ("ML", (231, 200, 146, 220)),
              ("HP", 0))
    outf.write(a)


#read_28833_29006_6954_rev
with pysam.AlignmentFile('tiny.r.bam', "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6954_r"
    a.query_sequence="AGCTTACGTAGCTACGCGTATATCTTGGTCTTGGCCG"
    #.query_sequence="AGCTT-CGTAGCT-C-CGTATATCTTGGTCTTGG-CG"
    a.flag = 16
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((4, 5), (7, 5), (2, 1), (7, 10), (8, 5), (7, 5), (1, 4),(7, 3))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<<<")
    a.tags = (("NM", 1),
              ("MM", "C+m?,0,4,0,1;"),
              ("ML", (231, 200, 146, 220)),
              ("HP", 0))
    outf.write(a)

with pysam.AlignmentFile('tiny.dup.bam', "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6954_2"
    a.query_sequence="AGCTTACGTAGCTACGCGTATATCTTGGTCTTGGCCG"
    #.query_sequence="AGCTT-CGTAGCT-C-CGTATATCTTGGTCTTGG-CG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((4, 5), (7, 5), (2, 1), (7, 10), (8, 5), (7, 5), (1, 4),(7, 3))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<<<")
    a.tags = (("NM", 1),
              ("MM", "C+m?,1,1,0,3;"),
              ("ML", (231, 200, 146, 220)),
              ("HP", 0))
    outf.write(a)

'''
samtools merge -f -o tiny.merged.bam tiny.bam tiny.rev.bam
samtools index tiny.merged.bam
 '''
# pysam.sort("-o", "output.bam", "ex1.bam")

pysam.merge("-f", "-o", "tiny.merged.2.bam", "tiny.bam", "tiny.dup.bam", 'tiny.r.bam',
                        "tiny.rev.bam")
pysam.index("tiny.merged.2.bam")

pysam.merge("-f", "-o", "tiny.merged.6.bam", "tiny.bam", "tiny.dup.bam",
                        "tiny.rev.bam", 'tiny.rev.2.bam', 'tiny.rev.3.bam')
pysam.index("tiny.merged.6.bam")

pysam.merge("-f", "-o", "tiny.merged.identical.bam", 'tiny.rev.2.bam', 'tiny.rev.3.bam')
pysam.index("tiny.merged.identical.bam")