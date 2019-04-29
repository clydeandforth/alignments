#!/scratch/b.bss81c/anaconda3/bin/python
#James Doonan
#arguments complete_script.py infile.gbk outfile.fasta 'gene of interest'


from Bio import SeqIO, Seq, AlignIO
import sys
import os
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio.Phylo import draw
from Bio import Phylo, AlignIO
from Bio.Phylo.Consensus import *
from Bio.Align.Applications import ClustalwCommandline
import subprocess
import matplotlib

import matplotlib.pyplot as plt
import networkx
import pylab

product = sys.argv[3]

def gene_find():
	with open(sys.argv[2], 'w') as t:
	        for rec in SeqIO.parse(sys.argv[1], "genbank"):
        	        if rec.features:
                	        for feature in rec.features:
                        	        if feature.type == "CDS":
                                	        if product in feature.qualifiers['product'][0]:
                                        	        t.write(
                                                	">%s|%s %s \n%s\n" % (
                                              		feature.qualifiers['locus_tag'][0],
                                                	feature.qualifiers['gene'][0],
                                                	rec.description,
                                                        #rec.name,
                                                	feature.location.extract(rec).seq))

gene_find()

input_file = 'outfile.fasta'

def make_equal():

	records = (rec.upper() for rec in SeqIO.parse(input_file, 'fasta'))
	records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen

	maxlen = max(len(record.seq) for record in records)

	# pad sequences so that they all have the same length
	for record in records:
    		if len(record.seq) != maxlen:
        		sequence = str(record.seq).ljust(maxlen, '.')
        		record.seq = Seq.Seq(sequence)
	assert all(len(record.seq) == maxlen for record in records)
	output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
	with open(output_file, 'w') as f:
    		SeqIO.write(records, f, 'fasta')

make_equal()

def align_seq():
        cline = ClustalwCommandline("clustalw2", infile="outfile_padded.fasta")
        stdout = cline()
align_seq()

def draw_tree():
        alignment = AlignIO.read('outfile_padded.aln', 'clustal') # reading the alignment file
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        msas = bootstrap(alignment, 100)
        calculator = DistanceCalculator('blosum62')

        constructor = DistanceTreeConstructor(calculator)

        trees = bootstrap_trees(alignment, 100, constructor)

        consensus_tree = bootstrap_consensus(alignment, 1000, constructor, majority_consensus)
        consensus_tree.ladderize()
        consensus_tree.root.color="green"
        #mrca = tree.common_ancestor({"name": "PC_00004"}, {"name": "BG_I_00594"})
        mrca = consensus_tree.common_ancestor({"name": "PC_00004|DNA"})
        mrca.color = "salmon"


        Phylo.write(consensus_tree,  'TreeToCutOff.xml', 'phyloxml')

        #plt.rc('font', size=10)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
        #plt.rc('axes', titlesize=14)     # fontsize of the axes title
        #plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
        #plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
        #plt.rc('figure', titlesize=18)   # fontsize of the figure title

#plt.savefig("TreeToCutOff_check.svg", format='svg', dpi=1200, bbox_inches='tight')


        Phylo.draw(consensus_tree,  show_confidence=True)
        pylab.gcf().set_dpi(300)
        pylab.savefig("phylo-dot.png")
        pylab.clf()
draw_tree()

