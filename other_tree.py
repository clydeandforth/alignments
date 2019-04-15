#!/home/b.bss81c/miniconda2/bin/python

# Modules to build the tree
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio.Phylo import draw
from Bio import Phylo, AlignIO
from Bio.Phylo.Consensus import *
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import pylab
import sys

alignment = AlignIO.read('gyrB.fa', 'fasta') # reading the alignment file
msas = bootstrap(alignment, 100)
calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(alignment) # distance matrix

constructor = DistanceTreeConstructor(calculator)

#tree = constructor.nj(dm) # build with neighbour joining algorithm a tree from dm
tree = constructor.upgma(dm)
#tree = bootstrap_tree(msa, 50, constructor)
trees = bootstrap_trees(dm, 100, constructor)
#tree.root_at_midpoint()

consensus_tree = bootstrap_consensus(alignment, 100, constructor, majority_consensus)
consensus_tree.ladderize()
consensus_tree.root.color="green"
mrca = consensus_tree.common_ancestor({"name": "PC_00004"}, {"name": "BG_I_00594"})
mrca.color = "salmon"

plt.xlabel("Branch length")
plt.ylabel("Taxa")

Phylo.write(consensus_tree,  'TreeToCutOff.xml', 'phyloxml')

plt.rc('font', size=10)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
plt.rc('axes', titlesize=14)     # fontsize of the axes title
#plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
#plt.rc('figure', titlesize=18)   # fontsize of the figure title


#draw(tree, do_show=False)
#plt.savefig("TreeToCutOff_check.svg", format='svg', dpi=1200, bbox_inches='tight')
#pylab.savefig("phylo-dot.png")

Phylo.draw(consensus_tree,  show_confidence=True)
pylab.savefig("phylo-dot.png", dpi=300)
#plt.show()
#aln = AlignIO.read(open('gyrB.phy'), 'phylip')
#starting_tree = Phylo.read('TreeToCutOff.nwk', 'newick')
#scorer = ParsimonyScorer()
#searcher = NNITreeSearcher(scorer)
#constructor = ParsimonyTreeConstructor(searcher, starting_tree)
#pars_tree = constructor.build_tree(aln)
#print pars_tree
