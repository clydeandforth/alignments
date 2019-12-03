#!/scratch/b.bss81c/anaconda3/envs/ETE-toolkit/bin/python

# Modules to build the tree
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio.Phylo import draw
from Bio import Phylo, AlignIO
from Bio.Phylo.Consensus import *
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.PAML import codeml
from ete3 import Tree, TreeStyle, faces
from ete3 import PhyloTree
from Bio import Nexus
import subprocess
import matplotlib

import matplotlib.pyplot as plt
import networkx
import pylab
import sys



def ML_tree():
        ete3 build -n concat_padded.fasta --cogs temp3.txt -o sptree1_results -m sptree_fasttree_100 -w standard_raxml_bootstrap --clearall --noimg
# 	tree = PhyloTree("/scratch/b.bss81c/alignments/sptree1_results/clustalo_default-none-none-raxml_default/outfile_stripped.fasta.final_tree.nw")
#	ts = TreeStyle()
#	tree.render("formatted.png", tree_style = ts, h=3000, w=3000)

ML_tree()
