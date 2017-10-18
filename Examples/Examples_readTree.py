
#!/usr/bin/env python
'''
To Do:
'''

-------------------------------------------------------
exit()



python


from __future__ import (division, print_function)
from readTree import Node
from readTree import Tree
import readTree
import dendropy
from cStringIO import StringIO
from Bio import Phylo

############################################################ 
#Scratch paper
############################################################



###############################################################################
#Examples of working parts
###############################################################################
#Notes
'''
Needs to be updated with new options written into readTree.py
'''


############################## 
#Input tree options
############################## 
# From string
doo="(P:0.09,(Q:0.07,(X:0.02,((Y:0.03,Z:0.01):0.02,W:0.08):0.06):0.03):0.04)"
boo="(P:0.01,(Q:0.01,(X:0.01,((Y:0.01,Z:0.01):0.01,W:0.01):0.01):0.01):0.01)"
# Create readTree.tree object
in_tree_object = Tree(doo)

# From file
in_tree_object=readTree.read_nexus('oneTree.t',0)

# Create random starting tree
in_tree_object = readTree.rand_tree(tips=10,brl_avg=1,brl_std=None,verbose='T')

############################## 
# Make moves
############################## 
# Print newick string
in_tree_object.newick(in_tree_object.root)

# Make NNI move 
new_tree_object = readTree.NNI(orig_tree=in_tree_object,node_choice='exponential')

# Make multiple moves on a single tree
new_tree_object = readTree.NNI_mult_moves(in_tree=in_tree,num_moves=1,node_choice='exponential',no_dup_start_tree='T')

# Make multiple trees, each with one NNI move from starting tree. Output as nexus file.
readTree.NNI_mult_trees(in_tree=in_tree,num_out_trees=10,num_nni_moves=2,out_file='outFile2.t',node_choice='random',no_dup_start_tree='T')

############################## 
# Look at and compare trees
############################## 

# Print newick string
new_tree_object.newick(new_tree_object.root)

# View tree in terminal along with newick string
readTree.view_phylo(new_tree_object)

# Check that all trees have the right number of moves from start tree
readTree.compare_tree_file(in_file='outFile1.t',total_trees=11,distance_metric="uRF")

# Compare RF distance between two trees.
readTree.rf_unweighted(in_tree_object,new_tree_object,normalized='T')
