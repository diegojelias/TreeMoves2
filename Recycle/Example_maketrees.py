
#!/usr/bin/env python
############################## 
### NEED TO REWRITE INTO FUNCTIONS 
# Need to change lots of details below. Then whole script can be run, or pieces can be run in python interactively.
# This script calls makeNexus.sh, which requires user to input a custom path to Paup* 
# This script is all based on unweighted RF distances, will need to edit normalized RF calculations for weighted RF I think. 

############################## 

from __future__ import (division, print_function)
from readTree import Node
from readTree import Tree
import readTree
import dendropy
import subprocess
from cStringIO import StringIO
from Bio import Phylo

############################################################  
#User input
############################################################

# Number of tips = 10
tips = 75

# Total out trees per cloud
cloud_size = 1000

# Desired normalized RF distance between starting trees and for each cloud of trees
RF_norm_start = 1.0
RF_norm_cloud = 0.125

# Title run
name =("%stip_%strees_%s_%sstart_taco" % (tips,cloud_size,RF_norm_cloud,RF_norm_start))

############################################################  
#Choose starting trees
############################################################
print("Making starting trees...")
# Create random tree. Average branch length of 1. No variation in branch length. Don't print out tree structure
t1 = readTree.rand_tree(tips=tips,brl_avg=1,brl_std=None,verbose='F')

# Calculate number of NNI moves based on desired normalized RF distance.
RF_max = 2*(tips-2)
NNI_moves_start = int((RF_max * RF_norm_start)/2)
if NNI_moves_start == 0:
	NNI_moves_start = 1
# Create second starting tree
t2 = readTree.NNI_mult_moves(in_tree=t1,num_moves=NNI_moves_start,node_choice='random',no_dup_start_tree='F', req_min_RF=RF_norm_start)

# Write out tree files
readTree.write_single_tree(t1,'%s_starting_tree_1.tree' % name)
readTree.write_single_tree(t2,'%s_starting_tree_2.tree' % name)

############################################################  
# Make cluster of trees around each starting tree
############################################################
# Calculate number of NNI moves based on desired normalized RF distance.
NNI_moves_cloud = int((RF_max * RF_norm_cloud)/2)
if NNI_moves_cloud == 0:
	NNI_moves_cloud = 1
# Make a cloud for each starting tree
c_size = int(cloud_size)-1 

print("Making clouds...")

# Make clouds
cluster1 = readTree.NNI_mult_trees(in_tree=t1,num_out_trees=c_size,num_nni_moves=NNI_moves_cloud,out='list')
cluster2 = readTree.NNI_mult_trees(in_tree=t2,num_out_trees=c_size,num_nni_moves=NNI_moves_cloud,out='list')
# Make a nexus file with starting trees and cloud trees
readTree.list_to_out(cluster1, cluster2,'%s_cloud.tree' % name)

############################################################  
# Print info to log file
############################################################
print("Calculating stats on trees...")
# Calculate emperical distance between two start trees
RF_emp = readTree.rf_unweighted(t1,t2,normalized='T')[1]
# Calculate average density of each cloud
RF_cloud1 = readTree.cluster_density_avg(in_file='%s_cloud.tree' % name, NNI_trees= int(cloud_size)-1, starting_tree_number=0)
RF_cloud2 = readTree.cluster_density_avg(in_file='%s_cloud.tree' % name, NNI_trees= int(cloud_size)-1, starting_tree_number=cloud_size)
# write to log file
with open('%s.log' % name, "w") as log_file:
	line1 = "File name: "+str(name)
	line1b = "Tips: "+str(tips)+", Trees per cloud: "+str(cloud_size)
	line2 = "Starting trees - RF_input: "+str(RF_norm_start)+", RF_calc: "+str(RF_emp)+", NNI_moves: "+str(NNI_moves_start)
	line3 = "Cloud of trees - RF_input: "+str(RF_norm_cloud)+", RF_calc1: "+str(RF_cloud1)+", RF_calc1: "+str(RF_cloud2)+", NNI_moves: "+str(NNI_moves_cloud)
	log_file.write("%s\n%s\n%s" % (line1, line2, line3))

print("Making into nexus files...")
# Turn script into Nexus that is readable by TreeScaper
# This script turns all *.tree files into *.nex files 
subprocess.call(['./makeNexus.sh'])
print("Done")

