
#!/usr/bin/env python
############################## 
# Edit main() as needed
# There are a bunch of other functions from a pre-evolution scramble to make data
# This script calls makeNexus.sh, which requires user to input a custom path to Paup* 

############################## 

from __future__ import (division, print_function)
from readTree import Node
from readTree import Tree
import readTree
import dendropy
import subprocess
from cStringIO import StringIO
from Bio import Phylo
import os
import glob
import shutil 
import itertools
import re 


def makeStarting(tips,RF_norm_start,name,RF_norm_cloud,cloud_size):
	#Choose starting trees
	print(name)
	print("Making starting trees...")
	# Create random tree. Average branch length of 1. No variation in branch length. Don't print out tree structure
	t1 = readTree.rand_tree(tips=tips,brl_avg=1,brl_std=None,verbose='F')

	# Calculate number of NNI moves based on desired normalized RF distance.
	RF_max = 2*(tips-2)
	NNI_moves_start = int((RF_max * RF_norm_start)/2)
	if NNI_moves_start == 0:
		NNI_moves_start = 1
	# Calc cloud for log file
	NNI_moves_cloud = int((RF_max * RF_norm_cloud)/2)
	# Create second starting tree
	t2 = readTree.NNI_mult_moves(in_tree=t1,num_moves=NNI_moves_start,node_choice='random',no_dup_start_tree='F', req_min_RF=RF_norm_start)
	readTree.rf_unweighted(t1,t2,"T")
	t3 = readTree.NNI_mult_moves(in_tree=t1,num_moves=NNI_moves_start,node_choice='random',no_dup_start_tree='F', req_min_RF=RF_norm_start)
	d12=readTree.rf_unweighted(t1,t2,"T")[1]
	d23=readTree.rf_unweighted(t2,t3,"T")[1]
	d13=readTree.rf_unweighted(t1,t3,"T")[1]
	print("first t2 to t3 RF"+str(d23))
	#t1 and t2 are guarenteed to be RF start apart, same with t2 and t3. So we only need to check t1 and t3.
	while d23 != RF_norm_start:
		print("redoing third start tree")
		t3 = readTree.NNI_mult_moves(in_tree=t1,num_moves=NNI_moves_start,node_choice='random',no_dup_start_tree='F', req_min_RF=RF_norm_start)
		d23=readTree.rf_unweighted(t2,t3,"T")[1]
		print(d23)

	# Write out tree files
	readTree.write_single_tree(t1,'%s_starting_tree_1.tree' % name)
	readTree.write_single_tree(t2,'%s_starting_tree_2.tree' % name)
	readTree.write_single_tree(t3,'%s_starting_tree_3.tree' % name)

	# Write out log file
	rfs=str(d12)+", "+str(d13)+", "+str(d23)
	with open('%s.log' % name, "w") as log_file:
		line1 = "File name: "+str(name)
		line1b = "Tips: "+str(tips)+", Trees per cloud: "+str(cloud_size)
		line2 = "Starting trees - RF_input: "+str(RF_norm_start)+", RF_calc: "+str(rfs)+", NNI_moves: "+str(NNI_moves_start)
		line3 = "Cloud of trees - RF_input: "+str(RF_norm_cloud)+", NNI_moves: "+str(NNI_moves_cloud)
		log_file.write("%s\n%s\n%s" % (line1, line2, line3))
	# Pass tree files to next function
	return t1,t2,t3


def makeCloud(tips,RF_norm_cloud,name,cloud_size,starting_trees): 
	# Make cluster of trees around each starting tree
	# Get starting trees parsed
	t1=starting_trees[0]
	t2=starting_trees[1]
	t3=starting_trees[2]
	# Calculate number of NNI moves based on desired normalized RF distance.
	RF_max = 2*(tips-2)
	NNI_moves_cloud = int((RF_max * RF_norm_cloud)/2)
	if NNI_moves_cloud == 0:
		NNI_moves_cloud = 1
	# Make a cloud for each starting tree
	c_size = int(cloud_size)-1 

	print("Making clouds...")

	# Make clouds
	cluster1 = readTree.NNI_mult_trees(in_tree=t1,num_out_trees=c_size,num_nni_moves=NNI_moves_cloud,out='list')
	cluster2 = readTree.NNI_mult_trees(in_tree=t2,num_out_trees=c_size,num_nni_moves=NNI_moves_cloud,out='list')
	cluster3 = readTree.NNI_mult_trees(in_tree=t3,num_out_trees=c_size,num_nni_moves=NNI_moves_cloud,out='list')
	tree_list=[cluster1,cluster2,cluster3]
	# Make a nexus file with starting trees and cloud trees
	readTree.list_to_out(tree_list, '%s_cloud.tree' % name)


############################################################  
#User input
############################################################
def main():
	# User input
	tips = 75
	cloud_size = 1000
	RF_norm_cloud = 0.125
	RF_norm_start = 1.0
	number_replicates = 100
	# Make a bunch of trees
	for num in range(1,number_replicates+1):
		name =("%stip_%strees_%s_%sstart_%s" % (tips,cloud_size,RF_norm_cloud,RF_norm_start, num))
		trees = makeStarting(tips,RF_norm_start,name,RF_norm_cloud,cloud_size)
		makeCloud(tips,RF_norm_cloud,name,cloud_size,trees)
	cloud_size = 1000
	# Make a bunch of trees
	for num in range(1,number_replicates+1):
		name =("%stip_%strees_%s_%sstart_%s" % (tips,cloud_size,RF_norm_cloud,RF_norm_start, num))
		trees = makeStarting(tips,RF_norm_start,name,RF_norm_cloud,cloud_size)
		makeCloud(tips,RF_norm_cloud,name,cloud_size,trees)

	print("Making into nexus files...")
	# Turn script into Nexus that is readable by TreeScaper
	# This script turns all *.tree files into *.nex files 
	subprocess.call(['./makeNexus.sh'])
	print("Done")

#clear
# rm *.tree
# rm *.log

if __name__=='__main__':
	main()

