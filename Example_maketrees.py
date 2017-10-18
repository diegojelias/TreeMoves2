
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
	NNI_moves_cloud = int((RF_max * RF_norm_cloud)/2)
	# Create second starting tree
	t2 = readTree.NNI_mult_moves(in_tree=t1,num_moves=NNI_moves_start,node_choice='random',no_dup_start_tree='F', req_min_RF=RF_norm_start)

	# Write out tree files
	readTree.write_single_tree(t1,'%s_starting_tree_1.tree' % name)
	readTree.write_single_tree(t2,'%s_starting_tree_2.tree' % name)

	# Write out log file
	# Calculate emperical distance between two start trees as a gut check
	RF_emp = readTree.rf_unweighted(t1,t2,normalized='T')[1]
	# write to log file
	with open('%s.log' % name, "w") as log_file:
		line1 = "File name: "+str(name)
		line1b = "Tips: "+str(tips)+", Trees per cloud: "+str(cloud_size)
		line2 = "Starting trees - RF_input: "+str(RF_norm_start)+", RF_calc: "+str(RF_emp)+", NNI_moves: "+str(NNI_moves_start)
		line3 = "Cloud of trees - RF_input: "+str(RF_norm_cloud)+", NNI_moves: "+str(NNI_moves_cloud)
		log_file.write("%s\n%s\n%s" % (line1, line2, line3))

	# Pass tree files to next function
	return t1,t2
def makeCloud(tips,RF_norm_cloud,name,cloud_size,starting_trees): 
	# Make cluster of trees around each starting tree
	# Get starting trees parsed
	t1=starting_trees[0]
	t2=starting_trees[1]
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
	tree_list=[cluster1,cluster2]
	# Make a nexus file with starting trees and cloud trees
	readTree.list_to_out(tree_list,'%s_cloud.tree' % name)
def makeFolders(current_folder):
	os.chdir(current_folder)
	mainDir = os.getcwd()
	for t in glob.glob('*.nex'):

		# Get base name of file
		treeSet=str(t)
		treeSetIndex = treeSet.find(".nex")
		fName = treeSet[:treeSetIndex]

		# Create path 
		dirPath = os.path.join(mainDir,fName)
		print(dirPath)

		# Make directory 
		if not os.path.exists(dirPath):
			#print("NEW DIR")
			os.mkdir(dirPath)	

		# Copy nexus to new folder
		os.system("cp %s %s" % (t,dirPath))
		os.system("cp %s %s" % ("CLVTreeScaper",dirPath))
def runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec"):

	treeSet=str(inNexus)
	treeSetIndex = treeSet.find(".nex")
	treeSetTrunc = treeSet[:treeSetIndex]




	# The following is specialized for my naming scheme
	tips = treeSetTrunc.split("_")[0].strip("tip")
	trees = treeSetTrunc.split("_")[1].strip("trees")
	cloudRF = treeSetTrunc.split("_")[2]
	startRF = treeSetTrunc.split("_")[3].strip("start")
	
	outLog = treeSetTrunc+"_"+model

	if network == 'Covariance':

		#os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t %s -cm Covariance -lm auto -hf %s -lf %s > %s_CovAuto.out" % (clvPath, treeSet, weighted, rooted, model, hf, lf, treeSet))
		os.system("./CLVTreeScaper -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm auto -hf %s -lf %s > %s_CovAuto.out" % (treeSet, weighted, rooted, model, hf, lf, outLog))


		print("./CLVTreeScaper -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm auto -hf %s -lf %s > %s_CovAuto.out" % (treeSet, weighted, rooted, model, hf, lf, outLog))

	if network == 'Affinity':

		os.system("./CLVTreeScaper -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm auto -dm %s -am %s > %s_AffAuto.out" % (treeSet, weighted, rooted, model, dm, am, outLog))
		print(("./CLVTreeScaper -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm auto -dm %s -am %s > %s_AffAuto.out" % (treeSet, weighted, rooted, model, dm, am, outLog)))
def runCLV(mainDir):
	# iterates through folders and runs cov and Aff for all models
	for t in glob.glob('*.nex'):
		os.chdir(maindDir)
		inNexus = t
		treeSet=str(inNexus)
		treeSetIndex = treeSet.find(".nex")
		treeSetTrunc = treeSet[:treeSetIndex]

		
		dirPath = os.path.join(mainDir,treeSetTrunc)
		os.chdir(dirPath)
		print("Changing to folder: "+str(dirPath))

		network = 'Covariance'
		model = 'CNM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")
		model = 'CPM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")
		model = 'ERNM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")
		model = 'NNM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")

		network = 'Affinity'
		model = 'CNM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")
		model = 'CPM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")
		model = 'ERNM'
		runTreescaper(inNexus, network, model, weighted=0, rooted=1, lf=0.05, hf=0.95, dm="URF", am="Rec")

		os.chdir(mainDir)

############################################################  
#User input
############################################################
def main():
	# User input
	tips = 10
	cloud_size = 10
	RF_norm_cloud = 0.125
	RF_norm_start = 1.0
	number_replicates = 2
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

	#I think this is to also make folders and run everything. Haven't tested recently. 
	#mainDir = os.getcwd()
	#makeFolders(mainDir)
	#mainDir = os.getcwd()
	#runCLV(mainDir)

if __name__=='__main__':
	main()

