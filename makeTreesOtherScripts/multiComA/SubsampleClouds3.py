#! /usr/bin/env python
import os
import glob
import subprocess

# File name example
# 10tip_1000trees_0.125_0.5start_10_cloud.nex
# 10 lines header
# + num taxa
# + line
# tress start
# don't forget to copy end


def parseNumber(string,word):
	string=str(string)
	findWord = string.find(word)
	number = string[:findWord]
	return float(number)

def subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName):
	for t in glob.glob(suffix):
		# Split off name, dependant on what your iterating handle files are
		listName = t.split("_")
		# Replace name with 500, remove existing file, open new empty file
		newName = listName[0]+"_"+str(ratioName)+"trees_"+listName[2]+"_"+listName[3]+"_"+listName[4]+"_"+listName[5]
		#os.system("rm %s" % newName)
		outFile = open(newName,'a')
		# Get number of tips from file name
		numTips = parseNumber(listName[0],"tip")
		# Count lines to nexus header plus new size of first cloud
		firstCloud = 11 + int(numTips) + cloud1
		# Count lines to second cloud and size of new second cloud
		originalSize = parseNumber(listName[1],"trees")
		secondStart = 11 + int(numTips) + originalSize
		secondCloud = secondStart + cloud2
		# Count lines to second cloud and size of new second cloud
		thirdStart = secondStart + originalSize
		thirdCloud = thirdStart + cloud3
		# Print out info 
		print(newName)
		print("originalSize, secondStart, secondCloud, thirdStart, thirdCloud")
		print(originalSize, secondStart, secondCloud, thirdStart, thirdCloud)
		# Iterate through file 
		with open(t, 'r') as nexus:
			for i, line in enumerate(nexus):
				if i < firstCloud:
					#print(line)
					outFile.write(line)
				elif secondCloud > i >= secondStart:
					#print(line)
					outFile.write(line)
				elif thirdCloud > i >= thirdStart:
					#print(line)
					outFile.write(line)
			outFile.write("End;")
		outFile.close()

suffix="*cloud.nex"
os.system("mkdir subsampledFiles")

ratioName="1v5"
cloud1=167
cloud2=833
cloud3=833
subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName)
os.system("mv *%s* subsampledFiles" % ratioName)

ratioName="1v10"
cloud1=91
cloud2=909
cloud3=909
subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName)
os.system("mv *%s* subsampledFiles" % ratioName)

ratioName="1v20"
cloud1=48
cloud2=952
cloud3=952
subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName)
os.system("mv *%s* subsampledFiles" % ratioName)

ratioName="1v50"
cloud1=20
cloud2=980
cloud3=980
subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName)
os.system("mv *%s* subsampledFiles" % ratioName)

ratioName="1v100"
cloud1=10
cloud2=990
cloud3=990
subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName)
os.system("mv *%s* subsampledFiles" % ratioName)

ratioName="1v500"
cloud1=2
cloud2=998
cloud3=998
subSetTreeFile(suffix,cloud1,cloud2,cloud3,ratioName)
os.system("mv *%s* subsampledFiles" % ratioName)




