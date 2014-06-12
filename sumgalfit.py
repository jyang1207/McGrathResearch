
#! /usr/bin/python


import os
import sys
import re
import fnmatch
import time
import argparse


def sum_galfit(resultFilename):
	'''
	returns a string summary of the result specified by the parameter result filename
	
	parameter resultFilename - the result filename from running galfit to be summarized
	return - a string summarizing the results, ending in a new line
	'''
	return resultFilename + "\n"

def parseDirectory(d):
	'''	
	raises an argument exception if the string d is not a directory
	modifies d to ensure that the directory ends with a forward slash
	
	parameter d - the string to be checked as a directory
	returns - parameter d with an appended forward slash
	'''
	if not os.path.isdir(d):
		msg = "directory {} either does not exist or in not accessible".format(d)
		raise argparse.ArgumentTypeError(msg)
	elif d[-1] != "/":
		d = d + "/"
		
	return d
	
	
def parseFile(f):
	'''	
	raises an argument exception if the string f is not a file
	
	parameter f - the string to be checked as a file
	returns - parameter f, unchanged
	'''
	if not os.path.isfile(f):
		msg = "file {} either does not exist or in not accessible".format(f)
		raise argparse.ArgumentTypeError(msg)
		
	return f


if __name__ == "__main__":

	singlePattern = "*_single_result.txt"	
	bulgePattern = "*_bulge_result.txt"	

	# used to parse command line arguments
	parser = argparse.ArgumentParser()
	
	# directory and file are mutually exclusive parameters
	group = parser.add_mutually_exclusive_group()
	
	# directory specifies the directory where the images are
	group.add_argument("-d","--directory", 
						help="set the directory containing the galfit results to summarize",
						type=parseDirectory, default="./")
	
	# file specifies the full path filename of the list of images to run
	group.add_argument("-f","--file", 
						help="set the file containing the list of full galfit result filenames",
						type=parseFile)
						
	# bulge is a boolean (true or false) specifying if the simulation should fit
	# an additional component after the initial fit from imexam results
	parser.add_argument("-b","--bulge", 
						help="turn on to summarize galfit bulge fit results (off for single)",
						action="store_true")
						
	# file specifies the full path filename of the list of images to run
	parser.add_argument("-o","--output", 
						help="set the output file where summary will be stored",
						type=parseFile)
						
	# Magnitude photometric zeropoint					
	# Plate scale
	# PSF
	
	# parse the command line using above parameter rules
	args = parser.parse_args()
	
	# set the list of results by the command line argument
	if args.file:
		resultListFilename = args.file
		
	# set the list of results by filtering in the specified directory (default current)
	elif args.bulge:
		resultListFilename = args.directory + "all_cam_bulge_results_" + time.strftime("%m-%d-%Y") + ".txt"
		os.system("ls " + args.directory + bulgePattern + " > " + resultListFilename)
	else:
		resultListFilename = args.directory + "all_cam_single_results_" + time.strftime("%m-%d-%Y") + ".txt"
		os.system("ls " + args.directory + singlePattern + " > " + resultListFilename)
		
	#this will be the file that will contain the images
	r = open(resultListFilename, 'r')		
	
	# the first line of file containing images
	resultFilenames = r.readlines()
	
	# close the file now that it has been read
	r.close()
	
	# the output file
	if not args.output:
		outFileStr = "galfit_result_summary.txt"
	else:
		outFileStr = args.output
		
	outFile = open(outFileStr, 'w')
	
	outFile.write("galfit result file\n")
	
	outFile.close()
	
	outFile = open(outFileStr, 'a')
	
	# this loops through every image in images file and
	for resultFilename in resultFilenames:
	
		# summarize galfit and write to output
		outFile.write( sum_galfit(resultFilename.strip()) )
		
	outFile.close()

######################### done ################################################