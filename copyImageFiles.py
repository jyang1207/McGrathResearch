
#! /usr/bin/python


import os
import argparse

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

	# used to parse command line arguments
	parser = argparse.ArgumentParser()
						
	parser.add_argument("input", 
						help="set the file containing the list of full path image filenames",
						type=parseFile)
	
	parser.add_argument("destination", 
						help="set the top level directory to where images will be copied",
						type=parseDirectory)
	
	# parse the command line using above parameter rules
	args = parser.parse_args()
	
	# read all the lines of the input file to get a list of full path image filenames
	inputFile = open(args.input, 'r')
	imageFilenames = inputFile.readlines()
	inputFile.close()
	
	# this loops through every image in images file and removes new line
	imageFilenames = [ imageFilename.strip() for imageFilename in imageFilenames ]
	
	# copy each image into the destination given, preserving folder structure
	for imageFilename in imageFilenames:

		# destDirectory is originally the top level given on command line
		destDirectory = args.destination
		
		# for images with surrounding folder structure, duplicate folder structure in destination
		if len(imageFilename.split("/")) != 1:
			for directory in imageFilename.split("/")[:-1]:
				destDirectory = destDirectory + directory + "/"
				if not os.path.isdir(destDirectory):
					os.system("mkdir " + destDirectory)
		
		# destDirectory contains all folder structure information, if any
		os.system("cp " + imageFilename + " " + destDirectory + imageFilename.split("/")[-1])

		
		
		
		
		