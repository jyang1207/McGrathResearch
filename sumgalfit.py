
#! /usr/bin/python


import os
import sys
import re
import math
import fnmatch
import time
import argparse
import redShiftConverter

def compute_distance(p0, p1):
	'''
	http://stackoverflow.com/questions/5407969/distance-formula-between-two-points-in-a-list
	'''
	
	return math.sqrt((float(p0[0]) - float(p1[0]))**2 + (float(p0[1]) - float(p1[1]))**2)

def sum_galfit(resultFilename):
	'''
	returns a string summary of the result specified by the parameter result filename
	
	parameter resultFilename - the result filename from running galfit to be summarized
	return - a string summarizing the results, ending in a new line
	'''
	
	resultFile = open(resultFilename, 'r')
	
	resultLines = resultFile.readlines()
	
	resultFile.close()
	
	skipSky = False
	galaxyID = ""
	px1 = ""
	px2 = ""
	mag1 = ""
	mag2 = ""
	rad1 = ""
	rad2 = ""
	sersIndex1 = ""
	sersIndex2 = ""
	ba1 = ""
	ba2 = ""
	pa1 = ""
	pa2 = ""
	
	for resultLine in resultLines:
		
		if not galaxyID and resultLine.strip()[:2] == "B)":
		
			# B) a0.220/VELA01_220_cam0_F160W_multi.fits      # Output data image block
			galaxyID = resultLine.split("_")[0].split("/")[1]
			timeStep = resultLine.split("_")[1]
			camera = resultLine.split("_")[2]
			filter = resultLine.split("_")[3]
	
					
		if resultLine.strip()[:6] == "0) sky":
			skipSky = True
		elif resultLine.strip()[:6] == "0) ser":
			skipSky = False
			
	
		if not skipSky:
			if not px1 and resultLine.strip()[:2] == "1)":
			
				#1) 301.6210 299.7872 1 1  #  Position x, y
				px1 = resultLine.strip().split(" ")[1]
				py1 = resultLine.strip().split(" ")[2]
	
			elif not px2 and resultLine.strip()[:2] == "1)":
				px2 = resultLine.strip().split(" ")[1]
				py2 = resultLine.strip().split(" ")[2]
	
			elif not mag1 and resultLine.strip()[:2] == "3)":
				mag1 = resultLine.strip().split(" ")[1]
				
			elif not mag2 and resultLine.strip()[:2] == "3)":
				mag2 = resultLine.strip().split(" ")[1]
	
			elif not rad1 and resultLine.strip()[:2] == "4)":
				rad1 = resultLine.strip().split(" ")[1]
				
			elif not rad2 and resultLine.strip()[:2] == "4)":
				rad2 = resultLine.strip().split(" ")[1]
				
			elif not sersIndex1 and resultLine.strip()[:2] == "5)":
				sersIndex1 = resultLine.strip().split(" ")[1]
				
			elif not sersIndex2 and resultLine.strip()[:2] == "5)":
				sersIndex2 = resultLine.strip().split(" ")[1]
				
			elif not ba1 and resultLine.strip()[:2] == "9)":
				ba1 = resultLine.strip().split(" ")[1]
				
			elif not ba2 and resultLine.strip()[:2] == "9)":
				ba2 = resultLine.strip().split(" ")[1]
				
			elif not pa1 and resultLine.strip()[:3] == "10)":
				pa1 = resultLine.strip().split(" ")[1]
				
			elif not pa2 and resultLine.strip()[:3] == "10)":
				pa2 = resultLine.strip().split(" ")[1]
	
	# a = 1/(1+z)
	# z = 1/a - 1
	timeStep = str(redShiftConverter.red_shift_to_gyr(1/(float(timeStep)/1000) - 1))
	
	
	errorFlag1 = ""
	errorFlag2 = ""
	
	# sersic index of 0.5 < si < 10 is good
	if float(sersIndex1) < 0.5 or float(sersIndex1) > 10:
		sersFlag1 = "*"
		errorFlag1 = "*"
	else:
		sersFlag1 = ""

	
	if not px2:
		if float(sersIndex1) > 2.5:
			type = "bulge"
		else:
			type = "disk"
			
		return (errorFlag1 + galaxyID + ", " + timeStep + ", " + camera + ", " + filter + ", " +
				type + ", " + px1 + ", " + py1 + ", " + sersIndex1 + ", " + mag1 + ", " + 
				rad1 + ", " + ba1 + ", " + pa1 + "\n")
				
	
	# component seperation of greater than 5 is bad 
	dist = compute_distance([px1, py1], [px2, py2])
	print dist
	if dist > 5:
		posFlag = "*"
		errorFlag1 = "*"
		errorFlag2 = "*"
	else:
		posFlag = ""
		
		
	# sersic index of 0.5 < si < 10 is good
	if float(sersIndex2) < 0.5 or float(sersIndex2) > 10:
		sersFlag2 = "*"
		errorFlag2 = "*"
	else:
		sersFlag2 = ""
		
	# test for type (bulge or disk)
	if float(sersIndex1) > 2.5 and float(sersIndex2) < 2.5:
		type1 = "bulge"
		type2 = "disk"
	elif float(sersIndex1) < 2.5 and float(sersIndex2) > 2.5:
		type2 = "bulge"
		type1 = "disk"
	else:
		if rad1 < rad2:
			type1 = "bulge"
			type2 = "disk"
		else:
			type2 = "bulge"
			type1 = "disk"
			
	result1 = (errorFlag1 + galaxyID + ", " + timeStep + ", " + camera + ", " + filter + ", " +
				type1 + ", " + posFlag + px1 + ", " + py1 + ", " + sersFlag1 + sersIndex1 + ", " + 
				mag1 + ", " + rad1 + ", " + ba1 + ", " + pa1 + "\n")
	result2 = (errorFlag2 + galaxyID + ", " + timeStep + ", " + camera + ", " + filter + ", " +
				type1 + ", " + posFlag + px2 + ", " + py2 + ", " + sersFlag2 + sersIndex2 + ", " + 
				mag2 + ", " + rad2 + ", " + ba2 + ", " + pa2 + "\n")
	return result1 + result2

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
		
	#outFile = open(outFileStr, 'w')
	
	#outFile.write("galfit result file\n" +
	#			"galaxy ID, time step, camera, filter, type, px, py, sersic, mag, rad, b/a, angle\n")
	
	#outFile.close()
	
	outFile = open(outFileStr, 'a')
	
	# this loops through every image in images file and
	for resultFilename in resultFilenames:
	
		# summarize galfit and write to output
		outFile.write( sum_galfit(resultFilename.strip()) )
		
	outFile.close()

######################### done ################################################