
#! /usr/bin/python


import os
import sys
import re
from pyraf import iraf
import math
import time
import argparse
import fnmatch


def run_imhead(imageFilename):
	'''
	invokes the iraf method imhead on the given image filename
	parses the result to be returned as a list of string values
	
	parameter imageFilename - the full path filename of the image on which to invoke imexam
	
	returns - 
		a string list of image info in the form 
		[galaxy_id, filter, cam_number, image_number, frame_height, frame_width]
	'''
	
	imhead_return = iraf.imhead(imageFilename, Stdout=1)[0].strip()

	# "VELA01_a0.110_0006317__skipir_CAMERA0-BROADBAND_F125W_simulation.fits[600,600][real]:" 
	
	
	frame_dimensions = imhead_return.split("[")[1].replace("]", "").split(",")
	# ["600","600"]

	frame_width = frame_dimensions[0]
	# "600"
	
	frame_height = frame_dimensions[1]
	# "600"
	
	galaxy_id = imhead_return.split("_")[0]
	#"VELA01"
	
	filter = imhead_return.split("_")[6]
	#"F125W"	
	
	image_number = imhead_return.split("_")[1].split(".")[1]
	# "110"


	
	cam_str = imhead_return.split("_")[5].split("-")[0]
	# "CAMERA0"

	
	#these statements accomodate for camera numbers greater than 9
	conditions = ["1","2","3","4","5","6","7","8","9"]			
	
	#determines if camera number is greater than 10, does not allow greater than 99
	# ex: "CAMERA0" -> cam_str[-2] = "M"
	# ex: "CAMERA12" -> cam_str[-2] = "1"
	if cam_str[-2] in conditions:										
		cam_number = cam_str[-2:]
	else:
		cam_number = cam_str[-1]					
	#print cam_number
	
	return [galaxy_id, filter, cam_number, image_number, frame_height, frame_width]
	
				
def run_imexam(image):
	'''
	invokes the iraf method imexam on the given image filename
	parses the result to be returned as a list of string values
	
	parameter image - the full path filename of the image on which to invoke imexam
	
	returns - list of string image info in the form [line, col, mag, radius, b_a, PA]
	'''
	
	#writes the output of the minmax function of iraf to a file for later use. 
	os.system('touch'+' coords.tmp')
	write_coords = open('coords.tmp', 'w')
	
	#call iraf's minmax function
	maxCoords = iraf.minmax(image, Stdout=1, update=0, force=1
		)[0].strip().split(" ")[3].replace("[", "").replace("]", "").split(",")
		
	#print maxCoords
	#['VELA01_a0.110_0006317__skipir_CAMERA0-BROADBAND_F125W_simulation.fits[1:600,1:300]', '[1,1]', 
	#'0.', '[363,285]', '0.1987768709659576']
	write_coords.write(maxCoords[0]) 				
	write_coords.write(" " + maxCoords[1])
	write_coords.close()	
	
	#run imexam passing coords.tmp
	imexam_out = str(iraf.imexam(image, use_display=0, imagecur="coords.tmp", Stdout=1))

	# "display frame (1:) (1): ['#   COL    LINE    COORDINATES', '#     " +
	# "R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT', " +
	# "'1283.40  692.16 1283.40 692.16', '  63.22  14.53  15461.  0.1202   19.53 0.52    " +
	# "1 1.75    21.58    18.85  21.39']"	
	
	# delete coords.tmp, no longer needed
	os.system('rm' + ' coords.tmp')	
	
	imexam_array = str.split(imexam_out, "'")
	# ["display frame (1:) (1): [", 
	# "#   COL    LINE    COORDINATES", 
	# ", ", "#     " +
	# "R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT", ", ",
	# "1283.40  692.16 1283.40 692.16", ", ", "  63.22  14.53  15461.  0.1202   19.53 0.52    " +
	# "1 1.75    21.58    18.85  21.39']"			
	
	#print imexam_out
	#print imexam_array	
	
	xy = imexam_array[5].strip()
	data = imexam_array[7].strip()
	data = str.split(data)
	xy = str.split(xy)
	
	col = xy[0] 									#this gives x value	
	line = xy[1]									#this gives y value
	mag = data[1]									#this prints the magnitude
	radius = data[0]								#this prints the radius
	b_a = int(1) - float(data[5])					#this prints the E variable that relates to b/a. note b/a=1-e
	PA = float(data[6]) - int(90)					#this gives the position angle. iraf measures up from x. Galfit down from y

# 	print col
# 	print line
# 	print mag
# 	print radius
# 	print b_a
# 	print PA
# 	exit()	
	
	return [line, col, mag, radius, b_a, PA]
	

def run_galfit(imageFilename, galfit_constraint_filename, includeBulgeComponent):
	'''
	opens file (parameter) containing list of images and loops over every image, 
	running galfit and writing the results to the same directory as the images
	
	calls methods to invoke iraf methods to write the initial galfit param file
	for each image
	
	parameter imageListFilename -
		the string of the full path filename of the text file holding the list
		of images on which galfit will be run
		
	parameter includeBulgeComponent - 
		boolean indicating if a bulge should be fit after first pass by galfit
	'''
	
	# run the method of dimensions.py
	[galaxy_id, filter, cam_number, image_number, height, width] = run_imhead(imageFilename)
	
	# calls a method of the imexam.py file, passing filename (with width and height) as a parameter
	# imexam gets the coordinates of the max pixel and uses iraf.imexam to set global
	# variables detailing the results of iraf's examination of the max coordinate
	[Y, X, magnitude, rad, BA, angle] = run_imexam(
		imageFilename + '[' + str(int(width)/2 - 75) + ':' + str(int(width)/2 + 75) + "," +
							  str(int(height)/2 - 75) + ':' + str(int(height)/2 + 75) + ']')
	
	Y = str(float(Y) + float(height)/2 - 75.0)
	X = str(float(X) + float(width)/2 - 75.0)
	
	# define filenames, galaxy_id has full path if any
	filename = galaxy_id + "_" + image_number + '_cam' + str(cam_number) + '_' + filter
	
	#galfit_output_filename = filename + "_multi.fits"
	galfit_single_parameter_filename = 	filename + '_single_param.txt'
	galfit_single_output_filename = 	filename + "_single_multi.fits"
	galfit_single_result_filename = 	filename + "_single_result.txt"
	galfit_bulge_parameter_filename = 	filename + '_bulge_param.txt'
	galfit_bulge_output_filename = 		filename + "_bulge_multi.fits"
	galfit_bulge_result_filename = 		filename + "_bulge_result.txt"
	
	os.system('touch ' + galfit_single_parameter_filename)								

######### writes gathered galfit parameters to file ###########################
	WP = open(galfit_single_parameter_filename,'w')
	WP.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
	WP.write("A) " + imageFilename + 
			"						#Input data image block\n")
	WP.write("B) " + galfit_single_output_filename +
			"						#Output data image block\n")
	WP.write("C)" + " none" + 
			"						#Sigma image name (made from data if blank or 'none')\n")
	WP.write("D)" + " none" + 
			"						#Input PSF image and (optional) diffusion kernel\n")#command line TODO
	WP.write("E)" + " 1" + 
			"						#PSF fine sampling factor relative to data\n")
	WP.write("F)" + " none" + 							
			"						#Bad pixel mask (FITS file or ASCIIcoord list)\n")
	WP.write("G)" + " none" + 
			"						#File with parameter constraints (ASCII file)\n")
	WP.write("H)" + " 1	" + width + " 1	" + height + 
			"						#Image region to fit (xmin xmax ymin ymax)\n")
	WP.write("I)" + " 200 200" + 
			"						#Size of the concolution box (x y)\n")
	WP.write("J)" + " 26.3" + 
			"						#Magnitude photometric zeropoint\n")#command line TODO 
	WP.write("K)" + " 0.06" + "  0.06" + 
			"						#Plate scale (dx dy)  [arcsec per pixel]\n")#command line TODO 
	WP.write("O)" + " regular" + 
			"						#display type (regular, curses, both\n")
	WP.write("P)" + " 0" + 
			"						#Options: 0=normal run; 1,2=make model/imgblock & quit\n")
	WP.write("S)" + " 0" + 
			"						#Modify/create components interactively?\n")
	WP.write("\n")
	WP.write("# INITIAL FITTING PARAMETERS\n")
	WP.write("#\n")
	WP.write("#	For component type, the allowed functions are:\n")
	WP.write("#		nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n")
	WP.write("#		ferrer, coresersic, sky and isophote.\n")
	WP.write("#\n")
	WP.write("#	Hidden parameters will only appear when they are specified:\n")
	WP.write("#		C0 (diskyness/boxyness),\n")
	WP.write("#		Fn (n=interger, Azimuthal Fourier Modes).\n")
	WP.write("#		R0-R10 (PA rotation, for creating spiral structures).\n")
	WP.write("#\n")
	WP.write("# ------------------------------------------------------------------------------\n")
	WP.write("#		par)	par value(s)	fit toggle(s)	# parameter description\n")
	WP.write("# ------------------------------------------------------------------------------\n")
	WP.write("\n")
	WP.write("# Componenet number: 1\n")
	WP.write(" 0) sersic					#Component type\n")
	WP.write(" 1) " + str(X) + "	" + str(Y) + "	1	1			#Position x,y\n")
	WP.write(" 3) " + str(magnitude) + "	1			#Integrated Magnitude\n")
	WP.write(" 4) " + str(rad) + "			1			#R_e (half-light radius)	[pix]\n")
	WP.write(" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
	WP.write(" 6) 0.0000		0			#	-----\n")
	WP.write(" 7) 0.0000		0			#	-----\n")
	WP.write(" 8) 0.0000		0			#	-----\n")
	WP.write(" 9) " + str(BA) + "			1			#Axis ratio (b/a)\n")
	WP.write(" 10) " + str(angle) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
	WP.write(" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
	WP.write("\n")


	WP.write("# Componenet number: 2\n")
	WP.write(" 0) sky						#Component type\n")
	WP.write(" 1) 0.0000		0			#	Sky background at center of fitting region [ADUs]\n")
	WP.write(" 2) 0.0000		0			#	dsky/dx (sky gradient in x) [ADUs/pix]\n")
	WP.write(" 3) 0.0000		0			#	dsky/dy (sky gradient in y) [ADUs/pix]\n")
	WP.write(" Z) 0						#Leave in [1] or subtract [0] this comp from data?\n")
	WP.write("\n")

	WP.close()


	# run galfit on paramter file
	os.system('galfit ' + galfit_single_parameter_filename)

	# remove temp files
	os.system("rm " + "fit.log")
	os.system("mv galfit.01 " + galfit_single_result_filename)
	#os.system("mv " + galfit_output_filename + " " + galfit_single_output_filename)
	
	if includeBulgeComponent:
		# here we would write the third component to the end of galfit_result_filename
		singleResultFile = open(galfit_single_result_filename, "r")		
		resultContents = singleResultFile.readlines()
		singleResultFile.close()
		
		# gather info about first fit while modifying output and coonstraint parameters
		
		bulgeParam = ""		
		positionLine = ""
		magnitudeLine = ""
		for resultLine in resultContents:
		
			if resultLine.strip()[:2] == "B)":
				bulgeParam = (bulgeParam + "B) " + galfit_bulge_output_filename +
							"						#Output data image block\n")
				
			elif resultLine.strip()[:2] == "G)":
				bulgeParam = (bulgeParam + "G) " + galfit_constraint_filename + 
							"						#File with parameter constraints (ASCII file)\n")
			
			elif not positionLine and resultLine.strip()[:2] == "1)":
				positionLine = resultLine.strip().split(" ")
				bulgeParam = bulgeParam + resultLine
				
			elif not magnitudeLine and resultLine.strip()[:2] == "3)":
				magnitudeLine = resultLine.strip().split(" ")
				bulgeParam = bulgeParam + resultLine
				
			else:
				bulgeParam = bulgeParam + resultLine
				
		
		resultX = positionLine[1]
		resultY = positionLine[2]
		
		resultMag = magnitudeLine[1]
		
		
		#append bulge component using above info
		
		bulgeParam = bulgeParam + ("# Componenet number: 3\n")
		bulgeParam = bulgeParam + (" 0) sersic					#Component type\n")
		bulgeParam = bulgeParam + (" 1) " + resultX + "	" + resultY + "	1	1			#Position x,y\n")
		bulgeParam = bulgeParam + (" 3) " + resultMag + "	1			#Integrated Magnitude\n")
		bulgeParam = bulgeParam + (" 4) " + str(rad) + "			1			#R_e (half-light radius)	[pix]\n")
		bulgeParam = bulgeParam + (" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
		bulgeParam = bulgeParam + (" 6) 0.0000		0			#	-----\n")
		bulgeParam = bulgeParam + (" 7) 0.0000		0			#	-----\n")
		bulgeParam = bulgeParam + (" 8) 0.0000		0			#	-----\n")
		bulgeParam = bulgeParam + (" 9) 1" + "			1			#Axis ratio (b/a)\n")
		bulgeParam = bulgeParam + (" 10) 0" + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
		bulgeParam = bulgeParam + (" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
		bulgeParam = bulgeParam + ("\n")


		bulgeParamFile = open(galfit_bulge_parameter_filename, "w")
		bulgeParamFile.write(bulgeParam)
		bulgeParamFile.close()
		
		# run galfit on paramter file
		os.system('galfit ' + galfit_bulge_parameter_filename)
	
		# remove temp files
		os.system("rm " + "fit.log")
		os.system("mv galfit.01 " + galfit_bulge_result_filename)
		#os.system('mv ' + galfit_output_filename + ' ' + galfit_bulge_output_filename)


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

	# directory and file are mutually exclusive parameters
	exclusive_group = parser.add_mutually_exclusive_group()
	
	# directory specifies the directory where the images are
	parser.add_argument("-d","--directory", 
						help="set the directory containing the images on which to run galfit",
						type=parseDirectory, default="./")
	
	# file specifies the full path filename of the list of images to run
	exclusive_group.add_argument("-f","--file", 
						help="set the file pattern to match image filenames in given directory",
						default="*_simulation.fits")
						
	exclusive_group.add_argument("-i","--input", 
						help="set the file containing the list of full path image filenames",
						type=parseFile)
						
	# bulge is a boolean (true or false) specifying if the simulation should fit
	# an additional component after the initial fit from imexam results
	parser.add_argument("-b","--bulge", 
						help="turn on to include a bulge fit after the initial galaxy fit",
						action="store_true")
						
	parser.add_argument("-c","--constraint", 
						help="set the file containing the galfit constraints",
						type=parseFile)
						
	# Magnitude photometric zeropoint					
	# Plate scale
	# PSF
	
	# parse the command line using above parameter rules
	args = parser.parse_args()
	if not args.input:
		imageFilenames = fnmatch.filter(os.listdir(args.directory), args.file)
		if len(imageFilenames) == 0:
			print "file pattern {} in directory {} does not match any files".format(
								args.file, args.directory)
			parser.print_help()
		
		# this loops through every image in images file and remove \n
		imageFilenames = [args.directory + imageFilename.strip() for imageFilename in imageFilenames ]
	else:
		inputFile = open(args.directory + args.input, 'r')
		imageFilenames = inputFile.readlines()
		inputFile.close()
		
		# this loops through every image in images file and prefixes with directory
		imageFilenames = [imageFilename.strip() for imageFilename in imageFilenames ]
		

	if not args.constraint:
		constraintFilename = "none"
	else:
		constraintFilename = args.constraint
	
	# this loops through every image in images file and
	for imageFilename in imageFilenames:
	
		# run galfit
		run_galfit(imageFilename, constraintFilename, args.bulge)
	
######################### done ################################################
