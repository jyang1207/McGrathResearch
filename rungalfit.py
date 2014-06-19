
#! /usr/bin/python


import os
import sys
import re
from pyraf import iraf
import time
import fnmatch
from optparse import OptionParser


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
	imhead_info = imhead_return.split("/")[-1]
	directory_location = imhead_return.replace(imhead_info, "")
	# "VELA01_a0.110_0006317__skipir_CAMERA0-BROADBAND_F125W_simulation.fits[600,600][real]:" 
	
	
	frame_dimensions = imhead_info.split("[")[1].replace("]", "").split(",")
	# ["600","600"]

	frame_width = frame_dimensions[0]
	# "600"
	
	frame_height = frame_dimensions[1]
	# "600"
	
	galaxy_id = imhead_info.split("_")[0]
	#"VELA01"
	
	filter = imhead_info.split("_")[6]
	#"F125W"	
	
	image_number = imhead_info.split("_")[1].split(".")[1]
	# "110"
	
	cam_str = imhead_info.split("_")[5].split("-")[0]
	# "CAMERA0"

	
	#these statements accomodate for camera numbers greater than 9
	digits = ["1","2","3","4","5","6","7","8","9"]			
	
	#determines if camera number is greater than 10, does not allow greater than 99
	# ex: "CAMERA0" -> cam_str[-2] = "M"
	# ex: "CAMERA12" -> cam_str[-2] = "1"
	if cam_str[-2] in digits:										
		cam_number = cam_str[-2:]
	else:
		cam_number = cam_str[-1]					
	#print cam_number
	
	return [directory_location, galaxy_id, filter, cam_number, image_number, frame_height, frame_width]
	
	
def run_gauss(imageFilename, sigma):
	'''
	using the image given, calls iraf's gauss method and returns the
	filename of the result, which is blurred using the sigma given
	
	parameter imageFilename - the original image's filename
	parameter sigma - the number of pixels to blur the image by
	
	returns - the blurred image's filename
	'''
	newImageFilename = imageFilename[:-5] + "_gauss.fits"
	iraf.gauss(imageFilename, newImageFilename, sigma)
	return newImageFilename
	
	
def run_minmax(imageFilename, xStart, xStop, yStart, yStop):
	'''
	method calls iraf's minmax method on the area defined by the four
	parameter values that define two points on the image
	start is the top left, stop is the bottom right
	
	parameter imageFilename - the full path filename of the image on which to invoke imexam
	
	returns - list of two strings defining the coordinates of the max pixel ['xMax', 'yMax']
	'''
	
	# string representing the part of the image to run minmax and imexam on
	areaStr = ('[' + str(int(xStart)) + ':' + str(int(xStop)) + "," +
					 str(int(yStart)) + ':' + str(int(yStop)) + ']')
					 
	#call iraf's minmax function
	return iraf.minmax(imageFilename + areaStr, Stdout=1, update=0, force=1
		)[0].strip().split(" ")[3].replace("[", "").replace("]", "").split(",")
		
		
def run_imexam(imageFilename, centerCoords, xStart, xStop, yStart, yStop):
	'''
	method calls iraf's imexam method on the area defined by the four
	
	parameter imageFilename - the full path filename of the image on which to invoke imexam
	parameter centerCoords - the coordinates of the initial guess of the center of the galaxy
	parameter xStart, xStop, yStart, yStop - 
		two points on the original image defining an area to run minmax on
		start is the top left, stop is the bottom right
	
	returns - list of string image info in the form [line, col, mag, radius, b_a, PA]
	'''
		
	
	# string representing the part of the image to run minmax and imexam on
	areaStr = ('[' + str(int(xStart)) + ':' + str(int(xStop)) + "," +
					 str(int(yStart)) + ':' + str(int(yStop)) + ']')
					 
	#writes the output of the minmax function of iraf to a file for later use. 
	coordsFilename = "coords"
	
	# this is t prevent parallel processes from overwriting the coords.tmp file
	while os.path.isfile(coordsFilename + ".tmp"):
		coordsFilename = coordsFilename + "x"
	coordsFilename = coordsFilename + ".tmp"
	
	# create coords[x*].tmp for writing max coordinate to pass to imexam
	os.system('touch '+ coordsFilename)
	write_coords = open(coordsFilename, 'w')
	write_coords.write(centerCoords[0] + " " + centerCoords[1])
	write_coords.close()	
	
	#run imexam passing coords.tmp
	imexam_out = str(iraf.imexam(imageFilename + areaStr, use_display=0, imagecur=coordsFilename, Stdout=1))

	# "display frame (1:) (1): ['#	 COL	LINE	COORDINATES', '#	 " +
	# "R	MAG	   FLUX		SKY	   PEAK	   E   PA BETA ENCLOSED	  MOFFAT DIRECT', " +
	# "'1283.40	 692.16 1283.40 692.16', '	63.22  14.53  15461.  0.1202   19.53 0.52	 " +
	# "1 1.75	 21.58	  18.85	 21.39']"	
	
	# delete coords.tmp, no longer needed
	os.system('rm ' + coordsFilename) 
	
	imexam_array = str.split(imexam_out, "'")		
	
	#print imexam_out
	#print imexam_array 
	
	xy = imexam_array[5].strip()
	data = imexam_array[7].strip()
	data = str.split(data)
	xy = str.split(xy)
	
	col = xy[0]										#this gives x value 
	line = xy[1]									#this gives y value
	mag = data[1]									#this prints the magnitude
	radius = data[0]								#this prints the radius
	b_a = int(1) - float(data[5])					#this prints the E variable that relates to b/a. note b/a=1-e
	PA = float(data[6]) - int(90)					#this gives the position angle. iraf measures up from x. Galfit down from y

#	print col
#	print line
#	print mag
#	print radius
#	print b_a
#	print PA
#	exit()	
	
	return [line, col, mag, radius, b_a, PA]
	

def run_galfit(imageFilename, galfit_constraint_filename, includeBulgeComponent, includeGaussianSmoothing):
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
	returnStr = ""
	
	# sigma [# pixels] is used for run gauss when blurring
	sigma = 15
	
	# run the method of dimensions.py
	[directory_location, galaxy_id, filter, cam_number, image_number, height, width
		] = run_imhead(imageFilename)

	# the x y location of the top left corner of the area to run minmax and imexam on
	# in the coordinate system of the original image (which would be 0, 0 for the original)
	xStart = float(width)/2 - 75.0
	xStop = float(width)/2 + 75.0
	yStart = float(height)/2 - 75.0
	yStop = float(height)/2 + 75.0
	
	# calls iraf's minmax method, passing image filename as a parameter
	# with image possible gaussian smoothed, as well as two points
	# on the image defining the area on which to run minmax 
	# returns as a list of two strings the coordinates of the max pixel
	if includeGaussianSmoothing:
		centerCoords = run_minmax(run_gauss(imageFilename, sigma), 
									xStart, xStop, yStart, yStop)
	else:
		centerCoords = run_minmax(imageFilename, 
									xStart, xStop, yStart, yStop)
	
	# calls iraf's imexam method, passing filename as a parameter
	# along with center coordinate and two points defining area
	# returns initial estimates of the returned paramters
	[Y, X, magnitude, rad, BA, angle
		] = run_imexam(imageFilename, centerCoords, xStart, xStop, yStart, yStop)
	
	Y = str(float(Y) + yStart)
	X = str(float(X) + xStart)
	
	# define filenames, galaxy_id has full path if any
	filename = (directory_location + galaxy_id + "_" + 
				image_number + '_cam' + str(cam_number) + '_' + filter)
	
	#galfit_output_filename = filename + "_multi.fits"
	galfit_single_parameter_filename =	filename + '_single_param.txt'
	galfit_single_output_filename =		filename + "_single_multi.fits"
	galfit_single_result_filename =		filename + "_single_result.txt"
	galfit_bulge_parameter_filename =	filename + '_bulge_param.txt'
	galfit_bulge_output_filename =		filename + "_bulge_multi.fits"
	galfit_bulge_result_filename =		filename + "_bulge_result.txt"
	
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
	WP.write("H)" + " 1 " + width + " 1 " + height + 
			"						#Image region to fit (xmin xmax ymin ymax)\n")
	WP.write("I)" + " 200 200" + 
			"						#Size of the concolution box (x y)\n")
	WP.write("J)" + " 26.3" + 
			"						#Magnitude photometric zeropoint\n")#command line TODO 
	WP.write("K)" + " 0.06" + "	 0.06" + 
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
	WP.write("# For component type, the allowed functions are:\n")
	WP.write("#		nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n")
	WP.write("#		ferrer, coresersic, sky and isophote.\n")
	WP.write("#\n")
	WP.write("# Hidden parameters will only appear when they are specified:\n")
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

	if os.path.isfile(galfit_single_output_filename):
		# remove temp files
		os.system("rm " + "fit.log")
		os.system("mv galfit.01 " + galfit_single_result_filename)
		#os.system("mv " + galfit_output_filename + " " + galfit_single_output_filename)
	else:
		return "galfit failed on single component, probably mushroom"
	
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
		bulgeParam = bulgeParam + (" 1) " + resultX + " " + resultY + " 1	1			#Position x,y\n")
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
	
	
		if os.path.isfile(galfit_bulge_output_filename):
			# remove temp files
			os.system("rm " + "fit.log")
			os.system("mv galfit.01 " + galfit_bulge_result_filename)
			#os.system('mv ' + galfit_output_filename + ' ' + galfit_bulge_output_filename)
		else:
			return "galfit failed on bulge component, probably mushroom"

	return "success"
	
	
if __name__ == "__main__":

	usage = "usage: %prog [-h for help] [options] input" 

	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# bulge is a boolean (true or false) specifying if the simulation should
	# fit an additional component after the initial fit from imexam results
	parser.add_option("-b","--bulge", 
						help="turn on to include a bulge fit after the initial galaxy fit",
						action="store_true")
						
	# gauss is a boolean (true or false) specifying if the simulation should
	# apply gaussian smoothing before running minmax
	parser.add_option("-g","--gauss", 
						help="turn on to include a gaussian smoothing before minmax is run",
						action="store_true")
						
	parser.add_option("-c","--constraint", 
						help="set the file containing the galfit constraints")
						
	# Magnitude photometric zeropoint					
	# Plate scale
	# PSF
	
	# parse the command line using above parameter rules
	[options, args] = parser.parse_args()
	
	if len(args) != 1:
		parser.error("incorrect number of arguments, " + 
					"must provide an input file containing the list of full path image filenames.")
	elif not os.path.isfile(args[0]):
		parser.error("input file " + args[0] + " either does not exist or is not accessible")
	
	
	# read list of image filenames from input file
	inputFile = open(args[0], 'r')
	imageFilenames = inputFile.readlines()
	inputFile.close()
	
	# this loops through every image in images file and removes new line
	imageFilenames = [ imageFilename.strip() for imageFilename in imageFilenames ]
		
	# set constraint to default unless one is given on command line
	if not options.constraint:
		constraintFilename = "none"
	elif not os.path.isfile(options.constraint):
		parser.error("constraint file " + options.constraint + " either does not exist or is not accessible")
	else:
		constraintFilename = options.constraint
	
	# set the log header
	logMsg = "run on " + time.strftime("%m-%d-%Y") + "\n"
	
	# this loops through every image in images file and writes the log, running galfit
	for imageFilename in imageFilenames:
		logMsg = logMsg + imageFilename + ": "
		
		# run galfit, preventing crashes but logging errors in log
		try:
			# galfit returns a string indicating success or some failure
			logMsg = logMsg + run_galfit(imageFilename, constraintFilename, options.bulge, options.gauss)
			
		# allow user to stop the program running altogether with ctrl-c
		except KeyboardInterrupt:
			print ("Escape character ctrl-c used to terminate rungalfit. Log file will still be written.")
			break
		
		# log all runtime errors, but move on to the next image regardless
		except:
			errorMsg = str(sys.exc_info()[0]) + str(sys.exc_info()[1])
			print errorMsg
			logMsg = logMsg + errorMsg
		
		# every image on its own line in the log file
		logMsg = logMsg + "\n"
	
	# create the log file in the same directory as python was called
	logFilename = ("rungalfit_log_" + 
					imageFilenames[0].split("/")[-1].split("_")[1].split(".")[1] + "_to_" + 
					imageFilenames[-1].split("/")[-1].split("_")[1].split(".")[1] +
					"_" + time.strftime("%m-%d-%Y") + ".txt")
	
	# write the log file
	print ("writing log file to " + logFilename)
	log = open(logFilename, 'w')
	log.write(logMsg)
	log.close()
	
	print ("Done!\nIn order to summarize results run sumgalfit.py")
######################### done ################################################
