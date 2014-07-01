
#! /usr/bin/python


import os
import sys
from pyraf import iraf
import time
import math
from optparse import OptionParser


def main(imageListFilename, galfit_constraint_filename, psf, mpZeropoint, 
		plateScale, includeBulgeComponent, includeGaussianSmoothing, 
		runSimSextractor, runRealSextractor, sextractorOptionsList):
	'''
	main method loops through all image filenames in image list, running galfit
	and logging errors to a log file, which is named according to the date and 
	the images on which galfit was run.
	Also runs sextractor first if so specified by the command line in order
	to provide galfit with a mask (pixel replaced segmentation file)
	
	parameter imageListFilename -
		the string of the full path filename of the list of images 
		on which galfit will be run
		
	parameter galfit_constraint_filename -
		the filename of the constraint file. if none, value is "none"
		
	parameter psf, mpZeropoint, plateScale - 
		optional command line inputs with defaults for galfit param file
	
	parameter includeBulgeComponent - 
		boolean indicating if a bulge should be fit after first pass by galfit
		
	parameter includeGaussianSmoothing - 
		boolean indicating if a gaussian smoothing should be applied
		
	parameter runSimSextractor, runRealSextractor - 
		booleans for sim or real version of sextractor config file
		
	parameter sextractorOptionsList - 
		list of keyword options for sextractor left over on command line
	'''	
	
	# read list of image filenames from input file
	inputFile = open(imageListFilename, 'r')
	imageFilenames = inputFile.readlines()
	inputFile.close()
	
	if len(imageFilenames) == 0:
		print("no images in given file (positional argument). Exiting.")
		exit()
		
	# convenience boolean for dealing with conditionally running sextractor
	runSextractor = runRealSextractor or runSimSextractor
		
	# sextractor specific tasks that can be done outside of image loop
	if runSextractor:
		
		#not based on images, ".".join(imageFilename.split(".")[:-1]) + ".sex"
		configFilename = "rungalfitSExtractorConfig.sex" 
		
		# these can be changed on the command line with flags
		outputCatFilename = "generic.cat" # -CATALOG_NAME <filename>
		segmentationMapFilename = "check.fits" # -CHECKIMAGE_NAME <filename>
		for sexIndex, sexOpt in enumerate(sextractorOptionsList):
			if sexOpt == "-CATALOG_NAME":
				outputCatFilename = sextractorOptionsList[sexIndex + 1]
			elif sexOpt == "-CHECKIMAGE_NAME":
				segmentationMapFilename = sextractorOptionsList[sexIndex + 1]
		
		# write the sextractor config file, which will be used for all images
		write_sextractor_config_file(configFilename, runSimSextractor)
		
		# prepend the -c <config file> option to the sextractor options list
		sextractorOptionsList = ["-c", configFilename] + sextractorOptionsList
		
		print("running sextractor with options: " + str(sextractorOptionsList))
	else:
		segmentationMapFilename = "none"

	# set the log header
	logMsg = "run on " + time.strftime("%m-%d-%Y") + "\n"
	
	# this loops through every image in images file and removes whitespace
	imageFilenames = [imageFilename.strip() for imageFilename in imageFilenames]
	
	# this loops through every image, runs galfit, and writes the log string 
	for imageFilename in imageFilenames:
		
		# run galfit, preventing crashes but printing and logging errors in log
		try:
					
			# run sextractor (if command line set to do so) for each image
			# galfit uses the resulting segmentation map
			if runSextractor:
				
				# smooth image for sextractor use if set to do so on command line
				if includeGaussianSmoothing:
					sexImageFilename = run_gauss(imageFilename, 15)
				else:
					sexImageFilename = imageFilename
					
				# run sextractor and saved return galaxy id for imreplace
				[gfitGalaxyID, gfitGalaxyDist] = run_sextractor(
					sexImageFilename, outputCatFilename, sextractorOptionsList)
				
				# use imreplace method to zero single galaxy in the segmentation map
				#run_imreplace(segmentationMapFilename, gfitGalaxyID, gfitGalaxyID)
				# rename outputs of sextractor so they will not be overwritten?
				# use sextractor galaxy (x, y) or continue to use minmax?
				os.system("mv " + outputCatFilename + " " +
				#		"/".join(imageFilename.split("/")[:-1]) + "/" +
						".".join(imageFilename.split("/")[-1
									].split(".")[:-1]) + ".cat")
				os.system("mv " + segmentationMapFilename + " " +
				#		"/".join(imageFilename.split("/")[:-1]) + "/" +
						".".join(imageFilename.split("/")[-1
									].split(".")[:-1]) + "_check.fits")
				logMsg = (logMsg + imageFilename + 
						": id dist = " +
						str(gfitGalaxyID) + " " + str(gfitGalaxyDist) + "\n")
				# TODO: so that galfit does not run, remove to run galfit
				continue
		
			# galfit returns a string indicating success or some failure
			logMsg = run_galfit(imageFilename, logMsg, 
							galfit_constraint_filename, psf,
							mpZeropoint, plateScale, segmentationMapFilename, 
							includeBulgeComponent, includeGaussianSmoothing)
		
		# allow user to stop the program running altogether with ctrl-c 
		# (might require multiple escapes)
		except KeyboardInterrupt:
			print ("Escape sequence ctrl-c used to terminate rungalfit. " + 
					"Log file will still be written.")
			break
		
		# catch, log, and ignore all runtime errors except explicit exits 
		# (for debugging). move on to the next image regardless
		except not SystemExit:
			errorMsg = str(sys.exc_info()[0]) + str(sys.exc_info()[1])
			print (errorMsg)
			logMsg = logMsg + errorMsg
		
		# every image on its own line in the log file
		logMsg = logMsg + "\n"
	
	# create the log file in the same directory as python was called
	try:
		logFilename = ("rungalfit_log_" + 
			imageFilenames[0].split("/")[-1].split("_")[1].split(".")[1] + 
			"_to_" + 
			imageFilenames[-1].split("/")[-1].split("_")[1].split(".")[1] +
			"_" + time.strftime("%m-%d-%Y") + ".txt")
	except:
		logFilename = ("rungalfit_log.txt")
	
	# write the log file
	print ("writing log file to " + logFilename)
	log = open(logFilename, 'w')
	log.write(logMsg)
	log.close()
	
	print ("Done!\nIn order to summarize results run sumgalfit.py")

	
def write_sextractor_config_file(sextractor_config_filename, runSimSextractor):
	'''
	writes the sextractor config file
	
	parameter sextractor_config_filename - 
		the filename of the sextractor config file being written
	parameter runSimSextractor - 
		boolean for sim version (real is default) of sextractor config file
	'''
	# variables describing the sextractor config file
	catalogName = "generic.cat" # Name of the output catalogue. If the name "STDOUT" is given and CATALOG TYPE is set to ASCII, ASCII HEAD, ASCII SKYCAT, or ASCII VOTABLE the catalogue will be piped to the standard output (stdout
	catalogType = "ASCII_HEAD" 	# Format of output catalog:
								# ASCII - ASCII table; the simplest, but space and time consuming,
								# ASCII HEAD - as ASCII, preceded by a header containing information about the content,
								# ASCII SKYCAT - SkyCat ASCII format (WCS coordinates required)
								# ASCII VOTABLE - XML-VOTable format, together with meta-data,
								# FITS 1.0 - FITS format as in SExtractor 1
								# FITS LDAC - FITS "LDAC" format (the original image header is copied).
	sextractorParamFilename = "sex.param"
	#------------------------------- Extraction ----------------------------------
	detectType = "CCD" 
	#fitsUnsigned = Force 16-bit FITS input data to be interpreted as unsigned integers.
	#flagImage = File name(s) of the flagimage(s)
	flagType = "OR"	# Combination method for flags on the same object:
					# OR arithmetical OR,
					# AND arithmetical AND,
					# MIN minimum of all flag values,
					# MAX maximum of all flag values,
					# MOST most common flag value.
					
	# these variable values change for simulation template vs real template
	if runSimSextractor:
		detectMinArea = "10000"	# Minimum number of pixels above threshold triggering detection
		detectThreshold = "20"	#Detection threshold (0-2). 1 argument: (ADUs or relative to Background RMS, see THRESH TYPE). 2 arguments: R (mag.arcsec 2 ), Zero-point (mag).
		deblendThreshold = "16"	# Minimum contrast parameter for deblending.
		deblendMinContrast = "0.005"	# Minimum contrast parameter for deblending.
	else:
		detectMinArea = "5"	# Minimum number of pixels above threshold triggering detection
		detectThreshold = "0.75"	#Detection threshold (0-2). 1 argument: (ADUs or relative to Background RMS, see THRESH TYPE). 2 arguments: R (mag.arcsec 2 ), Zero-point (mag).
		deblendThreshold = "16"	# Minimum contrast parameter for deblending.
		deblendMinContrast = "0.0001"
		
	analysisThreshold = "5"	#Threshold (in surface brightness) at which CLASS STAR and FWHM operate. 1 argument: relative to Background RMS. 2 arguments: mu (mag/arcsec 2 ), Zero-point (mag).
	filterBool = "Y"	#  If true,filtering is applied to the data before extraction.
	filterName = "sex.conv"	# Name and path of the file containing the filter definition
	#filterThresh = ? 	# Lower and higher thresholds (in back-ground standard deviations) for a pix-el
						#to be consideredin filtering (used for retinafiltering only).
	#threshType = ?	# Meaning of the DETECT_THRESH and ANALYSIS_THRESH parameters:
						# RELATIVE - scaling factor to the background RMS,
						# ABSOLUTE - absolute level (in ADUs or in surface brightness)
		
	cleanBool = "Y"	# If true, a cleaning of the catalog is done before being written to disk.
	cleanParam = "1.0"	# Efficiency of cleaning.
	maskType = "CORRECT"	#CORRECT - replace by values of pixels symmetric with respect to the source center.
							#BLANK --put detected pixels belonging to neighbors to zero,
							#NONE no masking,
	#------------------------------ Photometry -----------------------------------
	photoApertureDiameter = "10." 	# these threes variables are related to the Kron radius, which is introduced as a accurate flexible aperture
									# that would capture most of the flux from an object.
	magAutoKronFact = "2.5"
	magAutoMinRadius = "3.5"
	saturLevel = "120"  #Pixel value above which it is considered saturated.
	magGamma = "4.0"    #gamma of the emulsion (slope of the response function). Takes effect in PHOTO mode 
						#only but NEEDS to be specified, even for CCD images.
	pixelScale = "0.06" #Pixel size in arcsec.
	#------------------------- Star/Galaxy Separation ----------------------------
	seeingFWHM = "0.18"# ????? Seeing_FWHM FWHM of stellar images in arcsec.This quantity is used only
						#for the neural network star/galaxy separation as expressed in the CLASS STAR output.
	starNNWFilename = "default.nnw" #Name of the file containing the neural network weights for star/galaxy separation.
	#------------------------------ Background -----------------------------------
	backSize = "600"	# Size, or Width, Height (inpixels) of a background mesh.
	backFilterSize = "9"	#Size, or Width, Height (inbackground meshes) of the background-filtering mask.
	backPhotoType = "LOCAL"	#Background used to compute magnitudes:
							# GLOBAL - taken directly from the background map
							# LOCAL - recomputed in a rectangular annulus around the object
	backPhotoThickness = "100"	#Thickness (in pixels) of the background LOCAL annulus.
	backType = ""	# What background is subtracted from the images: 
					# AUTO - The internal interpolated background-map. In the manual it says "INTERNAL" here but the keyword is AUTO. 
					# MANUAL A user-supplied constant value provided in BACK VALUE.
	backValue = ""	# in BACK TYPE MANUAL mode, the constant value to be subtracted from the images.
	#------------------------------ Check Image ----------------------------------
	checkImageName = "check.fits" # CHECKIMAGE NAME - File name for each "check-image"
	checkImageType = "SEGMENTATION" #display patches corresponding to pixels attributed to each object
						#APERTURES-- MAG APER and MAG AUTO integration limits
						#-OBJECTS-- background-subtracted image with detected objects blanked,
						#OBJECTS detected objects,
						#FILTERED background-subtracted filtered image (requires FILTER = Y),
						#-BACKGROUND background-subtracted image,
						#MINIBACK RMS low-resolution background noise map, 
						#MINIBACKGROUND low-resolution background map,
						#BACKGROUND RMS full-resolution interpolated background noise map,
						#BACKGROUND full-resolution interpolated background map,
						#IDENTICAL identical to input image (useful for converting formats),
						#NONE no check-image,
	#--------------------- Memory (change with caution!) -------------------------
	memoryObjStack = "50000"
	memoryPixStack = "1000000"
	memoryBufferSize = "8500"
	#----------------------------- Miscellaneous ---------------------------------
	verboseType = "QUIET"	#run silently,
							#NORMAL display warnings and limited info concerning the work in progress,
							#EXTRA WARNINGS like NORMAL, plus a few more warnings if necessary,
							#FULL display a more complete information and the principal parameters of all the objects extracted.
	#------------------------------- New Stuff -----------------------------------
	weightType = "NONE" 	#variance-map derived from an external weight-map,
							#MAP_VAR external variance-map,
							#MAP RMS variance-map derived from an external RMS-map,
							#BACKGROUND variance-map derived from the image itself,
							#NONE no weighting,
	#WEIGHT IMAGE --File name of the detection and measurement weight image , respectively.
	#WEIGHT GAIN -- If true, weight maps are considered as gain maps.
	magZeropoint = "25.96"     #Zero-point offset to be applied to magnitudes.
	#------------------------- End of variable definitions -----------------------
	
	
	# use above variables to write the config file
	os.system('touch ' + sextractor_config_filename)
	sextractorConfigFile = open(sextractor_config_filename,'w')
	sextractorConfigFile.write(
'''
# Default configuration file for SExtractor V1.2b14 - > 2.0
# EB 23/07/98
# (*) indicates parameters which can be omitted from this config file.

#-------------------------------- Catalog ------------------------------------

''')
	if catalogName:
		#generic.fits
		sextractorConfigFile.write(
			"CATALOG_NAME    " + catalogName + "\n")
	
	if catalogType:
		# ASCII_HEAD
		sextractorConfigFile.write(
			"CATALOG_TYPE    " + catalogType + 
			"      # \"NONE\",\"ASCII_HEAD\",\"ASCII\",\"FITS_1.0\", or \"FITS_LDAC\"\n")
	
	if sextractorParamFilename:
		#WFC3.morphWG.param
		sextractorConfigFile.write(
			"PARAMETERS_NAME " + sextractorParamFilename + 
			"  # name of the file containing catalog contents\n")
		
	sextractorConfigFile.write(
'''

#------------------------------- Extraction ----------------------------------

''')
	
	if detectType:
		#CCD
		sextractorConfigFile.write(
			"DETECT_TYPE     " + detectType + 
			"             # \"CCD\" or \"PHOTO\" (*)\n")
	
	if flagType:
		#OR
		sextractorConfigFile.write(
			"FLAG_TYPE       " + flagType + "\n")
	
	if detectMinArea:
		#5
		sextractorConfigFile.write(
			"DETECT_MINAREA  " + detectMinArea +
			"               # minimum number of pixels above threshold\n")
	
	if detectThreshold:
		#0.75
		sextractorConfigFile.write(
			"DETECT_THRESH   " + detectThreshold +
			"            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n")
	
	if analysisThreshold:
		#5
		sextractorConfigFile.write(
			"ANALYSIS_THRESH " + analysisThreshold +
			"               # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n")
	
	if filterBool:
		#Y
		sextractorConfigFile.write(
			"FILTER          " + filterBool +
			"               # apply filter for detection (\"Y\" or \"N\")?\n")
	
	if filterName:
		#tophat_9.0_9x9.conv
		sextractorConfigFile.write(
			"FILTER_NAME     " + filterName + 
			"     # name of the file containing the filter\n")
	
	if deblendThreshold:
		#16
		sextractorConfigFile.write(
			"DEBLEND_NTHRESH " + deblendThreshold + 
			"              # Number of deblending sub-thresholds\n")
	
	if deblendMinContrast:
		#0.0001
		sextractorConfigFile.write(
			"DEBLEND_MINCONT " + deblendMinContrast +
			"  # Minimum contrast parameter for deblending\n\n")
	
	if cleanBool:
		#Y
		sextractorConfigFile.write(
			"CLEAN           " + cleanBool +
			"               # Clean spurious detections? (Y or N)?\n")
	
	if cleanParam:
		#1.0
		sextractorConfigFile.write(
			"CLEAN_PARAM     " + cleanParam + 
			"             # Cleaning efficiency\n\n")
	
	if maskType:
		#CORRECT
		sextractorConfigFile.write(
			"MASK_TYPE       " + maskType + 
			"         # type of detection MASKing: can be one of\n" +
			"                                # \"NONE\", \"BLANK\" or \"CORRECT\"\n")
	
	sextractorConfigFile.write(
'''

#------------------------------ Photometry -----------------------------------

''')
	
	if photoApertureDiameter:
		#10.
		sextractorConfigFile.write(
			"PHOT_APERTURES   " + photoApertureDiameter +
			"   # MAG_APER aperture diameter(s) in pixels\n")
	
	if magAutoKronFact and magAutoMinRadius:
		#2.5, 3.5
		sextractorConfigFile.write(
			"PHOT_AUTOPARAMS " + magAutoKronFact + ", " + magAutoMinRadius +
			"        # MAG_AUTO parameters: <Kron_fact>,<min_radius>\n")
	
	if saturLevel:
		#120
		sextractorConfigFile.write(
			"SATUR_LEVEL     " + saturLevel + 
			"             # level (in ADUs) at which arises saturation\n")
	
	if magGamma:
		#4.0
		sextractorConfigFile.write(
			"MAG_GAMMA       " + magGamma +
			"             # gamma of emulsion (for photographic scans)\n")
	
	if pixelScale:
		#0.06
		sextractorConfigFile.write(
			"PIXEL_SCALE     " + pixelScale +
			"            # size of pixel in arcsec (0=use FITS WCS info).\n")
	
	sextractorConfigFile.write(
'''

#------------------------- Star/Galaxy Separation ----------------------------

''')
	
	if seeingFWHM:
		#0.18
		sextractorConfigFile.write(
			"SEEING_FWHM     " + seeingFWHM + 
			"            # stellar FWHM in arcsec\n")
	
	if starNNWFilename:
		#default.nnw
		sextractorConfigFile.write(
			"STARNNW_NAME    " + starNNWFilename + 
			"     # Neural-Network_Weight table filename\n")
	
	sextractorConfigFile.write(
'''

#------------------------------ Background -----------------------------------

''')
	
	if backSize:
		#256
		sextractorConfigFile.write(
			"BACK_SIZE       " + backSize + 
			"             # Background mesh: <size> or <width>,<height>\n")
	
	if backType:
		#not in original
		sextractorConfigFile.write(
			"BACK_TYPE       " + backType + 
			"             # What background is subtracted from the images:" +
'''
				# AUTO - The internal interpolated background-map. In the manual it says "INTERNAL" here but the keyword is AUTO. 
				# MANUAL A user-supplied constant value provided in BACK VALUE.")
''')
	
	if backValue:
		#not in original
		sextractorConfigFile.write(
			"BACK_VALUE       " + backValue + 
			"             # in BACK TYPE MANUAL mode, the constant value to be subtracted from the images.\n")
	
	if backFilterSize:
		#9
		sextractorConfigFile.write(
			"BACK_FILTERSIZE " + backFilterSize +
			"               # Background filter: <size> or <width>,<height>\n\n")
	
	if backPhotoType:
		#LOCAL
		sextractorConfigFile.write(
			"BACKPHOTO_TYPE  " + backPhotoType + 
			"           # can be \"GLOBAL\" or \"LOCAL\" (*)\n")
	
	if backPhotoThickness:
		#100
		sextractorConfigFile.write(
			"BACKPHOTO_THICK " + backPhotoThickness + 
			"             # thickness of the background LOCAL annulus (*)\n")
	
	sextractorConfigFile.write(
'''

#------------------------------ Check Image ----------------------------------

''')
	
	if checkImageName:
		#check.fits
		sextractorConfigFile.write(
			"CHECKIMAGE_NAME  " + checkImageName + 
			"	# CHECKIMAGE NAME - File name for each \"check-image\"\n")
	
	if checkImageType:
		#SEGMENTATION
		sextractorConfigFile.write(
			"CHECKIMAGE_TYPE  " + checkImageType + 
			"   # can be one of \"NONE\", \"BACKGROUND\"," + 
'''
				# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
				# "-OBJECTS", "SEGMENTATION", "APERTURES",
				# or "FILTERED" (*)
''')
	
	sextractorConfigFile.write(
'''

#--------------------- Memory (change with caution!) -------------------------

# MEMORY_OBJSTACK   50000       # number of objects in stack
# MEMORY_PIXSTACK   1000000     # number of pixels in stack
# MEMORY_BUFSIZE    8500        # number of lines in buffer

''')
	
	if memoryObjStack:
		#2000
		sextractorConfigFile.write(
			"MEMORY_OBJSTACK " + memoryObjStack +
			"            # number of objects in stack\n")
	
	if memoryPixStack:
		#200000
		sextractorConfigFile.write(
			"MEMORY_PIXSTACK " + memoryPixStack + 
			"          # number of pixels in stack\n")
	
	if memoryBufferSize:
		#2048
		sextractorConfigFile.write(
			"MEMORY_BUFSIZE  " + memoryBufferSize + 
			"            # number of lines in buffer\n")
	
	sextractorConfigFile.write(
'''
 
#----------------------------- Miscellaneous ---------------------------------

''')
	
	if verboseType:
		#QUIET
		sextractorConfigFile.write(
			"VERBOSE_TYPE    " + verboseType + 
			"           # can be \"QUIET\", \"NORMAL\" or \"FULL\" (*)\n")
	
	sextractorConfigFile.write(
'''

#------------------------------- New Stuff -----------------------------------

''')
	
	if weightType:
		#MAP_WEIGHT
		sextractorConfigFile.write(
			"WEIGHT_TYPE   " + weightType + "\n")
	
	sextractorConfigFile.write(
'''
# WEIGHT_TYPE     MAP_RMS,MAP_RMS
# WEIGHT_THRESH   10000.0


''')
	
	if magZeropoint:
		#25.96
		sextractorConfigFile.write(
			"MAG_ZEROPOINT   " + magZeropoint + 
			"           # H-band magnitude zero-point\n")
	

def run_sextractor(imageFilename, outputCatFilename, sextractorOptionsList):
	'''
	runs sextractor on the given image, returning the galaxy id of
	the galaxy closest to the center of the image. Also produces
	generic.cat and check.fits in the calling directory, which 
	are galaxy info and the segmentation map, respectively.
	
	parameter imageFilename -
		the full path filename of the image on which sextractor will be run
	
	parameter outputCatFilename -
		the name of the catalog output from running sextractor, 
		as specified by the config file unless on in options list
		
	parameter sextractorOptionsList - 
		list of options for sextractor
		
	returns - the [ID, dist] of the galaxy closest to the center of the image
	'''

	# SYNTAX: sex <image> [<image2>][-c <configuration_file>][-<keyword> <value>]
	os.system("sex " + imageFilename + " " + " ".join(sextractorOptionsList))
	
	# get galaxyID of the galaxy closest to center of image from outputCatFile
	outputCatFile = open(outputCatFilename, 'r')
	outputCatContents = outputCatFile.readlines()
	outputCatFile.close()
	
	# start on the last line
	lineIndex = -1
	
	# set to a high value initially, just needs to be bigger than best
	prevBestDist = 1000.0
	
	# value to indicate if no galaxy is found closer than above distance
	prevBestID = -1
	
	# read backwards until comment character indicates end of galaxies
	while outputCatContents[lineIndex].strip()[0] != "#":
		
		# TODO: might be able to have indices collected from commented header
		# gather galaxy information
		galaxyOutputList = outputCatContents[lineIndex].strip().split()
		galaxyID = galaxyOutputList[0]
		galaxyXimage = galaxyOutputList[26]
		galaxyYimage = galaxyOutputList[27]
		
		# TODO: generalize center
		# compute distance from image center 
		dx = 300.0 - float(galaxyXimage)
		dy = 300.0 - float(galaxyYimage)
		curDist = math.sqrt( dx*dx + dy*dy )
		
		# store galaxyID if closer than prev closest galaxy to center of image
		if curDist < prevBestDist:
			prevBestID = galaxyID
			prevBestDist = curDist
		
		# go to the next (prev in file) line
		lineIndex = lineIndex - 1
		
	if prevBestID == -1:
		print ("sextractor did not yield a galaxy closer than " + 
				str(prevBestDist))
	else:
		print ("closest galaxy id is " + str(prevBestID) + 
				" at distance of " + str(prevBestDist))
		
	# return the id of the galaxy that was identified as the center galaxy
	return [prevBestID, prevBestDist]


def run_galfit(imageFilename, logMsg, galfit_constraint_filename, 
				psf, mpZeropoint, plateScale, segmentationMapFilename, 
				includeBulgeComponent, includeGaussianSmoothing):
	'''
	opens file (parameter) containing list of images and loops over every image, 
	running galfit and writing the results to the same directory as the images
	
	calls methods to invoke iraf methods to write the initial galfit param file
	for each image
	
	parameter imageFilename -
		the full path filename of the image on which galfit will be run
		
	parameter logMsg - 
		the current state of the log at the time of invocation
		so that this method can append to and return it to caller
		
	parameter galfit_constraint_filename -
		the filename of the constraint file. if none, value is "none"
		
	parameter psf, mpZeropoint, plateScale - 
		optional command line inputs with defaults for galfit param file
		
	parameter segmentationMapFilename -
		the filename of the segmentation map that masks galfit
	
	parameter includeBulgeComponent - 
		boolean indicating if a bulge should be fit after first pass by galfit
		
	parameter includeGaussianSmoothing - 
		boolean indicating if a gaussian smoothing should be applied
		
	returns - updated log message indicating success or failure
	'''
	
	# sigma [# pixels] is used for run gauss when blurring
	sigma = 15
	
	# run iraf's imhead method to get image information
	[directory_location, galaxy_id, filt, cam_number, image_number, height, width
		] = run_imhead(imageFilename)

	# the x y location of the top left corner of the area on which 
	# to run minmax and imexam in the coordinate system of the original image 
	# (which would be 0, 0 for the original)
	xStart = float(width)/2 - 75.0
	xStop = float(width)/2 + 75.0
	yStart = float(height)/2 - 75.0
	yStop = float(height)/2 + 75.0
	
	# calls iraf's minmax method, passing image filename as a parameter
	# with image possible gaussian smoothed, as well as two points
	# on the image defining the area on which to run minmax 
	# returns as a list of two strings (the coordinates of the max pixel)
	if includeGaussianSmoothing:
		centerCoords = run_minmax(run_gauss(imageFilename, sigma), 
									xStart, yStart, xStop, yStop)
	else:
		centerCoords = run_minmax(imageFilename, 
									xStart, yStart, xStop, yStop)
	
	# calls iraf's imexam method, passing filename as a parameter
	# along with center coordinates and two points defining area
	# returns initial estimates of the returned paramters
	[Y, X, magnitude, rad, BA, angle
		] = run_imexam(imageFilename, centerCoords, xStart, yStart, xStop, yStop)
	
	# write to the log if default values are being used
	if float(BA) == 1.0:
		logMsg = logMsg + "Default b/a used. "
	if float(angle) == 0.0:
		logMsg = logMsg + "Default position angle used. "
	
	# transform coordinates back into coordinates for original image
	Y = str(float(Y) + yStart)
	X = str(float(X) + xStart)
	
	# define filenames
	filename = (directory_location + galaxy_id + "_" + 
				image_number + '_cam' + str(cam_number) + '_' + filt)
	galfit_single_parameter_filename =	filename + '_single_param.txt'
	galfit_single_output_filename =		filename + "_single_multi.fits"
	galfit_single_result_filename =		filename + "_single_result.txt"
	galfit_bulge_parameter_filename =	filename + '_bulge_param.txt'
	galfit_bulge_output_filename =		filename + "_bulge_multi.fits"
	galfit_bulge_result_filename =		filename + "_bulge_result.txt"
							
	# writes the single component parameter file, given filenames and galxy parameters
	write_galfit_single_parameter(imageFilename, galfit_single_parameter_filename,
									galfit_single_output_filename, 
									psf, segmentationMapFilename,
									mpZeropoint, plateScale, 
									1, 1, width, height, 
									X, Y, magnitude, rad, BA, angle)

	# run galfit on paramter file
	os.system('galfit ' + galfit_single_parameter_filename)

	# detects atomic galfit error
	# If not failure then removes temp files, otherwise returns
	if os.path.isfile(galfit_single_output_filename):
		# remove temp files
		os.system("rm " + "fit.log")
		os.system("mv galfit.01 " + galfit_single_result_filename)
	else:
		return (logMsg + "galfit failed on single component, " + 
							"probably mushroom (atomic galfit error)")
	
	# done unless command line specified that a second galfit run
	# should be done by adding a bulge component to the result of the first run
	if includeBulgeComponent:
	
		# reads the results of the first run and returns it as a long string
		# with the output and constraint modified for bulge run and the
		# new bulge component appended with some intitial guess parameters
		bulgeParamStr = get_galfit_bulge_parameter_str(galfit_single_result_filename, 
							galfit_bulge_output_filename, galfit_constraint_filename, rad)

		# write the bulge parameter file using the modified contents of the single results
		galfitBulgeParamFile = open(galfit_bulge_parameter_filename, "w")
		galfitBulgeParamFile.write(bulgeParamStr)
		galfitBulgeParamFile.close()
		
		# run galfit on paramter file
		os.system('galfit ' + galfit_bulge_parameter_filename)
	
		# detects atomic galfit error
		# If not failure then removes temp files, otherwise returns
		if os.path.isfile(galfit_bulge_output_filename):
			# remove temp files
			os.system("rm " + "fit.log")
			os.system("mv galfit.01 " + galfit_bulge_result_filename)
			#os.system('mv ' + galfit_output_filename + ' ' + galfit_bulge_output_filename)
		else:
			return (logMsg + "galfit failed on bulge component, " + 
								"probably mushroom (atomic galfit error)")
	
	# if we get here then nothing went wrong!
	return logMsg + "success"
	

def run_imhead(imageFilename):
	'''
	invokes the iraf method imhead on the given image filename
	parses the result to be returned as a list of string values
	
	parameter imageFilename - 
		the full path filename of the image on which to invoke imexam
	
	returns - 
		a string list of image info in the form 
		[galaxy_id, filt, cam_number, image_number, frame_height, frame_width]
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
	
	filt = imhead_info.split("_")[6]
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
	
	return [directory_location, galaxy_id, filt, cam_number, image_number, frame_height, frame_width]
	
	
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
	
	
def run_minmax(imageFilename, xStart, yStart, xStop, yStop):
	'''
	method calls iraf's minmax method on the area defined by the four
	parameter values that define two points on the image
	start is the top left, stop is the bottom right
	
	parameter imageFilename - 
		the full path filename of the image on which to invoke minmax
	
	returns - list of two strings defining the coordinates of the max pixel 
				['xMax', 'yMax']
	'''
	
	# string representing the part of the image to run minmax and imexam on
	areaStr = ('[' + str(int(xStart)) + ':' + str(int(xStop)) + "," +
					 str(int(yStart)) + ':' + str(int(yStop)) + ']')

	#call iraf's minmax function
	return iraf.minmax(imageFilename + areaStr, Stdout=1, update=0, force=1
		)[0].strip().split(" ")[3].replace("[", "").replace("]", "").split(",")
		
		
def run_imexam(imageFilename, centerCoords, xStart, yStart, xStop, yStop):
	'''
	method calls iraf's imexam method on the area defined by the four
	
	parameter imageFilename - 
		the full path filename of the image on which to invoke imexam
	parameter centerCoords - 
		the coordinates of the initial guess of the center of the galaxy
	parameter xStart, xStop, yStart, yStop - 
		two points on the original image defining an area to run minmax on
		start is the top left, stop is the bottom right
	
	returns - list of string image info in the form 
				[line, col, mag, radius, b_a, PA]
	'''
		
	
	# string representing the part of the image to run minmax and imexam on
	areaStr = ('[' + str(int(xStart)) + ':' + str(int(xStop)) + "," +
					 str(int(yStart)) + ':' + str(int(yStop)) + ']')

	# writes the output of the minmax function of iraf to a file for later use. 
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
	
	# run imexam passing coords.tmp, 
	# data is returned in the last two elements of the return array of strings
	imexam_out = iraf.imexam(imageFilename + areaStr, 
						use_display=0, imagecur=coordsFilename, Stdout=1)[-2:]

	# ['COL	LINE X Y', 
	# 'R MAG FLUX SKY PEAK E PA BETA ENCLOSED MOFFAT DIRECT']
	
	# delete coords.tmp, no longer needed
	os.system('rm ' + coordsFilename) 
	
	xy = imexam_out[0].strip().split()
	data = imexam_out[1].strip().split()
	
	col = float(xy[0])							#this gives x value 
	line = float(xy[1])							#this gives y value
	radius = float(data[0])						#this gives the radius
	mag = float(data[1])						#this gives the magnitude
	
	# these might be indef, if so then set to default values
	try:
		b_a = 1.0 - float(data[5])				# b/a=1-e
	except:
		print("using default values for b/a")
		b_a = 1.0
	
	#this gives the position angle. iraf measures up from x. Galfit down from y
	try:
		PA = float(data[6]) - 90.0				
	except:
		print("using default values for position angle")
		PA = 0.0
	
	return [line, col, mag, radius, b_a, PA]
	
	
def run_imreplace(imageFilename, lowPixVal, uppPixVal):
	'''
	replace all pixels in image between lower and upper with zero

	parameter imageFilename - the full path filename of the image on which to invoke imexam
	parameter lowPixVal - the lower bound of pixels to be replaced with zero
	parameter uppPixVal - the upper bound of pixels to be replaced with zero
	'''
	# TODO: this could be inlines easily, only called once in main()
	iraf.imreplace(imageFilename, 0, lower=lowPixVal, upper=uppPixVal)


def write_galfit_single_parameter(imageFilename, 
									galfit_single_parameter_filename,
									galfit_single_output_filename, 
									psf, segmentationMapFilename,
									mpZeropoint, plateScale, 
									xMin, yMin, xMax, yMax, 
									X, Y, magnitude, rad, BA, angle):
	'''
	writes the galfit parameter file for the single component
	
	parameter imageFilename - 
		the filename of the image on which galfit is being run
	parameter galfit_single_parameter_filename - 
		the filename of the galfit parameter file being written
	parameter galfit_single_output_filename - 
		the filename where the output of the galfit run will be written
	parameter psf - the filename of the psf to give galfit
	parameter segmentationMapFilename - 
		the filename of the Bad pixel mask to give galfit
	parameter xMin, yMin, xMax, yMax - 
		the area on the image on which to run galfit
	parameter X, Y, magnitude, rad, BA, angle - 
		initial guess at the description of the single galaxy component
	'''
	os.system('touch ' + galfit_single_parameter_filename)
	galfitSingleParamFile = open(galfit_single_parameter_filename,'w')
	galfitSingleParamFile.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
	galfitSingleParamFile.write("A) " + imageFilename + 
		"						#Input data image block\n")
	galfitSingleParamFile.write("B) " + galfit_single_output_filename +
		"						#Output data image block\n")
	galfitSingleParamFile.write("C)" + " none" + 
		"						#Sigma image name (made from data if blank or 'none')\n")
	galfitSingleParamFile.write("D) " + psf + 
		"						#Input PSF image and (optional) diffusion kernel\n")
	galfitSingleParamFile.write("E)" + " 1" + 
		"						#PSF fine sampling factor relative to data\n")
	galfitSingleParamFile.write("F) " + segmentationMapFilename +							
		"						#Bad pixel mask (FITS file or ASCIIcoord list)\n")
	galfitSingleParamFile.write("G)" + " none" + 
		"						#File with parameter constraints (ASCII file)\n")
	galfitSingleParamFile.write("H) " + str(xMin) + " " + str(xMax) + " " + 
										str(yMin) + " " + str(yMax) + 
		"						#Image region to fit (xmin xmax ymin ymax)\n")
	galfitSingleParamFile.write("I)" + " 200 200" + 
		"						#Size of the convolution box (x y)\n")
	galfitSingleParamFile.write("J) " + str(mpZeropoint) + 
		"						#Magnitude photometric zeropoint\n")
	galfitSingleParamFile.write("K) " + str(plateScale) + "	 " + str(plateScale) + 
		"						#Plate scale (dx dy)  [arcsec per pixel]\n")
	galfitSingleParamFile.write("O)" + " regular" + 
		"						#display type (regular, curses, both\n")
	galfitSingleParamFile.write("P)" + " 0" + 
		"						#Options: 0=normal run; 1,2=make model/imgblock & quit\n")
	galfitSingleParamFile.write("S)" + " 0" + 
		"						#Modify/create components interactively?\n")
	galfitSingleParamFile.write("\n")
	galfitSingleParamFile.write("# INITIAL FITTING PARAMETERS\n")
	galfitSingleParamFile.write("#\n")
	galfitSingleParamFile.write("# For component type, the allowed functions are:\n")
	galfitSingleParamFile.write("#		nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n")
	galfitSingleParamFile.write("#		ferrer, coresersic, sky and isophote.\n")
	galfitSingleParamFile.write("#\n")
	galfitSingleParamFile.write("# Hidden parameters will only appear when they are specified:\n")
	galfitSingleParamFile.write("#		C0 (diskyness/boxyness),\n")
	galfitSingleParamFile.write("#		Fn (n=interger, Azimuthal Fourier Modes).\n")
	galfitSingleParamFile.write("#		R0-R10 (PA rotation, for creating spiral structures).\n")
	galfitSingleParamFile.write("#\n")
	galfitSingleParamFile.write("# ------------------------------------------------------------------------------\n")
	galfitSingleParamFile.write("#		par)	par value(s)	fit toggle(s)	# parameter description\n")
	galfitSingleParamFile.write("# ------------------------------------------------------------------------------\n")
	galfitSingleParamFile.write("\n")
	galfitSingleParamFile.write("# Componenet number: 1\n")
	galfitSingleParamFile.write(" 0) sersic					#Component type\n")
	galfitSingleParamFile.write(" 1) " + str(X) + "	" + str(Y) + "	1	1			#Position x,y\n")
	galfitSingleParamFile.write(" 3) " + str(magnitude) + "	1			#Integrated Magnitude\n")
	galfitSingleParamFile.write(" 4) " + str(rad) + "			1			#R_e (half-light radius)	[pix]\n")
	galfitSingleParamFile.write(" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
	galfitSingleParamFile.write(" 6) 0.0000		0			#	-----\n")
	galfitSingleParamFile.write(" 7) 0.0000		0			#	-----\n")
	galfitSingleParamFile.write(" 8) 0.0000		0			#	-----\n")
	galfitSingleParamFile.write(" 9) " + str(BA) + "			1			#Axis ratio (b/a)\n")
	galfitSingleParamFile.write(" 10) " + str(angle) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
	galfitSingleParamFile.write(" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
	galfitSingleParamFile.write("\n")


	galfitSingleParamFile.write("# Componenet number: 2\n")
	galfitSingleParamFile.write(" 0) sky						#Component type\n")
	galfitSingleParamFile.write(" 1) 0.0000		0			#	Sky background at center of fitting region [ADUs]\n")
	galfitSingleParamFile.write(" 2) 0.0000		0			#	dsky/dx (sky gradient in x) [ADUs/pix]\n")
	galfitSingleParamFile.write(" 3) 0.0000		0			#	dsky/dy (sky gradient in y) [ADUs/pix]\n")
	galfitSingleParamFile.write(" Z) 0						#Leave in [1] or subtract [0] this comp from data?\n")
	galfitSingleParamFile.write("\n")

	galfitSingleParamFile.close()
	
	
def get_galfit_bulge_parameter_str(galfit_single_result_filename, 
				galfit_bulge_output_filename, galfit_constraint_filename, rad):
	'''
	reads the results of the first run and returns it as a long string
	with the output and constraint modified for bulge run and the
	new bulge component appended
	
	parameter galfit_single_result_filename - the parameter file being read
	parameter galfit_bulge_output_filename - the image output filename
	parameter galfit_constraint_filename - the constraint filename
	parameter rad - the half light radius from original imexam estimate
	
	returns - the string of the result file given, 
			modified to be used as parameter file for two component fit 
	'''
	singleResultFile = open(galfit_single_result_filename, "r")		
	resultContents = singleResultFile.readlines()
	singleResultFile.close()
	
	# gather info about first fit while modifying output and coonstraint parameters
	bulgeParamStr = ""		
	positionLine = ""
	magnitudeLine = ""
	for resultLine in resultContents:
	
		if resultLine.strip()[:2] == "B)":
			bulgeParamStr = (bulgeParamStr + "B) " + galfit_bulge_output_filename +
						"						#Output data image block\n")
			
		elif resultLine.strip()[:2] == "G)":
			bulgeParamStr = (bulgeParamStr + "G) " + galfit_constraint_filename + 
						"						#File with parameter constraints (ASCII file)\n")
		
		elif not positionLine and resultLine.strip()[:2] == "1)":
			positionLine = resultLine.strip().split(" ")
			bulgeParamStr = bulgeParamStr + resultLine
			
		elif not magnitudeLine and resultLine.strip()[:2] == "3)":
			magnitudeLine = resultLine.strip().split(" ")
			bulgeParamStr = bulgeParamStr + resultLine
			
		else:
			bulgeParamStr = bulgeParamStr + resultLine
			
	# store results of single component run into variables
	resultX = positionLine[1]
	resultY = positionLine[2]
	resultMag = magnitudeLine[1]
	
	
	# write the third component to the end of galfit_single_result_filename
	bulgeParamStr = bulgeParamStr + ("# Componenet number: 3\n")
	bulgeParamStr = bulgeParamStr + (" 0) sersic					#Component type\n")
	bulgeParamStr = bulgeParamStr + (" 1) " + resultX + " " + resultY + " 1	1			#Position x,y\n")
	bulgeParamStr = bulgeParamStr + (" 3) " + resultMag + "	1			#Integrated Magnitude\n")
	bulgeParamStr = bulgeParamStr + (" 4) " + str(rad) + "			1			#R_e (half-light radius)	[pix]\n")
	bulgeParamStr = bulgeParamStr + (" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
	bulgeParamStr = bulgeParamStr + (" 6) 0.0000		0			#	-----\n")
	bulgeParamStr = bulgeParamStr + (" 7) 0.0000		0			#	-----\n")
	bulgeParamStr = bulgeParamStr + (" 8) 0.0000		0			#	-----\n")
	bulgeParamStr = bulgeParamStr + (" 9) 1" + "			1			#Axis ratio (b/a)\n")
	bulgeParamStr = bulgeParamStr + (" 10) 0" + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
	bulgeParamStr = bulgeParamStr + (" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
	bulgeParamStr = bulgeParamStr + ("\n")

	return bulgeParamStr
		

if __name__ == "__main__":
	'''
	parses the command line using the optparse package
	NOTE: argparse is the current version as of python 2.7, 
			but optparse is used to maintain better backwards compatibility
	'''
	
	usage = ("\n%prog inputFile [-h help] [options (with '-'|'--' prefix)]" +
			"  [sextractor options (no '-' prefix)]\n" +
			"Ex:\n" +
			"%prog images.txt\n" +
			"%prog images.txt -c pos.constraint -b -g\n" +
			"%prog images.txt -s parameter_name sex.param filter_name sex.conv")
			
	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# bulge is a boolean (true or false) specifying if the simulation should
	# fit an additional component after the initial fit from imexam results
	parser.add_option("-b","--bulge", 
				help="turn on to include a bulge fit after the initial galaxy fit",
				action="store_true")
						
	# gauss is a boolean (true or false) specifying if the simulation should
	# apply gaussian smoothing before running sextractor and/or iraf.minmax
	parser.add_option("-g","--gauss", 
				help="turn on to include pre gaussian smoothing of image",
				action="store_true")
						
	# run sextractor
	parser.add_option("-s","--simSextractor", 
				help="turn on to run sextractor for simulation images",
				action="store_true")
	parser.add_option("-r","--realSextractor", 
				help="turn on to run sextractor for real images",
				action="store_true")
	
	# the constraint file. verified after parsing
	parser.add_option("-c","--constraint", 
				help="set the file containing the galfit constraints")
			
	# Magnitude photometric zeropoint	
	parser.add_option("--mpz", metavar="MagnitudePhotometricZeropoint",
				type="float", default=26.3, 
				help="set the magnitude photometric zeropoint for" +
					" galfit to use [default: %default]")
						
	# Plate scale
	parser.add_option("--plate", metavar="PlateScale",
				type="float", default=0.06, 
				help="set the plate scale for galfit to use" +
						"[default: %default]")
								
	# PSF file. verified after parsing
	parser.add_option("--psf",
				help="set the file for galfit to use as a PSF")
						
	# parse the command line using above parameter rules
	# options - list with everthing defined above, 
	# args - anything left over after parsing options
	[options, args] = parser.parse_args()
	
	# verify that there is at least one positional argument
	if len(args) < 1:
		parser.error("incorrect number of arguments, must provide an input" + 
					" file containing the list of full path image filenames.")
	# verify that the positional argument is an existing file
	elif not os.path.isfile(args[0]):
		parser.error("input file " + args[0] + 
					" either does not exist or is not accessible")
		
	# for sextractor, if enabled on the command line
	if options.simSextractor and options.realSextractor:
		parser.error("options -s (--simSextractor) and -r (--realSextractor)" +
					" are mutually exclusive")
	sextractorOptionsList = []
	if (options.simSextractor or options.realSextractor) and len(args) > 1:
		for index, sexOpt in enumerate(args[1:]):
			if index % 2 == 0: 	# odd, consider a keyword
				sextractorOptionsList.append("-" + sexOpt.upper())
			else:				# even, consider the value for previous keyword
				sextractorOptionsList.append(sexOpt)
			
		print("command line extra input: " + " ".join(sextractorOptionsList) + 
			"\nAbove are being interpreted as sextractor keyword options")
		
	# galfit constraint default none unless one is given on command line
	if not options.constraint:
		constraintFilename = "none"
	# verify that file exists
	elif not os.path.isfile(options.constraint):
		parser.error("constraint file " + options.constraint + 
					" either does not exist or is not accessible")
	else:
		constraintFilename = options.constraint
		
	# galfit psf default none unless one is given on command line
	if not options.psf:
		psf = "none"
	# verify that the file exists
	elif not os.path.isfile(options.psf):
		parser.error("psf file " + options.psf +
					" either does not exist or is not accessible")
	else:
		psf = options.psf

	# pass pertinent (and verified) info to the main method
	main(args[0], constraintFilename, psf, options.mpz, options.plate, 
		options.bulge, options.gauss, options.simSextractor, 
		options.realSextractor, sextractorOptionsList)