#! /usr/bin/python

import os
import sys
from pyraf import iraf
import time
import math
from optparse import OptionParser

'''
class SimGalaxy:
	
	class to hold the methods and data for a particular 
	simulated galaxy
	
	
	def __init__(self, id):
		'''galaxy constructor'''
		
		# The unique identifier of the galaxy, e.g. "VELA02"
		self.id = id
		
		# a dictionary where the galaxies time steps are the keys
		# and the values are a list of the particular images
		self.timeSteps = {}
'''
class SimImage:
	'''
	class to hold the methods and data for a particular 
	simulated image of a galaxy
	'''
	
	def __init__(self, filename, galaxyID="", timeStep="", filter="", 
					camera="", height=0, width=0, model=None):
		'''image constructor'''
		
		# the full path filename of the image
		self.filename = filename
		
		# the galaxy that the image is of
		self.galaxyID = galaxyID
		
		# the time step of the image
		self.timeStep = timeStep
		
		# the filter of the image, e.g. "F160W"
		self.filter = filter
		
		# the camera of the image, e.g. "0"
		self.camera = camera
		
		# the height of the image, e.g. 600
		self.height = height
		
		# the width of the image, e.g. 600
		self.width = width
		
		# the model of the galaxy in the image
		self.model = model

class SimModel:
	'''
	class to hold the methods and data for a particular 
	model of a simulated galaxy from an image
	'''
	
	def __init__(self, centerCoords=[0.0,0.0], magnitude=0.0, 
					radius=0.0, ba=0.0, angle=0.0):
		'''model constructor'''
		
		# the center coordinates of the galaxy
		self.centerCoords = centerCoords
	
		# the magnitude of the galaxy
		self.magnitude = magnitude
	
		# the half light radius of the galaxy
		self.radius = radius
	
		# the b/a of the galaxy
		self.ba = ba
	
		# the position angle of the galaxy
		self.angle = angle
		
class SimulationProcessor:
	'''
	The controller class, which holds methods for analyzing simulations
	'''
		
	def __init__(self):
		'''
		
		'''
		self.outputCatFilename = "generic.cat" # -CATALOG_NAME <filename>
		self.segmentationMapFilename = "check.fits" # -CHECKIMAGE_NAME <filename>
		self.sextractorOptionsList = []
		self.sextractorConfigFilename = "config.sex"
		
		
	def parseGalfitOptions(self, parser, options):
		'''
		
		'''
		self.mpZeropoint = options.mpz
		self.plateScale = options.plate
		self.includeBulgeComponent = options.bulge
		self.includeGaussianSmoothing = options.gauss
		# galfit constraint default none unless one is given on command line
		if not options.constraint:
			self.constraintFilename = "none"
		# verify that file exists
		elif not os.path.isfile(options.constraint):
			parser.error("constraint file " + options.constraint + 
						" either does not exist or is not accessible")
		else:
			self.constraintFilename = options.constraint
			
		# galfit psf default none unless one is given on command line
		if not options.psf:
			self.psf = "none"
		# verify that the file exists
		elif not os.path.isfile(options.psf):
			parser.error("psf file " + options.psf +
						" either does not exist or is not accessible")
		else:
			self.psf = options.psf
		
		
	def parseSextractorOptions(self, parser, simSextractor, realSextractor, sextractorOptions):
		'''
		
		'''
		self.simSextractor = simSextractor
		self.realSextractor = realSextractor
		if self.simSextractor or self.realSextractor:
			if len(sextractorOptions) % 2 != 0:
				parser.error("Sextractor options parsing error, odd number of options. " + 
					"Expecting sextractor options to be <keyword> <value> pairs")
			else:
				try:
					self.sextractorOptionsList = []
					for sexIndex, sexOpt in enumerate(sextractorOptions):
						if sexIndex % 2 == 0: 	# even, consider a keyword
							self.sextractorOptionsList.append("-" + sexOpt.upper())
							if sexOpt.upper() == "CATALOG_NAME":
								self.outputCatFilename = sextractorOptions[sexIndex + 1]
							elif sexOpt.upper() == "CHECKIMAGE_NAME":
								self.segmentationMapFilename = sextractorOptions[sexIndex + 1]
						else:				# even, consider the value for previous keyword
							self.sextractorOptionsList.append(sexOpt)
				except:
					parser.error("Sextractor options parsing error. " + 
						"Expecting sextractor options to be <keyword> <value> pairs")
						
			print("command line extra input: " + " ".join(self.sextractorOptionsList) + 
				"\nAbove are being interpreted as sextractor keyword options")
		else:
			self.segmentationMapFilename = "none"
			
		
	def write_sextractor_config_file(self):
		'''
		writes the sextractor config file
		
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
		if self.simSextractor:
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
		os.system('touch ' + self.sextractorConfigFilename)
		sextractorConfigFile = open(self.sextractorConfigFilename,'w')
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
	
	
	def run_sextractor(image):
		'''
		Updates image with model containing galaxy position and 
		updates seg map to remove that galaxy by
		running sextractor on the given image, producing
		*.cat and *.fits in the calling directory, which 
		are galaxy info and the segmentation map, respectively.
		
		parameter image -
			the image on which sextractor will be run
		'''
	
		# SYNTAX: sex <image> [<image2>][-c <configuration_file>][-<keyword> <value>]
		os.system("sex " + image.filename + " -c " + self.sextractorConfigFilename + 
					" " + " ".join(self.sextractorOptionsList))
		
		# get galaxyID of the galaxy closest to center of image from outputCatFile
		outputCatFile = open(self.outputCatFilename, 'r')
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
			galaxyXimage = float(galaxyOutputList[26])
			galaxyYimage = float(galaxyOutputList[27])
			
			# TODO: generalize center
			# compute distance from image center 
			dx = 300.0 - galaxyXimage
			dy = 300.0 - galaxyYimage
			curDist = math.sqrt( dx*dx + dy*dy )
			
			# store galaxyID if closer than prev closest galaxy to center of image
			if curDist < prevBestDist:
				prevBestID = galaxyID
				prevBestDist = curDist
				prevBestX = galaxyXimage
				prevBestY = galaxyYimage
			
			# go to the next (prev in file) line
			lineIndex = lineIndex - 1
			
		if prevBestID == -1:
			print ("sextractor did not yield a galaxy closer than " + 
					str(prevBestDist))
		else:
			print ("closest galaxy id is " + str(prevBestID) + 
					" at distance of " + str(prevBestDist))
			
		# store (x,y) of the galaxy as the image's model
		image.model = SimModel(centerCoords=[prevBestX,prevBestY])
		
		# use the id of the galaxy that was identified as the center galaxy
		# to replace the appropriate pixels in segmentation map for galfit mask
		iraf.imreplace(self.segmentationMapFilename, 0, prevBestID, prevBestID)


	def run_imhead(self, image):
		'''
		invokes the iraf method imhead on the given image's filename
		parses the result to be stored in image instance given
		
		parameter image - 
			the image on which to invoke imexam
		'''
		
		imhead_return = iraf.imhead(image.filename, Stdout=1)[0].strip()
		imhead_info = imhead_return.split("/")[-1]
		# "VELA01_a0.110_0006317__skipir_CAMERA0-BROADBAND_F125W_simulation.fits[600,600][real]:" 
		
		
		frame_dimensions = imhead_info.split("[")[1].replace("]", "").split(",")
		# ["600","600"]
	
		image.width = frame_dimensions[0]
		# "600"
		
		image.height = frame_dimensions[1]
		# "600"
		
		image.galaxyID = imhead_info.split("_")[0]
		#"VELA01"
		
		image.filter = imhead_info.split("_")[6]
		#"F125W"	
		
		image.timeStep = imhead_info.split("_")[1].split(".")[1]
		# "110"
		
		cam_str = imhead_info.split("_")[5].split("-")[0]
		# "CAMERA0"
		
		# these statements accomodate for camera numbers
		image.camera = ""
		for c in cam_str:
			if c.isDigit():
				image.camera = image.camera + c
	
	
	def run_imexam(image):
		'''
		method calls iraf's imexam method on the area defined by the four
		
		parameter image - 
			the image on which to invoke imexam
		'''
	
		# writes the center coordinates to a file for imexam
		coordsFilename = "coords"
		
		# this is t prevent parallel processes from overwriting the coords.tmp file
		while os.path.isfile(coordsFilename + ".tmp"):
			coordsFilename = coordsFilename + "x"
		coordsFilename = coordsFilename + ".tmp"
		
		# TODO: should we run minmax to get them if the user doesnt use sextractor?
		if not image.model.centerCoords[0] or not image.model.centerCoords[1]:
			print("need center coordinates to run imexam., exiting")
			exit()
			
		# create coords[x*].tmp for writing max coordinate to pass to imexam
		os.system('touch '+ coordsFilename)
		write_coords = open(coordsFilename, 'w')
		write_coords.write(	image.model.centerCoords[0] + " " + 
							image.model.centerCoords[1])
		write_coords.close()	
		
		# run imexam passing coords.tmp, 
		# data is returned in the last two elements of the return array of strings
		imexam_out = iraf.imexam(imageFilename, use_display=0, 
								imagecur=coordsFilename, Stdout=1)[-2:]
	
		# ['COL	LINE X Y', 
		# 'R MAG FLUX SKY PEAK E PA BETA ENCLOSED MOFFAT DIRECT']
		
		# delete coords.tmp, no longer needed
		os.system('rm ' + coordsFilename) 
		
		xy = imexam_out[0].strip().split()
		data = imexam_out[1].strip().split()
		
		image.model.centerCoords[0] = float(xy[0])	
		image.model.centerCoords[1] = float(xy[1])
		image.model.radius = float(data[0])						
		image.model.magnitude = float(data[1])						
		
		# these might be indef, if so then set to default values
		try:
			image.model.ba = 1.0 - float(data[5])				# b/a=1-e
		except:
			print("using default values for b/a")
			image.model.ba = 1.0
		
		#this gives the position angle. iraf measures up from x. Galfit down from y
		try:
			image.model.angle = float(data[6]) - 90.0				
		except:
			print("using default values for position angle")
			image.model.angle = 0.0
		
		
	def defineGalfitFilenames(image):
		directory_location = "/".join(image.filename.split("/")[:-1]) + "/"
		filename = (directory_location + image.galaxyID + "_" + 
					image.timeStep + '_cam' + str(image.camera) + 
					'_' + image.filter)
		self.galfit_single_parameter_filename =	filename + '_single_param.txt'
		self.galfit_single_output_filename =	filename + "_single_multi.fits"
		self.galfit_single_result_filename =	filename + "_single_result.txt"
		self.galfit_bulge_parameter_filename =	filename + '_bulge_param.txt'
		self.galfit_bulge_output_filename =		filename + "_bulge_multi.fits"
		self.galfit_bulge_result_filename =		filename + "_bulge_result.txt"

	
	def write_galfit_single_parameter(image):
		'''
		writes the galfit parameter file for the single component
		
		parameter image - 
			the image on which galfit is being run
		'''
		os.system('touch ' + self.galfit_single_parameter_filename)
		galfitSingleParamFile = open(self.galfit_single_parameter_filename,'w')
		galfitSingleParamFile.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
		galfitSingleParamFile.write("A) " + image.filename + 
			"						#Input data image block\n")
		galfitSingleParamFile.write("B) " + self.galfit_single_output_filename +
			"						#Output data image block\n")
		galfitSingleParamFile.write("C)" + " none" + 
			"						#Sigma image name (made from data if blank or 'none')\n")
		galfitSingleParamFile.write("D) " + self.psf + 
			"						#Input PSF image and (optional) diffusion kernel\n")
		galfitSingleParamFile.write("E)" + " 1" + 
			"						#PSF fine sampling factor relative to data\n")
		galfitSingleParamFile.write("F) " + self.segmentationMapFilename +							
			"						#Bad pixel mask (FITS file or ASCIIcoord list)\n")
		galfitSingleParamFile.write("G)" + " none" + 
			"						#File with parameter constraints (ASCII file)\n")
		galfitSingleParamFile.write("H) " + "1" + " " + str(image.width) + " " + 
											"1" + " " + str(image.height) + 
			"						#Image region to fit (xmin xmax ymin ymax)\n")
		galfitSingleParamFile.write("I)" + " 200 200" + 
			"						#Size of the convolution box (x y)\n")
		galfitSingleParamFile.write("J) " + str(mpZeropoint) + 
			"						#Magnitude photometric zeropoint\n")
		galfitSingleParamFile.write("K) " + str(self.plateScale) + "	 " + str(self.plateScale) + 
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
		galfitSingleParamFile.write(" 1) " + str(image.model.centerCoords[0]) + "	" + str(image.model.centerCoords[1]) + "	1	1			#Position x,y\n")
		galfitSingleParamFile.write(" 3) " + str(image.model.magnitude) + "	1			#Integrated Magnitude\n")
		galfitSingleParamFile.write(" 4) " + str(image.model.radius) + "			1			#R_e (half-light radius)	[pix]\n")
		galfitSingleParamFile.write(" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
		galfitSingleParamFile.write(" 6) 0.0000		0			#	-----\n")
		galfitSingleParamFile.write(" 7) 0.0000		0			#	-----\n")
		galfitSingleParamFile.write(" 8) 0.0000		0			#	-----\n")
		galfitSingleParamFile.write(" 9) " + str(image.model.ba) + "			1			#Axis ratio (b/a)\n")
		galfitSingleParamFile.write(" 10) " + str(image.model.angle) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
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
	
	
	def write_galfit_bulge_parameter():
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
		singleResultFile = open(self.galfit_single_result_filename, "r")		
		resultContents = singleResultFile.readlines()
		singleResultFile.close()
		
		# gather info about first fit while modifying output and coonstraint parameters
		bulgeParamStr = ""		
		positionLine = ""
		magnitudeLine = ""
		radiusLine = ""
		for resultLine in resultContents:
		
			if resultLine.strip()[:2] == "B)":
				bulgeParamStr = (bulgeParamStr + "B) " + self.galfit_bulge_output_filename +
							"						#Output data image block\n")
				
			elif resultLine.strip()[:2] == "G)":
				bulgeParamStr = (bulgeParamStr + "G) " + self.galfit_constraint_filename + 
							"						#File with parameter constraints (ASCII file)\n")
			
			elif not positionLine and resultLine.strip()[:2] == "1)":
				positionLine = resultLine.strip().split(" ")
				bulgeParamStr = bulgeParamStr + resultLine
				
			elif not magnitudeLine and resultLine.strip()[:2] == "3)":
				magnitudeLine = resultLine.strip().split(" ")
				bulgeParamStr = bulgeParamStr + resultLine
				
			elif not radiusLine and resultLine.strip()[:2] == "4)":
				radiusLine = resultLine.strip().split(" ")
				bulgeParamStr = bulgeParamStr + resultLine
				
			else:
				bulgeParamStr = bulgeParamStr + resultLine
				
		# store results of single component run into variables
		resultX = positionLine[1]
		resultY = positionLine[2]
		resultMag = magnitudeLine[1]
		resultRad = radiusLine[1]
		
		
		# write the third component to the end of galfit_single_result_filename
		bulgeParamStr = bulgeParamStr + ("# Componenet number: 3\n")
		bulgeParamStr = bulgeParamStr + (" 0) sersic					#Component type\n")
		bulgeParamStr = bulgeParamStr + (" 1) " + resultX + " " + resultY + " 1	1			#Position x,y\n")
		bulgeParamStr = bulgeParamStr + (" 3) " + resultMag + "	1			#Integrated Magnitude\n")
		bulgeParamStr = bulgeParamStr + (" 4) " + resultRad + "			1			#R_e (half-light radius)	[pix]\n")
		bulgeParamStr = bulgeParamStr + (" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
		bulgeParamStr = bulgeParamStr + (" 6) 0.0000		0			#	-----\n")
		bulgeParamStr = bulgeParamStr + (" 7) 0.0000		0			#	-----\n")
		bulgeParamStr = bulgeParamStr + (" 8) 0.0000		0			#	-----\n")
		bulgeParamStr = bulgeParamStr + (" 9) 1" + "			1			#Axis ratio (b/a)\n")
		bulgeParamStr = bulgeParamStr + (" 10) 0" + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
		bulgeParamStr = bulgeParamStr + (" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
		bulgeParamStr = bulgeParamStr + ("\n")
	
		# write the bulge parameter file using the modified contents of the single results
		galfitBulgeParamFile = open(self.galfit_bulge_parameter_filename, "w")
		galfitBulgeParamFile.write(bulgeParamStr)
		galfitBulgeParamFile.close()
		
		
	def run_galfit(image):
		'''
		opens file (parameter) containing list of images and loops over every image, 
		running galfit and writing the results to the same directory as the images
		
		calls methods to invoke iraf methods to write the initial galfit param file
		for each image
		
		parameter imageFilename -
			the image on which galfit will be run
		'''
		
		# run iraf's imhead method to get image information
		self.run_imhead(image)
		
		# calls iraf's imexam method, passing filename as a parameter
		# along with center coordinates and two points defining area
		# stores initial estimates in image model
		self.run_imexam(image)
		
		# write to the log if default values are being used
		if float(image.model.ba) == 1.0:
			self.logMsg = self.logMsg + "Default b/a used. "
		if float(image.model.angle) == 0.0:
			self.logMsg = self.logMsg + "Default position angle used. "
		
		# define filenames
		self.defineGalfitFilenames()
		
		# writes the single component parameter file, given filenames and galxy parameters
		self.write_galfit_single_parameter(image)
	
		# run galfit on paramter file
		os.system('galfit ' + self.galfit_single_parameter_filename)
	
		# detects atomic galfit error
		# If not failure then removes temp files, otherwise returns
		if os.path.isfile(self.galfit_single_output_filename):
			# remove temp files
			os.system("rm " + "fit.log")
			os.system("mv galfit.01 " + self.galfit_single_result_filename)
		else:
			self.logMsg = self.logMsg + ("galfit failed on single component, " + 
								"probably mushroom (atomic galfit error)")
		
		# done unless command line specified that a second galfit run
		# should be done by adding a bulge component to the result of the first run
		if self.includeBulgeComponent:
		
			# reads the results of the first run and returns it as a long string
			# with the output and constraint modified for bulge run and the
			# new bulge component appended with some intitial guess parameters
			self.write_galfit_bulge_parameter()

			# run galfit on paramter file
			os.system('galfit ' + self.galfit_bulge_parameter_filename)
		
			# detects atomic galfit error
			# If not failure then removes temp files, otherwise returns
			if os.path.isfile(self.galfit_bulge_output_filename):
				# remove temp files
				os.system("rm " + "fit.log")
				os.system("mv galfit.01 " + self.galfit_bulge_result_filename)
				#os.system('mv ' + galfit_output_filename + ' ' + galfit_bulge_output_filename)
			else:
				self.logMsg = self.logMsg + ("galfit failed on bulge component, " + 
									"probably mushroom (atomic galfit error)")
		
		# if we get here then nothing went wrong!
		self.logMsg = self.logMsg + "success"
		

	def modelImages(self, imageFilenames):
		'''
		takes the list of images the command line and begins analysis
		'''	
		
		# verify there are images to model 
		if len(imageFilenames) == 0:
			print("no images in given file (positional argument). Exiting.")
			exit()
			
		# convenience boolean for dealing with conditionally running sextractor
		runSextractor = self.realSextractor or self.simSextractor
			
		# sextractor specific tasks that can be done outside of image loop
		if runSextractor:
		
			# write the sextractor config file, which will be used for all images
			self.write_sextractor_config_file()
	
		# set the log header
		self.logMsg = "run on " + time.strftime("%m-%d-%Y") + "\n"
		
		# this loops through every image in images file and removes whitespace
		imageFilenames = [imageFilename.strip() for imageFilename in imageFilenames]
		
		# this loops through every image, runs galfit, and writes the log string 
		for imageFilename in imageFilenames:
			
			curImage = SimImage(imageFilename)
			
			# run galfit, preventing crashes but printing and logging errors in log
			try:
						
				# run sextractor (if command line set to do so) for each image
				# galfit uses the resulting segmentation map
				if runSextractor:
						
					# run sextractor, which updates image and seg map
					self.run_sextractor(curImage)
					
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
			
				# galfit returns a string indicating success or some failure
				run_galfit(image)
			
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
				self.logMsg = self.logMsg + errorMsg
			
			# every image on its own line in the log file
			self.logMsg = self.logMsg + "\n"
		
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
		log.write(self.logMsg)
		log.close()
		
		print ("Done!\nIn order to summarize results run sumgalfit.py")
	
	
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
	else:		
		# read list of image filenames from input file
		inputFile = open(args[0], 'r')
		imageFilenames = inputFile.readlines()
		inputFile.close()
		
		# holds methods for analyzing simulations
		sp = SimulationProcessor()
		
		# parse the command line options
		sp.parseGalfitOptions(parser, options)
		sp.parseSextractorOptions(parser, options.simSextractor, 
									options.realSextractor, args[1:])
		
		# pass list of images to method that orchestrates modeling
		sp.modelImages(imageFilenames)
		
	
		
		
		
		
		
		
		
		
		
	