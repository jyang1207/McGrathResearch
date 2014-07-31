#!/usr/bin/env python

import os
import sys
from pyraf import iraf
import pyfits
import multiprocessing
import time
import math
from optparse import OptionParser
		
class ModelGenerator:
	'''
	The controller class, which holds methods for analyzing simulations
	'''
		
	def __init__(self, destDirectory=""):
		'''constructor sets some initial field defaults'''
		self.outputCatFilename = "generic.cat" # -CATALOG_NAME <filename>
		self.segmentationMapFilename = "check.fits" # -CHECKIMAGE_NAME <filename>
		self.sextractorOptionsList = []
		# "/".join(image.filename.split("/")[:-1]) + "/"
		
		# verify and set the destination directory
		if destDirectory and not destDirectory.endswith("/"):
			destDirectory = destDirectory + "/"
		self.destDirectory = destDirectory
		if self.destDirectory and not os.path.isdir(self.destDirectory):
			os.mkdir(self.destDirectory)
			
		# set the names of the sextractor configuration files
		self.sextractorConfigFilename = destDirectory + "configInit.sex"
		self.sextractorReduceComponentConfigFilename = destDirectory + "configFewerComp.sex"
		
		
	def parseGalfitOptions(self, parser, options):
		'''
		parses the galfit options from the command line and uses to
		initialize fields of this class instance with options or defaults
		'''

		# store command line options into fields of this class instance
		self.mpZeropoint = options.mpz
		self.plateScale = options.plate
		self.includeBulgeComponent = options.bulge
		self.galfitOff = options.galfitOff

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
		
		
	def parseSextractorOptions(self, parser, realSextractor, sextractorOptions):
		'''
		parses the sextractor options from the command line and uses to
		initialize fields with options or defaults
		'''
		
		# all positional arguments after image list filename must come in pairs
		if len(sextractorOptions) % 2 != 0:
			parser.error("Sextractor options parsing error, odd number of options. " + 
				"Expecting sextractor options to be <keyword> <value> pairs")
				
		# store running sextractor booleans
		self.realSextractor = realSextractor
		
		# if running sextractor, gather options in a list and update fields
		self.sextractorOptionsList = []
		for sexIndex, sexOpt in enumerate(sextractorOptions):
			if sexIndex % 2 == 0:	# even, consider a keyword
				if sexOpt.upper() == "CATALOG_NAME":
					continue#self.outputCatFilename = sextractorOptions[sexIndex + 1]
				elif sexOpt.upper() == "CHECKIMAGE_NAME":
					continue#self.segmentationMapFilename = sextractorOptions[sexIndex + 1]
				self.sextractorOptionsList.append("-" + sexOpt.upper())
			else:				# even, consider the value for previous keyword
				self.sextractorOptionsList.append(sexOpt)
			
		
	def write_sextractor_config_file(self, fname):
		'''
		writes the sextractor configuration file
		'''
		
		requiredFilesDirectory = "requiredFiles/"
		# variables describing the sextractor config file
		catalogName = "generic.cat" # Name of the output catalogue. If the name "STDOUT" is given and CATALOG TYPE is set to ASCII, ASCII HEAD, ASCII SKYCAT, or ASCII VOTABLE the catalogue will be piped to the standard output (stdout
		catalogType = "ASCII_HEAD"	# Format of output catalog:
									# ASCII - ASCII table; the simplest, but space and time consuming,
									# ASCII HEAD - as ASCII, preceded by a header containing information about the content,
									# ASCII SKYCAT - SkyCat ASCII format (WCS coordinates required)
									# ASCII VOTABLE - XML-VOTable format, together with meta-data,
									# FITS 1.0 - FITS format as in SExtractor 1
									# FITS LDAC - FITS "LDAC" format (the original image header is copied).
		sextractorParamFilename = requiredFilesDirectory + "sex.param"
		#------------------------------- Extraction ----------------------------------
		detectType = "CCD" 
		#fitsUnsigned = Force 16-bit FITS input data to be interpreted as unsigned integers.
		#flagImage = File name(s) of the flagimage(s)
		flagType = "OR" # Combination method for flags on the same object:
						# OR arithmetical OR,
						# AND arithmetical AND,
						# MIN minimum of all flag values,
						# MAX maximum of all flag values,
						# MOST most common flag value.
		
		# these variable values change for simulation template vs real template
		if self.realSextractor:
			detectMinArea = "5" # Minimum number of pixels above threshold triggering detection
			detectThreshold = "0.75"	#Detection threshold (0-2). 1 argument: (ADUs or relative to Background RMS, see THRESH TYPE). 2 arguments: R (mag.arcsec 2 ), Zero-point (mag).
			deblendThreshold = "16" # Minimum contrast parameter for deblending.
			deblendMinContrast = "0.0001"
			cleanBool = "Y" # If true, a cleaning of the catalog is done before being written to disk.
			cleanParam = "1.0"	# Efficiency of cleaning.
			analysisThreshold = "5" #Threshold (in surface brightness) at which CLASS STAR and FWHM operate. 1 argument: relative to Background RMS. 2 arguments: mu (mag/arcsec 2 ), Zero-point (mag)
		elif fname == self.sextractorConfigFilename:
			detectMinArea = "10000" # Minimum number of pixels above threshold triggering detection
			detectThreshold = "20"	#Detection threshold (0-2). 1 argument: (ADUs or relative to Background RMS, see THRESH TYPE). 2 arguments: R (mag.arcsec 2 ), Zero-point (mag).
			deblendThreshold = "16" # Minimum contrast parameter for deblending.
			deblendMinContrast = "0.005"	# Minimum contrast parameter for deblending.
			cleanBool = "Y" # If true, a cleaning of the catalog is done before being written to disk.
			cleanParam = "0.5"	# Efficiency of cleaning.
			analysisThreshold = "5" #Threshold (in surface brightness) at which CLASS STAR and FWHM operate. 1 argument: relative to Background RMS. 2 arguments: mu (mag/arcsec 2 ), Zero-point (mag).
		else:# fname == self.sextractorReduceComponentConfigFilename:
			detectMinArea = "10000" # Minimum number of pixels above threshold triggering detection
			detectThreshold = "50"	#Detection threshold (0-2). 1 argument: (ADUs or relative to Background RMS, see THRESH TYPE). 2 arguments: R (mag.arcsec 2 ), Zero-point (mag).
			deblendThreshold = "16" # Minimum contrast parameter for deblending.
			deblendMinContrast = "0.02"	# Minimum contrast parameter for deblending.
			cleanBool = "Y" # If true, a cleaning of the catalog is done before being written to disk.
			cleanParam = "0.5"	# Efficiency of cleaning.
			analysisThreshold = "20" #Threshold (in surface brightness) at which CLASS STAR and FWHM operate. 1 argument: relative to Background RMS. 2 arguments: mu (mag/arcsec 2 ), Zero-point (mag).
			
		filterBool = "Y"	#  If true,filtering is applied to the data before extraction.
		filterName = requiredFilesDirectory + "sex.conv" # Name and path of the file containing the filter definition
		#filterThresh = ?	# Lower and higher thresholds (in back-ground standard deviations) for a pix-el
							#to be consideredin filtering (used for retinafiltering only).
		#threshType = ? # Meaning of the DETECT_THRESH and ANALYSIS_THRESH parameters:
							# RELATIVE - scaling factor to the background RMS,
							# ABSOLUTE - absolute level (in ADUs or in surface brightness)
			
		maskType = "CORRECT"	#CORRECT - replace by values of pixels symmetric with respect to the source center.
								#BLANK --put detected pixels belonging to neighbors to zero,
								#NONE no masking,
		#------------------------------ Photometry -----------------------------------
		photoApertureDiameter = "10."	# these threes variables are related to the Kron radius, which is introduced as a accurate flexible aperture
										# that would capture most of the flux from an object.
		magAutoKronFact = "2.5"
		magAutoMinRadius = "3.5"
		saturLevel = "120"	#Pixel value above which it is considered saturated.
		magGamma = "4.0"	#gamma of the emulsion (slope of the response function). Takes effect in PHOTO mode 
							#only but NEEDS to be specified, even for CCD images.
		pixelScale = "0.06" #Pixel size in arcsec.
		#------------------------- Star/Galaxy Separation ----------------------------
		seeingFWHM = "0.18"# ????? Seeing_FWHM FWHM of stellar images in arcsec.This quantity is used only
							#for the neural network star/galaxy separation as expressed in the CLASS STAR output.
		starNNWFilename = requiredFilesDirectory + "default.nnw" #Name of the file containing the neural network weights for star/galaxy separation.
		#------------------------------ Background -----------------------------------
		backSize = "600"	# Size, or Width, Height (inpixels) of a background mesh.
		backFilterSize = "9"	#Size, or Width, Height (inbackground meshes) of the background-filtering mask.
		backPhotoType = "LOCAL" #Background used to compute magnitudes:
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
		weightType = "NONE"		#variance-map derived from an external weight-map,
								#MAP_VAR external variance-map,
								#MAP RMS variance-map derived from an external RMS-map,
								#BACKGROUND variance-map derived from the image itself,
								#NONE no weighting,
		#WEIGHT IMAGE --File name of the detection and measurement weight image , respectively.
		#WEIGHT GAIN -- If true, weight maps are considered as gain maps.
		magZeropoint = "25.96"	   #Zero-point offset to be applied to magnitudes.
		#------------------------- End of variable definitions -----------------------
		
		
		# use above variables to write the config file
		os.system('touch ' + fname)
		sextractorConfigFile = open(fname,'w')
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
				"CATALOG_NAME	 " + catalogName + "\n")
		
		if catalogType:
			# ASCII_HEAD
			sextractorConfigFile.write(
				"CATALOG_TYPE	 " + catalogType + 
				"	   # \"NONE\",\"ASCII_HEAD\",\"ASCII\",\"FITS_1.0\", or \"FITS_LDAC\"\n")
		
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
				"DETECT_TYPE	 " + detectType + 
				"			  # \"CCD\" or \"PHOTO\" (*)\n")
		
		if flagType:
			#OR
			sextractorConfigFile.write(
				"FLAG_TYPE		 " + flagType + "\n")
		
		if detectMinArea:
			#5
			sextractorConfigFile.write(
				"DETECT_MINAREA	 " + detectMinArea +
				"				# minimum number of pixels above threshold\n")
		
		if detectThreshold:
			#0.75
			sextractorConfigFile.write(
				"DETECT_THRESH	 " + detectThreshold +
				"			 # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n")
		
		if analysisThreshold:
			#5
			sextractorConfigFile.write(
				"ANALYSIS_THRESH " + analysisThreshold +
				"				# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n")
		
		if filterBool:
			#Y
			sextractorConfigFile.write(
				"FILTER			 " + filterBool +
				"				# apply filter for detection (\"Y\" or \"N\")?\n")
		
		if filterName:
			#tophat_9.0_9x9.conv
			sextractorConfigFile.write(
				"FILTER_NAME	 " + filterName + 
				"	  # name of the file containing the filter\n")
		
		if deblendThreshold:
			#16
			sextractorConfigFile.write(
				"DEBLEND_NTHRESH " + deblendThreshold + 
				"			   # Number of deblending sub-thresholds\n")
		
		if deblendMinContrast:
			#0.0001
			sextractorConfigFile.write(
				"DEBLEND_MINCONT " + deblendMinContrast +
				"  # Minimum contrast parameter for deblending\n\n")
		
		if cleanBool:
			#Y
			sextractorConfigFile.write(
				"CLEAN			 " + cleanBool +
				"				# Clean spurious detections? (Y or N)?\n")
		
		if cleanParam:
			#1.0
			sextractorConfigFile.write(
				"CLEAN_PARAM	 " + cleanParam + 
				"			  # Cleaning efficiency\n\n")
		
		if maskType:
			#CORRECT
			sextractorConfigFile.write(
				"MASK_TYPE		 " + maskType + 
				"		  # type of detection MASKing: can be one of\n" +
				"								 # \"NONE\", \"BLANK\" or \"CORRECT\"\n")
		
		sextractorConfigFile.write(
'''

#------------------------------ Photometry -----------------------------------

''')
		
		if photoApertureDiameter:
			#10.
			sextractorConfigFile.write(
				"PHOT_APERTURES	  " + photoApertureDiameter +
				"	# MAG_APER aperture diameter(s) in pixels\n")
		
		if magAutoKronFact and magAutoMinRadius:
			#2.5, 3.5
			sextractorConfigFile.write(
				"PHOT_AUTOPARAMS " + magAutoKronFact + ", " + magAutoMinRadius +
				"		 # MAG_AUTO parameters: <Kron_fact>,<min_radius>\n")
		
		if saturLevel:
			#120
			sextractorConfigFile.write(
				"SATUR_LEVEL	 " + saturLevel + 
				"			  # level (in ADUs) at which arises saturation\n")
		
		if magGamma:
			#4.0
			sextractorConfigFile.write(
				"MAG_GAMMA		 " + magGamma +
				"			  # gamma of emulsion (for photographic scans)\n")
		
		if pixelScale:
			#0.06
			sextractorConfigFile.write(
				"PIXEL_SCALE	 " + pixelScale +
				"			 # size of pixel in arcsec (0=use FITS WCS info).\n")
		
		sextractorConfigFile.write(
'''

#------------------------- Star/Galaxy Separation ----------------------------

''')
		
		if seeingFWHM:
			#0.18
			sextractorConfigFile.write(
				"SEEING_FWHM	 " + seeingFWHM + 
				"			 # stellar FWHM in arcsec\n")
		
		if starNNWFilename:
			#default.nnw
			sextractorConfigFile.write(
				"STARNNW_NAME	 " + starNNWFilename + 
				"	  # Neural-Network_Weight table filename\n")
		
		sextractorConfigFile.write(
'''

#------------------------------ Background -----------------------------------

''')
		
		if backSize:
			#256
			sextractorConfigFile.write(
				"BACK_SIZE		 " + backSize + 
				"			  # Background mesh: <size> or <width>,<height>\n")
		
		if backType:
			#not in original
			sextractorConfigFile.write(
				"BACK_TYPE		 " + backType + 
				"			  # What background is subtracted from the images:" +
'''
				# AUTO - The internal interpolated background-map. In the manual it says "INTERNAL" here but the keyword is AUTO. 
				# MANUAL A user-supplied constant value provided in BACK VALUE.")
''')
		
		if backValue:
			#not in original
			sextractorConfigFile.write(
				"BACK_VALUE		  " + backValue + 
				"			  # in BACK TYPE MANUAL mode, the constant value to be subtracted from the images.\n")
		
		if backFilterSize:
			#9
			sextractorConfigFile.write(
				"BACK_FILTERSIZE " + backFilterSize +
				"				# Background filter: <size> or <width>,<height>\n\n")
		
		if backPhotoType:
			#LOCAL
			sextractorConfigFile.write(
				"BACKPHOTO_TYPE	 " + backPhotoType + 
				"			# can be \"GLOBAL\" or \"LOCAL\" (*)\n")
		
		if backPhotoThickness:
			#100
			sextractorConfigFile.write(
				"BACKPHOTO_THICK " + backPhotoThickness + 
				"			  # thickness of the background LOCAL annulus (*)\n")
		
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
				"	# can be one of \"NONE\", \"BACKGROUND\"," + 
'''
				# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
				# "-OBJECTS", "SEGMENTATION", "APERTURES",
				# or "FILTERED" (*)
''')
		
		sextractorConfigFile.write(
'''

#--------------------- Memory (change with caution!) -------------------------

# MEMORY_OBJSTACK	50000		# number of objects in stack
# MEMORY_PIXSTACK	1000000		# number of pixels in stack
# MEMORY_BUFSIZE	8500		# number of lines in buffer

''')
		
		if memoryObjStack:
			#2000
			sextractorConfigFile.write(
				"MEMORY_OBJSTACK " + memoryObjStack +
				"			 # number of objects in stack\n")
		
		if memoryPixStack:
			#200000
			sextractorConfigFile.write(
				"MEMORY_PIXSTACK " + memoryPixStack + 
				"		   # number of pixels in stack\n")
		
		if memoryBufferSize:
			#2048
			sextractorConfigFile.write(
				"MEMORY_BUFSIZE	 " + memoryBufferSize + 
				"			 # number of lines in buffer\n")
		
		sextractorConfigFile.write(
'''
 
#----------------------------- Miscellaneous ---------------------------------

''')
		
		if verboseType:
			#QUIET
			sextractorConfigFile.write(
				"VERBOSE_TYPE	 " + verboseType + 
				"			# can be \"QUIET\", \"NORMAL\" or \"FULL\" (*)\n")
		
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
# WEIGHT_TYPE	  MAP_RMS,MAP_RMS
# WEIGHT_THRESH	  10000.0


''')
		
		if magZeropoint:
			#25.96
			sextractorConfigFile.write(
				"MAG_ZEROPOINT	 " + magZeropoint + 
				"			# H-band magnitude zero-point\n")
	
	
	def run_sextractor(self, image, configFilename):
		'''
		runs sextractor on the given image, producing
		*.cat and *.fits in the destination directory, which 
		are galaxy info and the segmentation map, respectively.
		
		Updates the model of the SimImage instance passed by parameter
		Updates segmentation map to remove the galaxy closest to (300, 300)
		
		parameter image -
			the image on which sextractor will be run
			
		returns - boolean indicating if a necessary parameter is missing
		'''
	
		# run sextractor
		# sex <image> [<image2>][-c <configuration_file>][-<keyword> <value>]
		self.outputCatFilename = (self.destDirectory + 
									".".join(image["filename"].split("/")[-1
										].split(".")[:-1]) + ".cat")
		self.segmentationMapFilename = (self.destDirectory + 
									".".join(image["filename"].split("/")[-1
										].split(".")[:-1]) + "_check.fits")
		os.system(	"sex " + image["filename"] + 
					" -c " + configFilename + 
					" -CATALOG_NAME " + self.outputCatFilename + 
					" -CHECKIMAGE_NAME " + self.segmentationMapFilename + 
					" " + " ".join(self.sextractorOptionsList))
		
		# gather galaxy info from .cat file into image["models"]
		with open(self.outputCatFilename, 'r') as outputCatFile:
			outputCatContents = outputCatFile.readlines()
			
		# count the number of components in the catalog output
		numGalaxies = 0
		for catalogOutputLine in reversed(outputCatContents):
			if catalogOutputLine[0] == "#":
				break
			else:
				numGalaxies = numGalaxies + 1
		
		# if there are too many galaxies after the first fit, recursively
		# run sextractor again, now with a different configuration file
		if (configFilename == self.sextractorConfigFilename) and (numGalaxies > 3):
			return self.run_sextractor(image, self.sextractorReduceComponentConfigFilename)
		
		# loop over every line in catalog to build list of image models				
		indexDict = {}
		errorStr = ""
		for catalogOutputLine in outputCatContents:
			defaultStr = ""
			galaxyOutputList = catalogOutputLine.strip().split()
			
			# if in the top comment section
			if galaxyOutputList[0] == "#":
				indexDict[galaxyOutputList[2].upper()] = int(galaxyOutputList[1]) - 1
			
			# if in the bottom components section
			else:
				
				# gather galaxy information
				if "X_IMAGE" in indexDict:
					galaxyX = float(galaxyOutputList[indexDict["X_IMAGE"]])
				else:
					errorStr = (errorStr + 
						"X_IMAGE is a required field of the parameter file\n")
				
				if "Y_IMAGE" in indexDict:
					galaxyY = float(galaxyOutputList[indexDict["Y_IMAGE"]])
				else:
					errorStr = (errorStr + 
						"Y_IMAGE is a required field of the parameter file\n")
				
				if errorStr:
					print(errorStr)
					return False

				if "MAG_AUTO" in indexDict:
					galaxyMag = float(galaxyOutputList[indexDict["MAG_AUTO"]])
				else:
					defaultStr = defaultStr + "MAG_AUTO "
					galaxyMag = 22.3
				
				if "FLUX_RADIUS" in indexDict:
					galaxyRad = float(galaxyOutputList[indexDict["FLUX_RADIUS"]])
					'''if "KRON_RADIUS" in indexDict:
						galaxySers = galaxyRad/float(galaxyOutputList[indexDict["KRON_RADIUS"]])
					else:
						defaultStr = defaultStr + "KRON_RADIUS(for sersic index = FLUX_RADIUS/KRON_RADIUS) "
						galaxySers = 2.5'''
				else:
					defaultStr = defaultStr + "FLUX_RADIUS "
					galaxyRad = 50.0
				galaxySers = 2.5 # need to have accurate KronRadius for above
								
				if "A_IMAGE" in indexDict:
					galaxyA = float(galaxyOutputList[indexDict["A_IMAGE"]])
				else:
					defaultStr = defaultStr + "A_IMAGE "
					galaxyA = 0.0	
					
				if "B_IMAGE" in indexDict:
					galaxyB = float(galaxyOutputList[indexDict["B_IMAGE"]])
				else:
					defaultStr = defaultStr + "B_IMAGE "
					galaxyB = 0.0	
					
				if "ELLIPTICITY" in indexDict:
					galaxyE = float(galaxyOutputList[indexDict["ELLIPTICITY"]])
				else:
					defaultStr = defaultStr + "ELLIPTICITY "
					galaxyE = 1.0
					
				if "THETA_IMAGE" in indexDict:
					galaxyAng = float(galaxyOutputList[indexDict["THETA_IMAGE"]])
				else:
					defaultStr = defaultStr + "THETA_IMAGE "
					galaxyAng = 0.0
					
				if defaultStr:
					print(defaultStr + "are missing from the sextractor" + 
								" parameter file, default(s) are being used")
					
				# if eith a or b are zero, use 1-ellipticity, otherwise use b/a
				if galaxyA and galaxyB:
					galaxyBA = galaxyB/galaxyA
				else:
					galaxyBA = 1.0-galaxyE

				# store (x,y) of the galaxy as the image's model
				image["models"].append({"centerCoords":[galaxyX,galaxyY],
										"magnitude":galaxyMag,
										"radius":galaxyRad,
										"ba":galaxyBA,
										"angle":galaxyAng,
										"sers":galaxySers})
				
		self.logMsg = (self.logMsg + " configFilename = " + configFilename + 
						" numGalaxies = " + str(numGalaxies))
		return True
				

	def write_ds9_reg_file(self, image):
		'''
		use the image to write a reg file for ds9, named according to the
		image filename and with a circle for each model in the image's list
		of models
		'''
		# join method is to handle the relative path names
		regFileName = (self.destDirectory + 
						".".join(image["filename"].split("/")[-1
							].split(".")[:-1]) + ".reg")
		
		# write the header
		regFileString = ("#Region file format: DS9 version 4.1\n" +
						"global color=green dashlist=8 3 width=1 " + 
						"font=\"helvetica 10 normal roman\" select=1 " + 
						"highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 " +
						"include=1 source=1\n" + "physical\n")
		
		# write a circle for each component of the image
		for model in image["models"]:
			regFileString = (regFileString + "circle(" + 
								str(model["centerCoords"][0]) + "," + 
								str(model["centerCoords"][1]) + "," + 
								str(model["radius"]) + ")\n")
		
		# write file
		with open(regFileName, "w") as regFile:
			regFile.write(regFileString)
			
			
	def run_imhead(self, image):
		'''
		invokes the iraf method imhead on the given image's filename
		parses the result to be stored in image instance given
		
		parameter image - 
			the image on which to invoke imexam
		'''
		
		imhead_return = iraf.imhead(image["filename"], Stdout=1)[0].strip()
		imhead_info = imhead_return.split("/")[-1]
		# "VELA01_a0.110_0006317__skipir_CAMERA0-BROADBAND_F125W_simulation.fits[600,600][real]:" 
		
		
		frame_dimensions = imhead_info.split("[")[1].replace("]", "").split(",")
		# ["600","600"]
	
		image["width"] = float(frame_dimensions[0])
		# "600"
		
		image["height"] = float(frame_dimensions[1])
		# "600"
		
		image["galaxyID"] = imhead_info.split("_")[0]
		#"VELA01"
		
		image["filter"] = imhead_info.split("_")[6]
		#"F125W"	
		
		image["timeStep"] = imhead_info.split("_")[1].split(".")[1]
		# "110"
		
		cam_str = imhead_info.split("_")[5].split("-")[0]
		# "CAMERA0"
		
		# these statements accomodate for camera numbers
		image["camera"] = ""
		for c in cam_str:
			if c.isdigit():
				image["camera"] = image["camera"] + c
		
		
	def defineGalfitFilenames(self, image):
		filename = (self.destDirectory + image["galaxyID"] + "_" + 
					image["timeStep"] + '_cam' + str(image["camera"]) + 
					'_' + image["filter"])
		self.galfit_single_constraint_filename =filename + '_single_constraint.txt'
		self.galfit_single_parameter_filename = filename + '_single_param.txt'
		self.galfit_single_output_filename =	filename + "_single_multi.fits"
		self.galfit_single_result_filename =	filename + "_single_result.txt"
		self.galfit_bulge_constraint_filename = filename + '_bulge_constraint.txt'
		self.galfit_bulge_parameter_filename =	filename + '_bulge_param.txt'
		self.galfit_bulge_output_filename =		filename + "_bulge_multi.fits"
		self.galfit_bulge_result_filename =		filename + "_bulge_result.txt"

	
	def write_galfit_single_parameter(self, image):
		'''
		writes the galfit parameter file for the single component
		
		parameter image - 
			the image on which galfit is being run
		'''
		os.system('touch ' + self.galfit_single_parameter_filename)
		galfitSingleParamFile = open(self.galfit_single_parameter_filename,'w')
		galfitSingleParamFile.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
		galfitSingleParamFile.write("A) " + image["filename"] + 
			"						#Input data image block\n")
		galfitSingleParamFile.write("B) " + self.galfit_single_output_filename +
			"						#Output data image block\n")
		galfitSingleParamFile.write("C)" + " none" + 
			"						#Sigma image name (made from data if blank or 'none')\n")
		galfitSingleParamFile.write("D) " + self.psf + 
			"						#Input PSF image and (optional) diffusion kernel\n")
		galfitSingleParamFile.write("E)" + " 1" + 
			"						#PSF fine sampling factor relative to data\n")
		galfitSingleParamFile.write("F) none" + 					
			"						#Bad pixel mask (FITS file or ASCIIcoord list)\n")
		galfitSingleParamFile.write("G) none" + #self.constraintFilename + 
			"						#File with parameter constraints (ASCII file)\n")
		galfitSingleParamFile.write("H) " + "1" + " " + str(image["width"]) + " " + 
											"1" + " " + str(image["height"]) + 
			"						#Image region to fit (xmin xmax ymin ymax)\n")
		galfitSingleParamFile.write("I)" + " 200 200" + 
			"						#Size of the convolution box (x y)\n")
		galfitSingleParamFile.write("J) " + str(self.mpZeropoint) + 
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
		
		compNum = 1
		for model in image["models"]:
			galfitSingleParamFile.write("# Componenet number: " + str(compNum) + "\n")
			galfitSingleParamFile.write(" 0) sersic					#Component type\n")
			galfitSingleParamFile.write(" 1) " + str(model["centerCoords"][0]) + "	" + str(model["centerCoords"][1]) + "	1	1			#Position x,y\n")
			galfitSingleParamFile.write(" 3) " + str(model["magnitude"]) + " 1			#Integrated Magnitude\n")
			galfitSingleParamFile.write(" 4) " + str(model["radius"]) + "			1			#R_e (half-light radius)	[pix]\n")
			galfitSingleParamFile.write(" 5) " + str(model["sers"]) + "		1			#Sersic index n (de Vaucouleurs n=4)\n")
			galfitSingleParamFile.write(" 6) 0.0000		0			#	-----\n")
			galfitSingleParamFile.write(" 7) 0.0000		0			#	-----\n")
			galfitSingleParamFile.write(" 8) 0.0000		0			#	-----\n")
			galfitSingleParamFile.write(" 9) " + str(model["ba"]) + "			1			#Axis ratio (b/a)\n")
			galfitSingleParamFile.write(" 10) " + str(model["angle"]) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
			galfitSingleParamFile.write(" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
			galfitSingleParamFile.write("\n")
			compNum = compNum + 1
	
		galfitSingleParamFile.write("# Componenet number: " + str(compNum) + "\n")
		galfitSingleParamFile.write(" 0) sky						#Component type\n")
		galfitSingleParamFile.write(" 1) 0.0000		0			#	Sky background at center of fitting region [ADUs]\n")
		galfitSingleParamFile.write(" 2) 0.0000		0			#	dsky/dx (sky gradient in x) [ADUs/pix]\n")
		galfitSingleParamFile.write(" 3) 0.0000		0			#	dsky/dy (sky gradient in y) [ADUs/pix]\n")
		galfitSingleParamFile.write(" Z) 0						#Leave in [1] or subtract [0] this comp from data?\n")
		galfitSingleParamFile.write("\n")
	
		galfitSingleParamFile.close()
	
	
	def getCentermostID(self, imageHeight, imageWidth, resultFilename):
		'''
		return the integer ID of the galaxy closest to the images center
		'''
		with open(resultFilename, "r") as resultFile:
			resultContents = resultFile.readlines()
		
		closestDist = 1000.0
		closestID = -1
		currentID = -1
		skipSky = False
		for resultLine in resultContents:
		
			# toggle flag on the line indicating sky for use by future lines
			if resultLine.strip()[:6] == "A) sky":
				skipSky = True
			elif resultLine.strip()[:6] == "B) ser":
				skipSky = False
			
			# create an array version of the current result line
			resultLineList = resultLine.split()
			
			# get the component number if on the component header comment 
			if (resultLineList
				and (resultLine.strip()[0] == "#")
				and (resultLine.strip()[-1].isdigit())):
				currentID = int(resultLineList[-1])
			
			# if on the position line of a non-sky component, check if it is
			# the closest yet to the center of the image and, if so, save ID
			if not skipSky and resultLine.strip()[:2] == "1)":
				px = resultLineList[1]
				py = resultLineList[2]
				dx = imageWidth/2.0 - float(px)
				dy = imageHeight/2.0 - float(py)
				dist = math.sqrt(math.pow(dx,2) + math.pow(dy,2))
				if dist < closestDist:
					closestDist = dist
					closestID = currentID
		
		# return the ID of the component closest to the center of the image
		return closestID
	
	
	def write_galfit_bulge_parameter(self, centerID):
		'''
		reads the results of the first run and writes it to a new file
		with the output and constraint modified for bulge run and the
		new bulge component appended
		'''
		with open(self.galfit_single_result_filename, "r") as singleResultFile:
			resultContents = singleResultFile.readlines()
			
		# gather info about first fit while modifying output and constraint parameters
		bulgeParamStr = ""
		currentID = -1
		for resultLine in resultContents:
			resultLineList = resultLine.split()
				
			if resultLine.strip()[:2] == "B)":
				bulgeParamStr = (bulgeParamStr + "B) " + self.galfit_bulge_output_filename +
							"						#Output data image block\n")
				
			elif resultLine.strip()[:2] == "G)":
				bulgeParamStr = (bulgeParamStr + "G) " + self.galfit_bulge_constraint_filename + 
							"						#File with parameter constraints (ASCII file)\n")
							
			else:
			
				if (resultLineList
					and (resultLine.strip()[0] == "#")
					and (resultLine.strip()[-1].isdigit())):
					currentID = int(resultLineList[-1])
					bulgeParamStr = bulgeParamStr + resultLine
					
				if centerID == currentID:
					if resultLine.strip()[:2] == "1)":
						resultX = resultLineList[1]
						resultY = resultLineList[2]
						
					elif resultLine.strip()[:2] == "3)":
						resultMag = resultLineList[1]
						
					elif resultLine.strip()[:2] == "4)":
						resultRad = resultLineList[1]
						
					# rewrite the sersic index to be 1 and fixed for disk comp
					if resultLine.strip()[:2] == "5)":
						bulgeParamStr = (bulgeParamStr + " 5) " + 
							"1.0		0			#Sersic index n (de Vaucouleurs n=4)\n")
					else:
						bulgeParamStr = bulgeParamStr + resultLine
				else:
					bulgeParamStr = bulgeParamStr + resultLine
						
		
		# write the third component to the end of galfit_single_result_filename
		bulgeComponentID = currentID + 1
		bulgeParamStr = bulgeParamStr + ("# Componenet number: " + str(bulgeComponentID) + "\n")
		bulgeParamStr = bulgeParamStr + (" 0) sersic					#Component type\n")
		bulgeParamStr = bulgeParamStr + (" 1) " + resultX + " " + resultY + " 1 1			#Position x,y\n")
		bulgeParamStr = bulgeParamStr + (" 3) " + resultMag + " 1			#Integrated Magnitude\n")
		bulgeParamStr = bulgeParamStr + (" 4) " + resultRad + "			1			#R_e (half-light radius)	[pix]\n")
		bulgeParamStr = bulgeParamStr + (" 5) " + "2.5		1			#Sersic index n (de Vaucouleurs n=4)\n")
		bulgeParamStr = bulgeParamStr + (" 6) 0.0000		0			#	-----\n")
		bulgeParamStr = bulgeParamStr + (" 7) 0.0000		0			#	-----\n")
		bulgeParamStr = bulgeParamStr + (" 8) 0.0000		0			#	-----\n")
		bulgeParamStr = bulgeParamStr + (" 9) 1" + "			1			#Axis ratio (b/a)\n")
		bulgeParamStr = bulgeParamStr + (" 10) 0" + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
		bulgeParamStr = bulgeParamStr + (" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
		bulgeParamStr = bulgeParamStr + ("\n")
	
		# write the bulge parameter file using the modified contents of the single results
		with open(self.galfit_bulge_parameter_filename, "w") as galfitBulgeParamFile:
			galfitBulgeParamFile.write(bulgeParamStr)
			
		self.write_galfit_bulge_constraint(centerID, bulgeComponentID)
		
		
	def write_galfit_bulge_constraint(self, centerID, bulgeComponentID):
		'''
		simple method to write the constraint for given components position
		'''
		# write the bulge constraint file for two given components
		with open(self.galfit_bulge_constraint_filename, "w") as constraintFile:
			constraintFile.write(str(centerID) + "," + str(bulgeComponentID) + 
						" x -5  5	# Soft constraint: Constrains x position\n" +
						str(centerID) + "," + str(bulgeComponentID) + 
						" y -5  5	# Soft constraint: Constrains y position\n")
		
		
	def getGalfitNNFilename(self, multiFitsFilename):
		'''
		use the output of GALFIT to get the galfit.NN filename
		
		parameter multiFitsFilename - GALFIT's output (*multi.fits)
		returns - [gNN, errorFlag]
			gNN - galfit.NN filename  
			errorFlag - a string indicating if a numerical error has occurred
						empty string if no numerical error detected by GALFIT
		'''
		print ("getting galfit.NN filename from " + multiFitsFilename)
		# use pyfits to gather info from output of galfit
		imageHeaders = pyfits.open(multiFitsFilename)
		print ("opened output fits using pyfits")
		# get the filename from the dictionary mapping header keywords to their values
		gNN = imageHeaders[2].header["LOGFILE"]
		print ("got galfit.NN: " + gNN)
		if "2" in imageHeaders[2].header["FLAGS"].split():
			errorFlag = "# numerical error detected *"
		else:
			errorFlag = ""
		print ("checked for error in FLAGS")
		imageHeaders.close()
		print ("closed output fits file, returning")
		return [gNN, errorFlag]
		
		
	def run_galfit(self, image):
		'''
		runs galfit on the given image, using the image filename
		to define the filenames of all of the resulting files
		
		parameter image - the image on which galfit will be run
		'''
		
		# define filenames
		self.defineGalfitFilenames(image)
		
		# writes the single component parameter file
		self.write_galfit_single_parameter(image)
	
		# run galfit on paramter file
		os.system('galfit ' + self.galfit_single_parameter_filename)
		#os.system("rm fit.log")
			
		# detects atomic galfit error
		if not os.path.isfile(self.galfit_single_output_filename):
			self.logMsg = (	self.logMsg + 
							" galfit failed on single component, " + 
							"probably mushroom (atomic galfit error)")
			return

		# write the result using the resulting *_multi.fits[2] header
		[galfitNNFilename, errorFlag] = self.getGalfitNNFilename(
											self.galfit_single_output_filename)
		os.system("mv " + galfitNNFilename + " " + 
					self.galfit_single_result_filename)
		if errorFlag:
			with open (self.galfit_single_result_filename, 'a') as resultFile:
				resultFile.write(errorFlag)
		
		# done unless command line specified that a second galfit run
		# should be done by adding a bulge component to the result of the first run
		if self.includeBulgeComponent:
		
			# determine the centermost galaxy, on which the bulge component will be run
			centerID = self.getCentermostID(image["height"], image["width"],
											self.galfit_single_result_filename)
		
			# reads the results of the first run and writes it 
			# with the output and constraint modified for bulge run and the
			# new bulge component appended with some intitial guess parameters
			self.write_galfit_bulge_parameter(centerID)

			# run galfit on paramter file
			os.system('galfit ' + self.galfit_bulge_parameter_filename)
			#os.system("rm fit.log")
						
			# detects atomic galfit error
			if not os.path.isfile(self.galfit_bulge_output_filename):
				self.logMsg = (	self.logMsg + 
								" galfit failed on bulge component, " + 
								"probably mushroom (atomic galfit error)")
				return
			
			# write the result using the resulting *_multi.fits[2] header
			[galfitNNFilename, errorFlag] = self.getGalfitNNFilename(
												self.galfit_bulge_output_filename)
			os.system("mv " + galfitNNFilename + " " + 
						self.galfit_bulge_result_filename)
			if errorFlag:
				with open (self.galfit_bulge_result_filename, 'a') as resultFile:
					resultFile.write(errorFlag)
		
		# if we get here then nothing went wrong!
		self.logMsg = self.logMsg + " success"
		

	def modelImage(self, imageFilename):
		'''
		models a single image
		
		parameter imageFilename - the filename of the image to be modeled
		
		returns - the log string of the modeling
		''' 
			
		curImage = {"filename":imageFilename, "models":[]}
		self.logMsg = imageFilename + ": "
		
		# run galfit, preventing crashes but printing and logging errors in log
		try:
	
			# run iraf's imhead method to populate image dimensions
			self.run_imhead(curImage)
					
			# run sextractor, exiting if a necessary parameter is missing
			if not self.run_sextractor(curImage, self.sextractorConfigFilename):
				exit()
		
			# create a reg file for ds9 using the updated image dictionary
			self.write_ds9_reg_file(curImage)
			
			# run galfit if not suppressed by command line
			if not self.galfitOff:
				self.run_galfit(curImage)
		
		# catch, log, and ignore all runtime errors except explicit exits 
		# (for debugging). move on to the next image regardless
		except not (KeyboardInterrupt):
			errorMsg = str(sys.exc_info()[0]) + str(sys.exc_info()[1])
			print (errorMsg)
			self.logMsg = self.logMsg + errorMsg
		
		# not sure how well this is working in parallel
		except KeyboardInterrupt:
			print("User cancelled execution with ctrl-c")
			exit()
			
		return self.logMsg
		

def runModelGenerator(parameterList):
	'''
	uses the command line inputs gathered in __main__ to create an
	instance of the model generator class and invoke its methods
	
	parameter parameterList -
		wierd way of receiving parameters needed to facilitate parallel
		first four elements are parser, options, sexOptions, and destDirectory,
		respectively, and the rest is the list of image filenames to model
	'''
	
	# handle parameters this way to enable parallelism
	parser = parameterList[0]
	options = parameterList[1]
	sextractorKeywordOptions = parameterList[2]
	destDirectory = parameterList[3]
	imageFilenames = parameterList[4:]
	
	# holds methods for analyzing simulations
	modelGen = ModelGenerator(destDirectory=destDirectory)
	
	# parse the command line options
	modelGen.parseGalfitOptions(parser, options)
	modelGen.parseSextractorOptions(parser, options.realSextractor, 
									sextractorKeywordOptions)
	
	# write the sextractor config file, which will be used for all images
	modelGen.write_sextractor_config_file(modelGen.sextractorConfigFilename)
	modelGen.write_sextractor_config_file(modelGen.sextractorReduceComponentConfigFilename)

	# store log result of modeling each image using modelGen instance
	results = []
	for imageFilename in imageFilenames:
		if os.path.isfile(imageFilename.strip()):
			results.append(modelGen.modelImage(imageFilename.strip()))
		else:
			results.append(imageFilename.strip() + ": image file does not exist")

	# compose the log
	log = "run on " + time.strftime("%m-%d-%Y") + "\n"
	for result in results:
		log = log + result + "\n"
	
	# create the log file in the same directory as python was called
	try:
		logFilename = (modelGen.destDirectory + "rungalfit_log_" + 
			imageFilenames[0].split("/")[-1].split("_")[1].split(".")[1] + 
			"_to_" + 
			imageFilenames[-1].split("/")[-1].split("_")[1].split(".")[1] +
			"_" + time.strftime("%m-%d-%Y") + ".txt")
	except:
		logFilename = (modelGen.destDirectory + "rungalfit_log.txt")
	
	# write the log file
	print ("writing log file to " + logFilename)
	with open(logFilename, 'w') as logFile:
		logFile.write(log)
	
	return results

		
if __name__ == "__main__":
	'''
	parses the command line using the optparse package
	NOTE: argparse is the current version as of python 2.7, 
			but optparse is used to maintain better backwards compatibility
	'''
	
	#define the command line interface with simUtility.py
	usage = ("\n%prog inputFile [-h help] [options (with '-'|'--' prefix)]" +
			"  [sextractor options (no '-' prefix)]\n")
			
	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# fit an additional component after the initial GALFIT fit
	parser.add_option("-b","--bulge", 
				help="include to run a GALFIT bulge fit after the initial galaxy fit",
				action="store_true")
						
	# suppress running galfit after sextractor
	parser.add_option("-g","--galfitOff", 
				help="include to suppress running galfit after Source-Extractor",
				action="store_true")
						
	# run sextractor
	parser.add_option("-r","--realSextractor", 
				help="include to run Source-Extractor for real images, otherwise sim images assumed",
				action="store_true")
	
	# run in parallel
	parser.add_option("-p","--parallel", 
				help="include to run images in parallel, otherwise series",
				action="store_true")
	
	# the GALFIT constraint file. verified after parsing
	# TODO: this might be deprecated, since the component numbering
	# 		is less predictable since the introduction of sextractor
	#		Also constraint files are not recommended in GALFIT generally
	parser.add_option("-c","--constraint", 
				help="set the file constraining the GALFIT results")
			
	# Magnitude photometric zeropoint	
	parser.add_option("--mpz", metavar="MagnitudePhotometricZeropoint",
				type="float", default=26.3, 
				help="set the magnitude photometric zeropoint for" +
					" GALFIT to use [default: %default]")
						
	# Plate scale
	parser.add_option("--plate", metavar="PlateScale",
				type="float", default=0.06, 
				help="set the plate scale for GALFIT to use" +
						"[default: %default]")
								
	# PSF file. verified after parsing
	parser.add_option("--psf",
				help="set the file for GALFIT to use as a PSF")
						
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

	# read list of image filenames from input file
	with open(args[0], 'r') as inputFile:
		imageFilenames = inputFile.readlines()
		
	# verify there are images to model 
	numImages = len(imageFilenames)
	if numImages == 0:
		parser.error("input file " + args[0] +
					" has no contents (full path image filenames).")
					
	# verify that the calling computer has the necessary commands in PATH
	if os.system("sex --help"):
		print("Must have Source-Extractor's 'sex' command available by PATH " + 
				"environment variable. Exiting execution")
		exit()
	if os.system("galfit -help"):
		print("Must have GALFIT's 'galfit' command available by PATH " + 
				"environment variable. Exiting execution")
		exit()
	
	# TODO: might solve windows vs unix/mac, not sure though
	newFilenames = []
	for imageFilename in imageFilenames:
		newFilename = os.path.normpath(imageFilename.strip())
		newFilenames.append(newFilename)
		if not os.path.isfile(newFilename):
			parser.error("input file has filename " + newFilename + "\n"
						"This file does not exist or is not accessible\n" +
						"input file " + args[0] + 
						" must be only full path image filenames, one per line.")
	imageFilenames = newFilenames
	# done with immediate verifying of comamnd line #
		
	# for parallel, only use half the cpus available
	numCPUs = int(multiprocessing.cpu_count())
	
	# create the destination directory for results
	collectiveDestDirectory = "results_" + time.strftime("%m-%d-%Y-%H-%M-%S") + "/"
	os.mkdir(collectiveDestDirectory)
	
	# time how long it takes to run the program
	startTime = time.time()
	
	# only do parallel if rerquested and if enough images to warrant
	if not (options.parallel):# and (numImages >= numCPUs)):
		runModelGenerator([parser, options, args[1:], collectiveDestDirectory] + 
							imageFilenames)
	
	# run the parallel version
	else:
		
		# construct list, each element is a list of arguments for separate cpu
		imageArgs = []
		destDirectories = []
		chunkSize = int(math.ceil(float(numImages)/float(numCPUs)))
		print ("running in parallel, the list of images is being divided among" + 
				" your available processors in the following chunk sizes")
		for i in range(0, numImages, chunkSize):
			parDestDirectory = "results" + str(i)
			destDirectories.append(parDestDirectory)
			print(str(len(imageFilenames[i:i+chunkSize])) +
				" images being stored in the directory: " + parDestDirectory)
			imageArgs.append([parser, options, args[1:], parDestDirectory] + 
								imageFilenames[i:i+chunkSize])
								
		# see documentation on multiprocessing pool and map function
		print ("passing job to " + str(numCPUs) + " out of " + 
				str(multiprocessing.cpu_count()) + " CPUs, logs will be written as completed")
		pool = multiprocessing.Pool(numCPUs)
		results = pool.map(runModelGenerator, imageArgs)
		
		# compose the log
		log = ("run on " + time.strftime("%m-%d-%Y") + 
				" with command line input " + " ".join(sys.argv) + "\n")
		for result in results:
			for line in result:
				log = log + line + "\n"
		
		# create the log file in a collective destination directory
		logFilename = (collectiveDestDirectory + "combined_rungalfit_log_" + 
                    imageFilenames[0].split("/")[-1].split("_")[1].split(".")[1] + 
                    "_to_" + 
                    imageFilenames[-1].split("/")[-1].split("_")[1].split(".")[1] +
                    "_" + time.strftime("%m-%d-%Y") + ".txt")

		
		# write the log file
		print ("writing log file to " + logFilename)
		with open(logFilename, 'w') as logFile:
			logFile.write(log)
			
		# move all individual results into the new collective results folder
		for subDir in destDirectories:
			os.system(" ".join(["mv", subDir.rstrip("/") + "/*", collectiveDestDirectory]))
			os.rmdir(subDir)
	
	# TODO: this relies on all and only configuration files to start with "config"
	# move sextractor configuration files into the collective destination
	#os.system(" ".join(["mv","config*",collectiveDestDirectory]))
	if not options.galfitOff:
		os.system(" ".join(["mv","fit.log",collectiveDestDirectory]))
	
	# end program run time, print total
	elapsed = time.time() - startTime
	print(" ".join(["time","elapsed","=",str(int(elapsed)),"seconds",
			u"\u2245",str(int(elapsed/60.0)),"minutes"]))
			
	# unused code to enable version-free collective of user input
	'''
	# allow user to choose to overwrite the 
	if os.path.isdir(collectiveDestDirectory):
		
		# try to get user input for pre python 3.x users
		try:
			if raw_input("dest directory " + collectiveDestDirectory +
					"already exists.\nType Y and press ENTER to overwrite" +
					" or N to terminate program execution: "
					).strip().upper() == "Y":
				os.system(" ".join(["rm","-r",collectiveDestDirectory]))
			else:
				print("terminating.")
				exit()
		except:
			
			# try to get user input from post python 3.x users
			try:
				if input("dest directory " + collectiveDestDirectory +
						"already exists.\nType Y and press ENTER to overwrite" +
						" or N to terminate program execution: "
						).strip().upper() == "Y":
					os.system(" ".join(["rm","-r",collectiveDestDirectory]))
				else:
					print("terminating.")
					exit()
			
			# if all failed then just quit and alert the user
			except:
				print("dest directory " + collectiveDestDirectory + "already exists.")
				exit()
	'''
	
	# code for writing the result without needing the galfit NN file at all
	'''
		# not using, getting info from multi fits instead
		# because auto naming is a timing problem in parallel
		os.system("rm galfit.*")
		
		# detects atomic galfit error
		if not os.path.isfile(multiFitsFilename):
			self.logMsg = (	self.logMsg + 
							" galfit failed on single component, " + 
							"probably mushroom (atomic galfit error)")
			return False
		
		# use pyfits to gather info from output of galfit
		imageHeaders = pyfits.open(multiFitsFilename)
		
		# get the dictionary mapping header keywords to their values
		try:
			modelHeader = imageHeaders[2]
		except KeyError:
			self.logMsg = (	self.logMsg + 
							" getGalfitNNFIlename must be called on" +
							" a multi-extension cube (which galfit outputs)")
			return False
						
		# update image models to reflect galfit results
		resultModels = []
		for mIndex, model in enumerate(image["models"]):
			modelX = float(modelHeader[str(mIndex+1) + "_XC"].split()[0])
			modelY = float(modelHeader[str(mIndex+1) + "_YC"].split()[0])
			modelMag = float(modelHeader[str(mIndex+1) + "_MAG"].split()[0])
			modelRad = float(modelHeader[str(mIndex+1) + "_RE"].split()[0])
			modelSers = float(modelHeader[str(mIndex+1) + "_N"].split()[0])
			modelBA = float(modelHeader[str(mIndex+1) + "_AR"].split()[0])
			modelAng = float(modelHeader[str(mIndex+1) + "_PA"].split()[0])
			resultModels.append({	"centerCoords":[modelX,modelY],
									"magnitude":modelMag,
									"radius":modelRad,
									"ba":modelBA,
									"angle":modelAng,
									"sers":modelSers})
		image["models"] = resultModels
									
		# write info to the galfit result filename, using existing method
		singleParamFilename = self.galfit_single_parameter_filename
		self.galfit_single_parameter_filename = resultFilename
		self.write_galfit_single_parameter(image)
		
		# restore field value
		self.galfit_single_parameter_filename = singleParamFilename
		
		return True
	'''