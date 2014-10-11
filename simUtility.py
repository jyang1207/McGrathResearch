#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath, Ariunjargal Bat-Erdene '15, Ryan Cole '15
Last Edited: 8/5/2104
Colby College Astrophysics Research
'''
import os
import sys
from pyraf import iraf # only used for imhead
import multiprocessing # to run in parallel using all available cpus
import time # to print runtime in log
from optparse import OptionParser
import pprint
		
class ModelGenerator:
	'''
	The controller class, which holds methods for analyzing simulations
	'''
		
	def __init__(self, callingDirectory=""):
		'''constructor sets some initial field defaults'''
		self.outputCatFilename = "generic.cat" # -CATALOG_NAME <filename>
		self.segmentationMapFilename = "check.fits" # -CHECKIMAGE_NAME <filename>
		self.sextractorOptionsList = []
		# "/".join(image.filename.split("/")[:-1]) + "/"
	
		# cd into desination directory to prevent processor collisions
		destDirectory = os.path.join(callingDirectory,"_" + str(os.getpid()))
		if not os.path.isdir(destDirectory):
			os.mkdir(destDirectory)
		os.chdir(destDirectory)
		self.callingDirectory = callingDirectory
		self.destDirectory = destDirectory
			
		# set the names of the sextractor configuration files
		self.sextractorConfigFilename = os.path.join(destDirectory,"configInit.sex")
		self.sextractorReduceComponentConfigFilename = os.path.join(destDirectory, 
														"configFewerComp.sex")
		
		# some other required sextractor files
		requiredFilesDirectory = os.path.join(callingDirectory,"requiredFiles")
		self.sextractorFilterFilename = os.path.join(requiredFilesDirectory,
													"sex.conv")
		self.sextractorParamFilename = os.path.join(requiredFilesDirectory,
												"sex.param")
		self.sextractorNNWFilename = os.path.join(requiredFilesDirectory,
												"default.nnw")
		if not os.path.isfile(self.sextractorFilterFilename):
			print("Program needs the required file " + 
				self.sextractorFilterFilename + " to exist.")
			exit()
		if not os.path.isfile(self.sextractorParamFilename):
			print("Program needs the required file " + 
				self.sextractorParamFilename + " to exist.")
			exit()
		if not os.path.isfile(self.sextractorNNWFilename):
			print("Program needs the required file " + 
				self.sextractorNNWFilename + " to exist.")
			exit()
		
		
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

		# galfit sigma image default none unless one is given on command line
		if not options.sigmaImage:
			self.sigmaImage = "none"
		else:
			# adjust relative path or make full path
			if options.sigmaImage[:2] == "..":
				sigmaImage = os.path.join("..",options.sigmaImage)
			else:
				sigmaImage = os.path.join(self.callingDirectory,options.sigmaImage)
				
			# verify that file exists
			if not os.path.isfile(sigmaImage):
				parser.error("sigmaImage " + sigmaImage + 
							" either does not exist or is not accessible")
			else:
				self.sigmaImage = sigmaImage

		# galfit psf default none unless one is given on command line
		if not options.psf:
			self.psf = "none"
		else:
			# adjust relative path or make full path
			if options.psf[:2] == "..":
				psf = os.path.join("..",options.psf)
			else:
				psf = os.path.join(self.callingDirectory,options.psf)
				
			# verify that file exists
			if not os.path.isfile(psf):
				parser.error("psf " + psf + 
							" either does not exist or is not accessible")
			else:
				self.psf = psf
		
		
	def parseSextractorOptions(self, parser, realSextractor, sextractorOptions):
		'''
		parses the sextractor options from the command line and uses to
		initialize fields with options or defaults
		'''
				
		# store running sextractor booleans
		self.realSextractor = realSextractor
		
		# all positional arguments after image list filename must come in pairs
		if len(sextractorOptions) % 2 != 0:
			parser.error("Sextractor options parsing error, odd number of options. " + 
				"Expecting sextractor options to be <keyword> <value> pairs")
		
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
		
		# variables describing the sextractor config file
		catalogName = "generic.cat" # Name of the output catalogue. If the name "STDOUT" is given and CATALOG TYPE is set to ASCII, ASCII HEAD, ASCII SKYCAT, or ASCII VOTABLE the catalogue will be piped to the standard output (stdout
		catalogType = "ASCII_HEAD"	# Format of output catalog:
									# ASCII - ASCII table; the simplest, but space and time consuming,
									# ASCII HEAD - as ASCII, preceded by a header containing information about the content,
									# ASCII SKYCAT - SkyCat ASCII format (WCS coordinates required)
									# ASCII VOTABLE - XML-VOTable format, together with meta-data,
									# FITS 1.0 - FITS format as in SExtractor 1
									# FITS LDAC - FITS "LDAC" format (the original image header is copied).
		sextractorParamFilename = self.sextractorParamFilename
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
		if self.realSextractor: #TODO: look at these numbers
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
		filterName = self.sextractorFilterFilename # Name and path of the file containing the filter definition
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
		starNNWFilename = self.sextractorNNWFilename #Name of the file containing the neural network weights for star/galaxy separation.
		#------------------------------ Background -----------------------------------
		backSize = "600"	# Size, or Width, Height (inpixels) of a background mesh.
		backFilterSize = "9"	#Size, or Width, Height (inbackground meshes) of the background-filtering mask.
		backPhotoType = "LOCAL" #Background used to compute magnitudes:
								# GLOBAL - taken directly from the background map
								# LOCAL - recomputed in a rectangular annulus around the object
		backPhotoThickness = "100"	#Thickness (in pixels) of the background LOCAL annulus.
		backType = "AUTO"	# What background is subtracted from the images: 
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
				# MANUAL A user-supplied constant value provided in BACK VALUE.
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
		self.outputCatFilename = os.path.join(self.destDirectory,
									".".join(image["filename"].split("/")[-1
										].split(".")[:-1]) + ".cat")
		self.segmentationMapFilename = os.path.join(self.destDirectory,
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
		
		# if sextractor does not find any galaxies, return
		if not numGalaxies:
			return False
		
		# if there are too many galaxies after the first fit, recursively
		# run sextractor again, now with a different configuration file
		if not self.realSextractor and (configFilename == self.sextractorConfigFilename) and (numGalaxies > 3):
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
				
				# exit program execution if the parameter file is incomplete
				if errorStr:
					print(errorStr)
					exit()

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
		regFileName = os.path.join(self.destDirectory,
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
		
		# detect an iraf error
		if imhead_return[-1] != ":":
			return False
		
		imhead_info = imhead_return.split("/")[-1]
		# "VELA01_a0.110_0006317__skipir_CAMERA0-BROADBAND_F125W_simulation.fits[600,600][real]:" 
		
		
		frame_dimensions = imhead_info.split("[")[1].replace("]", "").split(",")
		# ["600","600"]
	
		image["width"] = float(frame_dimensions[0])
		# "600"
		
		image["height"] = float(frame_dimensions[1])
		# "600"
		
		image["models"] = []
		
		return True
	
		
	def defineGalfitFilenames(self, image):
		filename = os.path.join(self.destDirectory, 
					".".join(image["filename"].split("/")[-1].split(".")[:-1]))
		self.galfit_single_constraint_filename =filename + '_single_constraint.txt'
		self.galfit_single_parameter_filename = filename + '_single_param.txt'
		self.galfit_single_output_filename =	filename + "_single_multi.fits"
		self.galfit_single_result_filename =	filename + "_single_result.txt"
		self.galfit_bulge_constraint_filename = filename + '_bulge_constraint.txt'
		self.galfit_bulge_parameter_filename =	filename + '_bulge_param.txt'
		self.galfit_bulge_output_filename =		filename + "_bulge_multi.fits"
		self.galfit_bulge_result_filename =		filename + "_bulge_result.txt"

	
	def write_galfit_parameter(self, image, paramFile, ouputFilename, sigmaFilename):
		'''
		writes the galfit parameter file for the single component
		
		parameter image - 
			the image on which galfit is being run
		parameter paramFile - 
			the open file object to write parameter file contents to
		parameter outputFilename - 
			the filename of the ouput that GALFIT will produce when running
			the parameter file written by this method call
		'''
		paramFile.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
		paramFile.write("A) " + image["filename"] + 
			"						#Input data image block\n")
		paramFile.write("B) " + ouputFilename +
			"						#Output data image block\n")
		paramFile.write("C) " + sigmaFilename + 
			"						#Sigma image name (made from data if blank or 'none')\n")
		paramFile.write("D) " + self.psf + 
			"						#Input PSF image and (optional) diffusion kernel\n")
		paramFile.write("E)" + " 1" + 
			"						#PSF fine sampling factor relative to data\n")
		paramFile.write("F) none" + 					
			"						#Bad pixel mask (FITS file or ASCIIcoord list)\n")
		paramFile.write("G) none" +  
			"						#File with parameter constraints (ASCII file)\n")
		paramFile.write("H) " + "1" + " " + str(image["width"]) + " " + 
											"1" + " " + str(image["height"]) + 
			"						#Image region to fit (xmin xmax ymin ymax)\n")
		paramFile.write("I)" + " 200 200" + 
			"						#Size of the convolution box (x y)\n")
		paramFile.write("J) " + str(self.mpZeropoint) + 
			"						#Magnitude photometric zeropoint\n")
		paramFile.write("K) " + str(self.plateScale) + "	 " + str(self.plateScale) + 
			"						#Plate scale (dx dy)  [arcsec per pixel]\n")
		paramFile.write("O)" + " regular" + 
			"						#display type (regular, curses, both\n")
		paramFile.write("P)" + " 0" + 
			"						#Options: 0=normal run; 1,2=make model/imgblock & quit\n")
		paramFile.write("S)" + " 0" + 
			"						#Modify/create components interactively?\n")
		paramFile.write("\n")
		paramFile.write("# INITIAL FITTING PARAMETERS\n")
		paramFile.write("#\n")
		paramFile.write("# For component type, the allowed functions are:\n")
		paramFile.write("#		nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n")
		paramFile.write("#		ferrer, coresersic, sky and isophote.\n")
		paramFile.write("#\n")
		paramFile.write("# Hidden parameters will only appear when they are specified:\n")
		paramFile.write("#		C0 (diskyness/boxyness),\n")
		paramFile.write("#		Fn (n=interger, Azimuthal Fourier Modes).\n")
		paramFile.write("#		R0-R10 (PA rotation, for creating spiral structures).\n")
		paramFile.write("#\n")
		paramFile.write("# ------------------------------------------------------------------------------\n")
		paramFile.write("#		par)	par value(s)	fit toggle(s)	# parameter description\n")
		paramFile.write("# ------------------------------------------------------------------------------\n")
		paramFile.write("\n")
		
		compNum = 1
		for model in image["models"]:
			paramFile.write("# Componenet number: " + str(compNum) + "\n")
			paramFile.write(" 0) sersic					#Component type\n")
			paramFile.write(" 1) " + str(model["centerCoords"][0]) + "	" + str(model["centerCoords"][1]) + "	1	1			#Position x,y\n")
			paramFile.write(" 3) " + str(model["magnitude"]) + " 1			#Integrated Magnitude\n")
			paramFile.write(" 4) " + str(model["radius"]) + "			1			#R_e (half-light radius)	[pix]\n")
			paramFile.write(" 5) " + str(model["sers"]) + "		1			#Sersic index n (de Vaucouleurs n=4)\n")
			paramFile.write(" 6) 0.0000		0			#	-----\n")
			paramFile.write(" 7) 0.0000		0			#	-----\n")
			paramFile.write(" 8) 0.0000		0			#	-----\n")
			paramFile.write(" 9) " + str(model["ba"]) + "			1			#Axis ratio (b/a)\n")
			paramFile.write(" 10) " + str(model["angle"]) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
			paramFile.write(" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
			paramFile.write("\n")
			compNum = compNum + 1
	
		paramFile.write("# Componenet number: " + str(compNum) + "\n")
		paramFile.write(" 0) sky						#Component type\n")
		if self.realSextractor:
			paramFile.write(" 1) 0.0000		1			#	Sky background at center of fitting region [ADUs]\n")
		else:
			paramFile.write(" 1) 0.0000		0			#	Sky background at center of fitting region [ADUs]\n")
		paramFile.write(" 2) 0.0000		0			#	dsky/dx (sky gradient in x) [ADUs/pix]\n")
		paramFile.write(" 3) 0.0000		0			#	dsky/dy (sky gradient in y) [ADUs/pix]\n")
		paramFile.write(" Z) 0						#Leave in [1] or subtract [0] this comp from data?\n")
		paramFile.write("\n")
	
	
	def getCentermostID(self, imageHeight, imageWidth, resultFilename):
		'''
		return the integer ID of the galaxy closest to the images center
		'''
		with open(resultFilename, "r") as resultFile:
			resultContents = resultFile.readlines()
		
		closestDist = 0.0 # same as false in python
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
				dist = dx*dx + dy*dy
				if not closestDist or (dist < closestDist):
					closestDist = dist
					closestID = currentID
		
		# return the ID of the component closest to the center of the image
		return closestID
	
	
	def write_galfit_bulge_parameter(self, centerID):
		'''
		reads the logResults of the first run and writes it to a new file
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
	
		# write the bulge parameter file using the modified contents of the single logResults
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
		
		
	def run_galfit(self, image):
		'''
		runs galfit on the given image, using the image filename
		to define the filenames of all of the resulting files
		
		parameter image - the image on which galfit will be run
		'''
		
		# define filenames
		self.defineGalfitFilenames(image)
		
		# writes the single component parameter file
		with open(self.galfit_single_parameter_filename, 'w') as paramFile:
			self.write_galfit_parameter(image, paramFile, 
										self.galfit_single_output_filename,
										self.sigmaImage)
	
		print(os.getcwd())
		# run galfit on paramter file
		os.system('galfit ' + self.galfit_single_parameter_filename)

		print(os.getcwd())
		# detects atomic galfit error
		if not os.path.isfile(self.galfit_single_output_filename):
			self.logMsg = (	self.logMsg + 
							" galfit failed on single component, " + 
							"probably mushroom (atomic galfit error)")
			return
		
		# rename galfit.01 to result file
		if not os.path.isfile("galfit.01"):
			self.logMsg = self.logMsg + " galfit.01 DNE"
			return

		os.system(" ".join(["mv","galfit.01",self.galfit_single_result_filename]))
		
		# done unless command line specified that a second galfit run
		# should be done by adding a bulge component to the result of the first run
		if self.includeBulgeComponent:
		
			# determine the centermost galaxy, on which the bulge component will be run
			centerID = self.getCentermostID(image["height"], image["width"],
											self.galfit_single_result_filename)
		
			# reads the logResults of the first run and writes it 
			# with the output and constraint modified for bulge run and the
			# new bulge component appended with some intitial guess parameters
			self.write_galfit_bulge_parameter(centerID)

			# run galfit on paramter file
			os.system('galfit ' + self.galfit_bulge_parameter_filename)
			
			# detects atomic galfit error
			if not os.path.isfile(self.galfit_bulge_output_filename):
				self.logMsg = (	self.logMsg + 
								" galfit failed on bulge component, " + 
								"probably mushroom (atomic galfit error)")
				return
		
			# rename galfit.01 to result file
			if not os.path.isfile("galfit.01"):
				self.logMsg = self.logMsg + " galfit.01 DNE"
				return
			
			# rename galfit.NN to result file
			os.system(" ".join(["mv","galfit.01",self.galfit_bulge_result_filename]))
		
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
			if not self.run_imhead(curImage):
				errorMsg = "IRAf had an error runnign imhead"
				print(errorMsg)
				self.logMsg = self.logMsg + errorMsg
				return self.logMsg
					
			# run sextractor, returning if it doesn't find any galaxies
			if not self.run_sextractor(curImage, self.sextractorConfigFilename):
				errorMsg = "skipping image because no galaxies found by Source Extractor"
				print(errorMsg)
				self.logMsg = self.logMsg + errorMsg
				return self.logMsg
		
			# create a reg file for ds9 using the updated image dictionary
			self.write_ds9_reg_file(curImage)
			
			# run galfit if not suppressed by command line
			if not self.galfitOff:
				self.run_galfit(curImage)
		
		# catch, log, and ignore all runtime errors except explicit exits 
		# (for debugging). move on to the next image regardless
		except not (KeyboardInterrupt or SystemExit):
			errorMsg = str(sys.exc_info()[0]) + str(sys.exc_info()[1])
			print (errorMsg)
			self.logMsg = self.logMsg + errorMsg
		
		# not sure how well this is working in parallel
		except KeyboardInterrupt:
			print("User cancelled execution with ctrl-c")
			exit()
			
		return self.logMsg
		

def runModelGeneratorSerial(parser, options, sextractorKeywordOptions, 
							callingDirectory, imageFilenames):
	'''
	uses the command line inputs gathered in __main__ to create an
	instance of the model generator class and invoke its methods
	
	parameter parser
	parameter options
	parameter sextractorKeywordOptions
	parameter callingDirectory
	parameter imageFilename
	
	returns - the list of lines to be written to the log
	'''
	
	# holds methods for analyzing simulations
	modelGen = ModelGenerator(callingDirectory=callingDirectory)
	
	# parse the command line options
	modelGen.parseGalfitOptions(parser, options)
	modelGen.parseSextractorOptions(parser, options.realSextractor, 
									sextractorKeywordOptions)
	
	# write the sextractor config file, which will be used for all images
	modelGen.write_sextractor_config_file(modelGen.sextractorConfigFilename)
	modelGen.write_sextractor_config_file(modelGen.sextractorReduceComponentConfigFilename)

	# store log result of modeling each image using modelGen instance
	logResults = []
	for imageFilename in imageFilenames:
			
		# run modeling method
		logResults.append(modelGen.modelImage(imageFilename.strip()))
	
	return [logResults, modelGen.destDirectory]

	
def runModelGeneratorParallel(parameterList):
	'''
	uses the command line inputs gathered in __main__ to create an
	instance of the model generator class and invoke its methods
	
	parameter parameterList -
		wierd way of receiving parameters needed to facilitate parallel
		the elements are parser, options, sextractorKeywordOptions, 
		callingDirectory, and the image filename to model, respectively
	'''
	
	# handle parameters this way to enable parallelism
	parser = parameterList[0]
	options = parameterList[1]
	sextractorKeywordOptions = parameterList[2]
	callingDirectory = parameterList[3]
	imageFilename = parameterList[4]
	
	# holds methods for analyzing simulations
	modelGen = ModelGenerator(callingDirectory=callingDirectory)
	
	# parse the command line options
	modelGen.parseGalfitOptions(parser, options)
	modelGen.parseSextractorOptions(parser, options.realSextractor, 
									sextractorKeywordOptions)
	
	# write the sextractor config file, which will be used for all images
	modelGen.write_sextractor_config_file(modelGen.sextractorConfigFilename)
	modelGen.write_sextractor_config_file(modelGen.sextractorReduceComponentConfigFilename)
		
	# run modeling method and return the resulting line in the log
	return [modelGen.modelImage(imageFilename.strip()), modelGen.destDirectory]


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
	
	# run in parallel
	parser.add_option("-p","--parallel", 
				help="include to run images in parallel, otherwise series",
				action="store_true")
						
	# run sextractor
	parser.add_option("-r","--realSextractor", 
				help="include to run Source-Extractor for real images, otherwise sim images assumed",
				action="store_true")
	
	# the GALFIT sigma image
	parser.add_option("-s","--sigmaImage", 
				help="set the file defining the GALFIT sigma image")
			
	# Magnitude photometric zeropoint	
	parser.add_option("--mpz", metavar="MagnitudePhotometricZeropoint",
				type="float", default=26.23, 
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
	pprint.pprint(vars(options))
	
	# real images need psfs, warn if no psf specified for real
	if options.realSextractor and not options.psf:
		print("Warning: no psf given for candelized images, are you sure you want to continue?")
		# TODO: only works for pre 3 python
		if raw_input("Type 'yes' to continue, otherwise program will end: ").upper() != "YES":
			exit()
	
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
	if not options.galfitOff and os.system("galfit -help"):
		print("Must have GALFIT's 'galfit' command available by PATH " + 
				"environment variable. Exiting execution")
		exit()
	
	# TODO: might solve windows vs unix/mac, not sure though
	newFilenames = []
	for imageFilename in imageFilenames:
		if not imageFilename.strip():
			continue
		newFilename = os.path.normpath(imageFilename.strip())
		if not os.path.isfile(newFilename):
			parser.error("input file has filename " + newFilename + "\n"
						"This file does not exist or is not accessible\n" +
						"input file " + args[0] + 
						" must be only full path image filenames, one per line.")
			
		# adjust for image filename if relative, otherwise make full path
		if newFilename[:2] == "..":
			newFilename = os.path.join("..",newFilename)
		else:
			newFilename = os.path.join(os.getcwd(),newFilename)
		newFilenames.append(newFilename)
	imageFilenames = newFilenames
	
	# done with immediate verifying of comamnd line #
		
	# for parallel, only use half the cpus available
	numCPUs = int(multiprocessing.cpu_count())
	
	# create the destination directory for logResults
	collectiveDestDirectory = os.path.join(os.getcwd(),
								"results_" + time.strftime("%m-%d-%Y-%T"))
	os.mkdir(collectiveDestDirectory)
	
	# time how long it takes to run the program
	startTime = time.time()
	
	# only do parallel if rerquested and if enough images to warrant
	logResults = []
	if not (options.parallel and (numImages >= numCPUs)):
		print("running in serial")
		[logResults,destDirectory] = runModelGeneratorSerial(parser, options, 
															args[1:], 
															os.getcwd(), 
															imageFilenames)
		os.system(" ".join(["mv", os.path.join(destDirectory,"*"), collectiveDestDirectory]))
		os.system(" ".join(["rmdir",destDirectory]))
	
	# run the parallel version
	else:
		
		# construct list, each element is a list of arguments for separate cpu
		imageArgs = []
		print ("running in parallel")
		for imageFilename in imageFilenames:
			imageArgs.append([	parser, options, args[1:], 
								os.getcwd(), imageFilename])
								
		# see documentation on multiprocessing pool and map function
		print ("passing job to " + str(numCPUs) + " out of " + 
				str(multiprocessing.cpu_count()))
		pool = multiprocessing.Pool(numCPUs)
		allResults = pool.imap_unordered(runModelGeneratorParallel, imageArgs)
		pool.close()
		pool.join()
		
		# separate out the log file results from the directories returned
		logResults = []
		destDirectories = []
		for result in allResults:
			logResults.append(result[0])
			destDirectories.append(result[1])
		destDirectories = set(destDirectories)
		
		# move all individual logResults into the new collective logResults folder
		for destDirectory in destDirectories:
			os.system(" ".join(["mv", os.path.join(destDirectory,"*"), collectiveDestDirectory]))
			os.system(" ".join(["rmdir",destDirectory]))
		
	# end program run time, print total
	elapsed = time.time() - startTime
	print(" ".join(["total","time","elapsed","=",str(int(elapsed)),"seconds",
			"about","equal","to",str(int(elapsed/60.0)),"minutes"]))

	# compose the log
	log = ("run on " + time.strftime("%m-%d-%Y-%T") + 
			" with command line input " + " ".join(sys.argv + 
			["total","time","elapsed","=",str(int(elapsed)),"seconds",
			"about","equal","to",str(int(elapsed/60.0)),"minutes",
			"about","equal","to",str(int((elapsed/60.0)/60.0)),"hours"]) + "\n")
	for logLine in logResults:
		log = log + logLine + "\n"
	
	# create the log file in a collective destination directory
	logFilename = os.path.join(collectiveDestDirectory, "rungalfit_log.txt")
	
	# write the log file
	print ("writing log file to " + logFilename)
	with open(logFilename, 'w') as logFile:
		logFile.write(log)
		
