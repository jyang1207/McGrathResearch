#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 3/9/2105
Colby College Astrophysics Research
'''

import os
import sys
import time
from optparse import OptionParser
from math import *
import pprint
try:
	import pyfits as fits
except ImportError:
	from astropy.io import fits
from PIL import Image
import numpy

# TODO: use this to make plots with age on the bottom x axis and z on the top
def ned_wright_cosmology_calculator(z):
	'''
	http://www.astro.ucla.edu/~wright/CC.python 
	
	parameter z - the redshift from which to calculate age in GYR and kpc per arcsec
	returns [zage_Gyr, kpc_DA]
	'''
	H0 = 69.6  # Hubble constant
	WM = 0.286  # Omega(matter)
	WV = 0.714  # Omega(vacuum) or lambda
	
	# initialize constants
	
	WR = 0.  # Omega(radiation)
	WK = 0.  # Omega curvaturve = 1-Omega(total)
	c = 299792.458  # velocity of light in km/sec
	Tyr = 977.8  # coefficent for converting 1/H into Gyr
	DTT = 0.5  # time from z to now in units of 1/H0
	DTT_Gyr = 0.0  # value of DTT in Gyr
	age = 0.5  # age of Universe in units of 1/H0
	age_Gyr = 0.0  # value of age in Gyr
	zage = 0.1  # age of Universe at redshift z in units of 1/H0
	zage_Gyr = 0.0  # value of zage in Gyr
	DCMR = 0.0  # comoving radial distance in units of c/H0
	DCMR_Mpc = 0.0 
	DCMR_Gyr = 0.0
	DA = 0.0  # angular size distance
	DA_Mpc = 0.0
	DA_Gyr = 0.0
	kpc_DA = 0.0
	DL = 0.0  # luminosity distance
	DL_Mpc = 0.0
	DL_Gyr = 0.0  # DL in units of billions of light years
	V_Gpc = 0.0
	a = 1.0  # 1/(1+z), the scale factor of the Universe
	az = 0.5  # 1/(1+z(object))
	
	h = H0 / 100.
	WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species, T0 = 2.72528
	WK = 1 - WM - WR - WV
	az = 1.0 / (1 + 1.0 * z)
	age = 0.
	n = 1000  # number of points in integrals
	for i in range(n):
		a = az * (i + 0.5) / n
		adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
		age = age + 1. / adot
	
	zage = az * age / n
	zage_Gyr = (Tyr / H0) * zage
	DTT = 0.0
	DCMR = 0.0
	
	# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
	for i in range(n):
		a = az + (1 - az) * (i + 0.5) / n
		adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
		DTT = DTT + 1. / adot
		DCMR = DCMR + 1. / (a * adot)
	
	DTT = (1. - az) * DTT / n
	DCMR = (1. - az) * DCMR / n
	age = DTT + zage
	age_Gyr = age * (Tyr / H0)
	DTT_Gyr = (Tyr / H0) * DTT
	DCMR_Gyr = (Tyr / H0) * DCMR
	DCMR_Mpc = (c / H0) * DCMR
	
	# tangential comoving distance
	
	ratio = 1.00
	x = sqrt(abs(WK)) * DCMR
	if x > 0.1:
		if WK > 0:
			ratio = 	 0.5 * (exp(x) - exp(-x)) / x 
		else:
			ratio = sin(x) / x
	else:
		y = x * x
		if WK < 0: y = -y
		ratio = 1. + y / 6. + y * y / 120.
	DCMT = ratio * DCMR
	DA = az * DCMT
	DA_Mpc = (c / H0) * DA
	kpc_DA = DA_Mpc / 206.264806
	DA_Gyr = (Tyr / H0) * DA
	DL = DA / (az * az)
	DL_Mpc = (c / H0) * DL
	DL_Gyr = (Tyr / H0) * DL
	
	# comoving volume computation
	
	ratio = 1.00
	x = sqrt(abs(WK)) * DCMR
	if x > 0.1:
		if WK > 0:
			ratio = (0.125 * (exp(2.*x) - exp(-2.*x)) - x / 2.) / (x * x * x / 3.)
		else:
			ratio = (x / 2. - sin(2.*x) / 4.) / (x * x * x / 3.)
	else:
		y = x * x
		if WK < 0: y = -y
		ratio = 1. + y / 5. + (2. / 105.) * y * y
	VCM = ratio * DCMR * DCMR * DCMR / 3.
	V_Gpc = 4.*pi * ((0.001 * c / H0) ** 3) * VCM
	
	return [zage_Gyr, kpc_DA]

									
def run_pyfits(multiFitsFilename):
	'''
	Read given fits file using pyfits and return data read from header and computed from image data
	
	parameter multiFitsFilename - filename of multiple extension cube fits file to extract data from
	returns [resultModels, imageWidth, imageHeight, kpcPerPixel, timeZ, wholeRFF, partialRFF]
	'''
	# use pyfits to gather info from output of galfit
	multiCubeSlices = fits.open(multiFitsFilename)
	sigmaImageName = "_".join(multiFitsFilename.split("_")[:-2]) + "_sigma.fits"
	sigmaSlices = fits.open(sigmaImageName) # only works if the sigma image is in the same directory
	
	# get the dictionary mapping header keywords to their values
	try:
		imageHeader = multiCubeSlices[1].header
		imageData = multiCubeSlices[1].data
		modelHeader = multiCubeSlices[2].header
		modelData = multiCubeSlices[2].data
		residualData = multiCubeSlices[3].data
		sigmaData = sigmaSlices[0].data
		# open a file for writing the pixels included in elliptical sum
		ellipImageName = multiFitsFilename[:-5] + "_ellip.png"
		ellipImage = Image.new("RGB", imageData.shape)
		multiCubeSlices.close()
		sigmaSlices.close()
	except KeyError:
		print("Result file " + 	multiFitsFilename + 
			" must be a .fits multi-extension cube with four slices")
		multiCubeSlices.close()
		exit()
		
	# get image dimensions
	imageWidth, imageHeight = imageData.shape
					
	# create list of dictionaries representing component models
	resultModels = []
	compNum = 1
	while ("COMP_" + str(compNum)) in modelHeader:
	
		# for the sky component grab its value and continue to next component
		if modelHeader["COMP_" + str(compNum)].strip() == "sky":
			resultModels.append({"sky": removeGalfitChars(str(modelHeader[str(compNum) + "_SKY"]).split()[0])})
			compNum = compNum + 1
			continue
		
		# initialize model to empty dictionary
		model = {}
			
		# get all model information for the current model
		compX = removeGalfitChars(str(modelHeader[str(compNum) + "_XC"].split()[0]))
		try:
			errX = removeGalfitChars(str(modelHeader[str(compNum) + "_XC"].split()[2]))
		except:
			errX = "0.0"
		model["px"] = [compX, errX]
		compY = removeGalfitChars(str(modelHeader[str(compNum) + "_YC"].split()[0]))
		try:
			errY = removeGalfitChars(str(modelHeader[str(compNum) + "_YC"].split()[2]))
		except:
			errY = "0.0"
		model["py"] = [compY, errY]
		compMag = removeGalfitChars(str(modelHeader[str(compNum) + "_MAG"].split()[0]))
		try:
			errMag = removeGalfitChars(str(modelHeader[str(compNum) + "_MAG"].split()[2]))
		except:
			errMag = "0.0"
		model["mag"] = [compMag, errMag]
		compRad = removeGalfitChars(str(modelHeader[str(compNum) + "_RE"].split()[0]))
		try:
			errRad = removeGalfitChars(str(modelHeader[str(compNum) + "_RE"].split()[2]))
		except:
			errRad = "0.0"
		model["rad"] = [compRad, errRad]
		compSer = removeGalfitChars(str(modelHeader[str(compNum) + "_N"].split()[0]))
		try:
			errSer = removeGalfitChars(str(modelHeader[str(compNum) + "_N"].split()[2]))
		except:
			errSer = "0.0"
		model["ser"] = [compSer, errSer]
		compBA = removeGalfitChars(str(modelHeader[str(compNum) + "_AR"].split()[0]))
		try:
			errBA = removeGalfitChars(str(modelHeader[str(compNum) + "_AR"].split()[2]))
		except:
			errBA = "0.0"
		model["ba"] = [compBA, errBA]
		compPA = removeGalfitChars(str(modelHeader[str(compNum) + "_PA"].split()[0]))
		try:
			errPA = removeGalfitChars(str(modelHeader[str(compNum) + "_PA"].split()[2]))
		except:
			errPA = "0.0"
		model["pa"] = [compPA, errPA]

		# use ellipe equation r^2 = (x/r)^2 + (Y/(r*ba))^2
		# rotated counterclockwise by position angle PA and translated to xc, yc
		# http://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
		PA = float(compPA)
		cosPA = cos(PA * 180.0 / pi)  # start at positive y not x
		sinPA = sin(PA * 180.0 / pi)  # start at positive y not x
		cossqPA = cosPA * cosPA
		sinsqPA = sinPA * sinPA
		a = 2.0 * float(compRad)
		b = a * float(compBA)
		invasq = 1 / (a * a)
		invbsq = 1 / (b * b)
		A = cossqPA * invasq + sinsqPA * invbsq  # rotation of ellipse by PA
		B = 2 * cosPA * sinPA * (invasq - invbsq)
		C = sinsqPA * invasq + cossqPA * invbsq
		h = float(compX)  # translating of ellipse center
		k = float(compY)

		# could compute bounding box for speed up http://stackoverflow.com/questions/87734/how-do-you-calculate-the-axis-aligned-bounding-box-of-an-ellipse
		'''
		t1 = atan2(-tan(PA),2)
		t2 = t1 + pi
		x1 = h + a*cos(t1)*cos(PA) - b*sin(t1)*sin(PA)
		x2 = h + a*cos(t2)*cos(PA) - b*sin(t2)*sin(PA)	
		if x1 < x2:
			xlow = x1
			xhigh = x2
		else:
			xlow = x2
			xhigh = x1
		t1 = atan2(1/tan(PA),2)
		t2 = t1 + pi
		y1 = k + b*sin(t1)*cos(PA) + a*cos(t1)*sin(PA)
		y2 = k + b*sin(t2)*cos(PA) + a*cos(t2)*sin(PA)
		if y1 < y2:
			ylow = y1
			yhigh = y2
		else:
			ylow = y2
			yhigh = y1
		'''
			
		# inside ellipse if:
		# A*x*x + B*x*y + C*y*y - (2*A*h + k*B)*x - (2*C*k + B*h)*y + (A*h*h + B*h*k + C*k*k - 1) < 0
		AhkB = 2 * A * h + k * B
		CkBh = 2 * C * k + B * h
		AhhBhkCkk = A * h * h + B * h * k + C * k * k - 1
		residualSum = 0
		sigmaSum = 0
		modelSum = 0
		pix = ellipImage.load()
		# http://mnras.oxfordjournals.org/content/419/3/2703.full.pdf for definition of rff
		for [y, x], imageVal in numpy.ndenumerate(imageData): # y, x because y is row and x is col
			if (# ((x>=xlow and x<=xhigh) and (y>=ylow and y<=yhigh)) and
				((A * x * x + B * x * y + C * y * y - AhkB * x - CkBh * y + AhhBhkCkk) < 0)):
				
				residualVal = residualData[y, x]
				modelVal = modelData[y, x]
				sigmaVal = sigmaData[y, x]
				
				residualSum += abs(imageVal-modelVal);  # or abs(residualVal)
				modelSum += modelVal
				sigmaSum += sigmaVal
				
				pix[x, imageHeight - y - 1] = 255  # flip so origin in bottom left

		# compute rff and store in model dictionary
		if modelSum:
			model["rff"] = str((residualSum - 0.8*sigmaSum) / modelSum)
		else:
			model["rff"] = "0"
			
		# append completed model to list of model dictionaries and advance loop
		resultModels.append(model)
		compNum = compNum + 1
			
	# compute the rff from the residual and image data
	ellipImage.save(ellipImageName)
	
	# return all extracted pyfits data 
	return [resultModels, imageWidth, imageHeight, imageHeader]


def getCentermostID(imageHeight, imageWidth, models):
	'''
	return the integer ID of the galaxy closest to the images center from the given models
	
	parameter imageHeight imageWidth - the dimensions of the image to compute center from
	parameter models - a dictionary of components from galfit for a single image
	returns the integer ID of the closest galaxy to the images center
	'''
	
	closestDist = 1000.0
	closestID = -1
	for currentID, model in enumerate(models):
		currentID = currentID + 1
		
		# skip sky
		if "sky" in model:
			continue
		# check if this model is the closest yet to the center of the image 
		px = model["px"][0]
		py = model["py"][0]
		dx = imageWidth / 2.0 - float(px)
		dy = imageHeight / 2.0 - float(py)
		dist = sqrt(pow(dx, 2) + pow(dy, 2))
		# if so, save ID
		if dist < closestDist:
			closestDist = dist
			closestID = currentID
	
	# return the ID of the component closest to the center of the image
	return closestID


def getNextCentermostID(imageHeight, imageWidth, models, centerID):
	'''
	return the integer ID of the second closest galaxy to the images center from given models
	
	parameter imageHeight imageWidth - the dimensions of the image to compute center from
	parameter models - a dictionary of components from galfit for a single image
	returns the integer ID of the second closest galaxy to the images center
	'''
	closestDist = 1000.0
	closestID = -1
	for currentID, model in enumerate(models):
		currentID = currentID + 1
		
		# skip sky
		if ("sky" in model) or (currentID == centerID):
			continue
		
		# check if this model is the closest yet to the center of the image 
		px = model["px"][0]
		py = model["py"][0]
		dx = imageWidth / 2.0 - float(px)
		dy = imageHeight / 2.0 - float(py)
		dist = sqrt(pow(dx, 2) + pow(dy, 2))
		# if so, save ID
		if dist < closestDist:
			closestDist = dist
			closestID = currentID
	
	# return the ID of the component closest to the center of the image
	return closestID


def removeGalfitChars(resultString):
	'''
	For a given galfit value as a string, strips galfit's annotations, leaving just the string number
	
	parameter resultString - the galfit string to strip annotations from
	returns given string without galfit annotation characters, jsut the number as a string
	'''
	return resultString.replace('*', '').replace('[', '').replace(']', ''
									).replace('{', '').replace('}', '')
											

def sum_galfit(resultFilename, models, imageHeader, delim, centerIDs, options):
	'''
	returns a string summary of the results in the given result filename, using remaining parameters
	to decorate the values given by GALFIT and as additional fields in each record
	
	parameter resultFilename - the result filename from running galfit to be summarized
	parameter models - 		a dictionary of components from galfit for a single image
	parameter imageHeader -	info gathered from the pyfits header
	parameter delim - 		the character(s) that delimit the fields in a single record (e.g. " ")
	parameter centerIDs - 	the list of ids of centermost and optionally next centermost galaxies
	parameter options - 	the dictionary of command line options specified at runtime
	returns a string of the results delimited by new lines, one line per component
	'''
	outputFilename = resultFilename.split("/")[-1].strip()
	# VELA02MRP_0.201015_0002949__skipir_CAMERA0-BROADBAND_F160W_simulation_bulge_multi.fits
	outputFilenameSplit = outputFilename.split("_")
	# [VELA02MRP, 0.201015, 0002949, , skipir, CAMERA0-BROADBAND, F160W, simulation, bulge, multi.fits
	#     0           1        2    3    4              5          6         7        8          9
	
	# variables to be parsed from the image header, default is zero
	kpcPerPixel = timeZ = mass = sfr = ssfr = 0.0
	if "SCALESIM" in imageHeader:
		kpcPerPixel = imageHeader["SCALESIM"]
	if "REDSHIFT" in imageHeader:
		timeZ = imageHeader["REDSHIFT"]
	if "MASS" in imageHeader:
		mass = imageHeader["MASS"]
	if "SFR" in imageHeader:
		sfr = imageHeader["SFR"]
	if "SSFR" in imageHeader:
		ssfr = imageHeader["SSFR"]
	if "GALAXYID" in imageHeader:
		galaxyID = imageHeader["GALAXYID"]
	else:
		try:
			galaxyID = outputFilenameSplit[0]
		except:
			galaxyID = ""
	# "VELA01"
	if "AVAL" in imageHeader:
		timeStep = imageHeader["AVAL"]
	else:
		try:
			timeStep = "0."+outputFilenameSplit[1].split(".")[1] #could have a0.#, so split
		except:
			timeStep = ""
	# "0.110"
	if "HALOID" in imageHeader:
		haloID = imageHeader["HALOID"]
	else:
		try:
			haloID = outputFilenameSplit[2]
		except:
			haloID = ""
	# "0002949"
	if "CAMERA" in imageHeader:
		camera = imageHeader["CAMERA"]
	else:
		try:
			camera = outputFilenameSplit[5].split("-")[0].replace("CAMERA", "")
		except:
			camera = ""
	# "0"
	if "FILTER" in imageHeader:
		filt = imageHeader["FILTER"]
	else:
		try:
			filt = outputFilenameSplit[6]
		except:
			filt = ""
	# "F125W"	
	
	# try to compute redshift from aval
	try:
		timeZ = float(timeZ)
	except:
		pass
	if timeStep and not timeZ:
		# a = 1/(1+z), z = 1/a - 1
		timeZ = 1.0 / float(timeStep) - 1.0
	
	# if the redshift is greater than zero, compute age gyr
	if timeZ:
		[age_gyr, kpcPerArcsec] = ned_wright_cosmology_calculator(timeZ)
		age_gyr = '%1.3f' % age_gyr
	else: # default values
		kpcPerArcsec = 0.0
		age_gyr = "0.0"
	timeZ = str(timeZ)
	
	# for computing the radius in kpc from the radius in pixels
	try:
		kpcPerPixel = float(kpcPerPixel)
	except:
		kpcPerPixel = 0.0
	# for candelized images, use ned wright and constant scale factor for kpc
	if options.candelized:
		kpcPerPixel = kpcPerArcsec * 0.06
		
	# the individual component properties
	componentList = []
	
	# the list of component lists
	componentLists = []
	
	# default sky value
	sky = "0.0"
	
	# for every model, make a record with GALFIT results and append to the string to be returned
	for currentID, model in enumerate(models, start=1):
		
		# sky component only yields one value, will be appended to all records after loop
		if "sky" in model:
			sky = str(model["sky"])
			continue
		
		# x position
		componentList.append(model["px"][0])  # value
		componentList.append(model["px"][1])  # error
		
		# y position
		componentList.append(model["py"][0])  # value
		componentList.append(model["py"][1])  # error
		
		# magnitude
		componentList.append(model["mag"][0])  # value
		componentList.append(model["mag"][1])  # error

		# radius (pixels)
		componentList.append(model["rad"][0])  # value
		componentList.append(model["rad"][1])  # error
		
		# radius (kpc)
		componentList.append(str(float(model["rad"][0]) * kpcPerPixel))  # value
		componentList.append(str(float(model["rad"][1]) * kpcPerPixel))  # error
			
		# sersic index
		componentList.append(model["ser"][0])  # value
		componentList.append(model["ser"][1])  # error
			
		# b/a
		componentList.append(model["ba"][0])  # value
		componentList.append(model["ba"][1])  # error
			
		# position angle
		componentList.append(model["pa"][0])  # value
		componentList.append(model["pa"][1])  # error

		# rff
		componentList.append(model["rff"])
		
		# type
		if currentID in centerIDs:
		
			# for single component fit the type is simply central or other
			if not options.bulge:
				galaxyType = "central"
				
			# when running bulge fit we fix disk component sersic index to 1
			elif model["ser"][0] == '1.0000':
				galaxyType = "disk"
				
			# for bulge fit non disk central component, conclude it is the bulge component
			else:
				galaxyType = "bulge"
				
		# non central components are of type other
		else:
			galaxyType = "other"
				
		# add this completed component as a line to the string of results
		componentLists.append([galaxyType] + componentList)
		
		# reset component list for any remaining components in this file
		componentList = []
		
	# add in the invariant fields after all components are done
	results = ""
	for component in componentLists:
		results += delim.join([outputFilename, galaxyID, haloID,
								timeStep, age_gyr, timeZ, camera, filt] + 
							component + [sky, sfr, ssfr, mass]) + "\n"
		
	# return resulting string
	return results


def main(args, pb=None):
	'''
	args - equivalent of sys.argv[1:]
	pb - optional tk.IntVar tracking progress, incremet after each file
	'''
	
	# define the command line interface with simUtility.py
	usage = ("\nsumSimUtility.py resultsFile [-h help] [options (with '-'|'--' prefix)]")

	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-b", "--bulge",
				help="if running on results of GALFIT bulge (two component) fit",
				action="store_true")
	
	parser.add_option("-d", "--delim",
				help=("set the delimiter to separate the fields of the summary file " + 
					"[default: %default]"),
				default=" ")

	# indicate that images are candelized
	parser.add_option("-r", "--candelized",
				help="if running on candelized images",
				action="store_true")
	
	parser.add_option("-o", "--output",
				help=("set the filename to write the output summary file " + 
					"[default: %default]"),
				default="summary_" + time.strftime("%m-%d-%Y") + ".csv")

	# indicate that the program should print actions to console
	parser.add_option("-v", "--verbose",
                help="to enable command line printouts of state",
                action="store_true")
	
	# parse the command line using above parameter rules
	# options - list with everthing defined above, 
	# args - anything left over after parsing options
	[options, args] = parser.parse_args(args)
	pprint.pprint(vars(options))
	print(args)

	if not (len(args) and os.path.isfile(args[0])):
		parser.error("must give an existing file with result filenames")
	else:
		inputFilename = args[0]
		
	# this will be the file that will contain the images
	with open(inputFilename, 'r') as r:
		resultFilenames = r.readlines()
	
	# open the specified filename for writing
	outFile = open(options.output, 'w')
	
	# the delimiter of the summary file
	delim = options.delim
	
	# the summary file header
	outFile.write(delim.join(["filename", "type", "galaxyID", "haloID"
							"timeStep", "age(GYr)", "redshift(z)", 
							"camera", "filter",
							"px(pixels)", "errpx(pixels)",
							"py(pixels)", "errpy(pixels)",
							"mag", "errmag",
							"rad(pixels)", "errrad(pixels)",
							"rad(kpc)", "errrad(kpc)",
							"sersicIndex(n)", "errsersicIndex(n)",
							"b/a", "errb/a",
							"angle(deg)", "errangle(deg)",
							"rff", "sky", "sfr", "ssfr", "mass"
							]) + "\n" + 
				delim.join(["string", "enum", "enum", "numeric", "numeric",
						"numeric", "enum", "enum"]) + "\n")
	
	# this loops through every result, writing summary to output
	for resultFilename in resultFilenames:
	
		# remove the new line character and any leading or trailing white space
		resultFilename = resultFilename.strip()
		if options.verbose: print("processing result file: " + resultFilename)
		
		# verify that result file exists
		if not os.path.isfile(resultFilename):
			parser.error("result file " + resultFilename + 
						" must be an existing file with result filenames")
	
		# collect info from result file header using pyfits
		resultInfo = run_pyfits(resultFilename)	
		
		# unpack the result into separate variables
		models, imageWidth, imageHeight, imageHeader = resultInfo
		
		# get the id of the centermost galaxy or galaxies if a bulge is included
		centerID = getCentermostID(imageWidth, imageHeight, models)
		nextCenterID = getNextCentermostID(imageWidth, imageHeight, models, centerID)
		if options.bulge:
			centerIDs = [centerID, nextCenterID]
		else:
			centerIDs = [centerID]
			
		# summarize galfit and write to output
		outFile.write(sum_galfit(resultFilename, models, imageHeader, delim, centerIDs, options))
		if pb: # progress bar in dashboard GUI increment
			pb.set(pb.get()+1)
		
	outFile.close()
	
	print ("Output written to " + options.output)

if __name__ == "__main__":
	main(sys.argv[1:])
