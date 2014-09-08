#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 8/5/2104
Colby College Astrophysics Research
'''

import os
import time
from optparse import OptionParser
from math import sqrt, exp, sin, pi, pow
import pprint
import pyfits

	
def ned_wright_cosmology_calculator(z):
	'''
	http://www.astro.ucla.edu/~wright/CC.python	
	'''
	H0 = 69.6                         # Hubble constant
	WM = 0.286                        # Omega(matter)
	WV = 0.714						  # Omega(vacuum) or lambda
	
	# initialize constants
	
	WR = 0.        # Omega(radiation)
	WK = 0.        # Omega curvaturve = 1-Omega(total)
	c = 299792.458 # velocity of light in km/sec
	Tyr = 977.8    # coefficent for converting 1/H into Gyr
	DTT = 0.5      # time from z to now in units of 1/H0
	DTT_Gyr = 0.0  # value of DTT in Gyr
	age = 0.5      # age of Universe in units of 1/H0
	age_Gyr = 0.0  # value of age in Gyr
	zage = 0.1     # age of Universe at redshift z in units of 1/H0
	zage_Gyr = 0.0 # value of zage in Gyr
	DCMR = 0.0     # comoving radial distance in units of c/H0
	DCMR_Mpc = 0.0 
	DCMR_Gyr = 0.0
	DA = 0.0       # angular size distance
	DA_Mpc = 0.0
	DA_Gyr = 0.0
	kpc_DA = 0.0
	DL = 0.0       # luminosity distance
	DL_Mpc = 0.0
	DL_Gyr = 0.0   # DL in units of billions of light years
	V_Gpc = 0.0
	a = 1.0        # 1/(1+z), the scale factor of the Universe
	az = 0.5       # 1/(1+z(object))
	
	h = H0/100.
	WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
	WK = 1-WM-WR-WV
	az = 1.0/(1+1.0*z)
	age = 0.
	n=1000         # number of points in integrals
	for i in range(n):
		a = az*(i+0.5)/n
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
		age = age + 1./adot
	
	zage = az*age/n
	zage_Gyr = (Tyr/H0)*zage
	DTT = 0.0
	DCMR = 0.0
	
	# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
	for i in range(n):
		a = az+(1-az)*(i+0.5)/n
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
		DTT = DTT + 1./adot
		DCMR = DCMR + 1./(a*adot)
	
	DTT = (1.-az)*DTT/n
	DCMR = (1.-az)*DCMR/n
	age = DTT+zage
	age_Gyr = age*(Tyr/H0)
	DTT_Gyr = (Tyr/H0)*DTT
	DCMR_Gyr = (Tyr/H0)*DCMR
	DCMR_Mpc = (c/H0)*DCMR
	
	# tangential comoving distance
	
	ratio = 1.00
	x = sqrt(abs(WK))*DCMR
	if x > 0.1:
		if WK > 0:
			ratio =  0.5*(exp(x)-exp(-x))/x 
		else:
			ratio = sin(x)/x
	else:
		y = x*x
		if WK < 0: y = -y
		ratio = 1. + y/6. + y*y/120.
	DCMT = ratio*DCMR
	DA = az*DCMT
	DA_Mpc = (c/H0)*DA
	kpc_DA = DA_Mpc/206.264806
	DA_Gyr = (Tyr/H0)*DA
	DL = DA/(az*az)
	DL_Mpc = (c/H0)*DL
	DL_Gyr = (Tyr/H0)*DL
	
	# comoving volume computation
	
	ratio = 1.00
	x = sqrt(abs(WK))*DCMR
	if x > 0.1:
		if WK > 0:
			ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
		else:
			ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
	else:
		y = x*x
		if WK < 0: y = -y
		ratio = 1. + y/5. + (2./105.)*y*y
	VCM = ratio*DCMR*DCMR*DCMR/3.
	V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM
	
	return ['%1.3f' % zage_Gyr, kpc_DA]

									
def run_pyfits(multiFitsFilename):
	# use pyfits to gather info from output of galfit
	imageHeaders = pyfits.open(multiFitsFilename)
	
	# get the dictionary mapping header keywords to their values
	try:
		imageHeader = imageHeaders[1].header
		modelHeader = imageHeaders[2].header
		imageHeaders.close()
	except KeyError:
		print(" method must be called on a multi-extension cube (which galfit outputs)")
		imageHeaders.close()
		exit()
					
	# create list of dictionaries representing component models
	resultModels = []
	compNum = 1
	while ("COMP_"+str(compNum)) in modelHeader:
		model={}
		if modelHeader["COMP_"+str(compNum)].strip() == "sky":
			resultModels.append({"sky": removeGalfitChars(str(modelHeader[str(compNum)+"_SKY"]).split()[0])})
			compNum = compNum + 1
			continue
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
		resultModels.append(model)
		compNum = compNum + 1
			
	return resultModels


def getCentermostID(imageHeight, imageWidth, models):
	'''
	return the integer ID of the galaxy closest to the images center
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
		dx = imageWidth/2.0 - float(px)
		dy = imageHeight/2.0 - float(py)
		dist = sqrt(pow(dx,2) + pow(dy,2))
		# if so, save ID
		if dist < closestDist:
			closestDist = dist
			closestID = currentID
	
	# return the ID of the component closest to the center of the image
	return closestID


def getNextCentermostID(imageHeight, imageWidth, models, centerID):
	'''
	return the integer ID of the galaxy closest to the images center
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
		dx = imageWidth/2.0 - float(px)
		dy = imageHeight/2.0 - float(py)
		dist = sqrt(pow(dx,2) + pow(dy,2))
		# if so, save ID
		if dist < closestDist:
			closestDist = dist
			closestID = currentID
	
	# return the ID of the component closest to the center of the image
	return closestID


def removeGalfitChars(resultString):
	return resultString.replace('*','').replace('[','').replace(']',''
									).replace('{','').replace('}','')
											

def sum_galfit(resultFilename, models, delim, centerIDs, options):
	'''
	returns a string summary of the result specified by the parameter result filename
	
	parameter resultFilename - the result filename from running galfit to be summarized
	return - a string summarizing the results, ending in a new line
	'''
		
	outputFilename = resultFilename.split("/")[-1].strip()
	# VELA02MRP_0.201015_0002949__skipir_CAMERA0-BROADBAND_F160W_simulation_bulge_multi.fits

	#TODO: this is dependent on the particular structure of the filename
	#		the below code works for VELA simulations
	galaxyID = outputFilename.split("_")[0]
	#"VELA01"
	
	filt = outputFilename.split("_")[6]
	#"F125W"	
	
	timeStep = outputFilename.split("_")[1].split(".")[1]
	# "110"
	
	camera = outputFilename.split("_")[5].split("-")[0].replace("CAMERA","")
	# "CAMERA0"
	
	# a = 1/(1+z)
	timeA = float("0."+timeStep)
	# z = 1/a - 1
	timeZ = 1.0/timeA - 1.0
	[age_gyr, kpcPerArcsec] = ned_wright_cosmology_calculator(timeZ)
	age_gyr = str(age_gyr)
	# convert radius from pixels to kpc
	if not options.candelized:
		kpcPerPixel = (119.556)*(1.0/1000.0)#(pc/pixel)(kpc/pc)
	else:
		kpcPerPixel = kpcPerArcsec * 0.06
		
	componentList = [galaxyID, timeStep, age_gyr, str(timeZ), camera, filt]
	
	componentResults = ""
	sky = "0.0"
	for currentID, model in enumerate(models):
		
		# skip sky
		if "sky" in model:
			sky = str(model["sky"])
			continue
		
		# to correct for indexing starting at zero
		currentID = currentID + 1
		
		#position
		componentList.append(model["px"][0])
		componentList.append(model["py"][0])
		
		# magnitude
		componentList.append(model["mag"][0])

		# radius
		componentList.append(model["rad"][0])
		componentList.append(str(float(model["rad"][0])*kpcPerPixel))
			
		# sersic index
		componentList.append(model["ser"][0])
		if currentID in centerIDs:
			if options.bulge:
				if componentList[-1] == '1.0000':
					galaxyType = "disk"
				else:
					galaxyType = "bulge"
			else:
				galaxyType = "central"
		else:
			galaxyType = "other"
			
		# b/a
		componentList.append(model["ba"][0])
			
		# position angle
		componentList.append(model["pa"][0])
				
		# add this completed component as a line to the string of results
		componentResults = (componentResults + 
						delim.join([galaxyType]+componentList) + "\n")
		
		# reset component list for any remaining components in this file
		componentList = [galaxyID, timeStep, age_gyr, str(timeZ), camera, filt]
		
	results = ""
	for component in componentResults.split("\n")[:-1]:
		results = results + component + delim + sky + "\n"
	return results


if __name__ == "__main__":
	
	#define the command line interface with simUtility.py
	usage = ("\n%prog resultsFile [-h help] [options (with '-'|'--' prefix)]")

	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-b","--bulge", 
				help="if running on results of GALFIT bulge (two component) fit",
				action="store_true")

	# indicate that images are candelized
	parser.add_option("-c","--candelized", 
                      help="if running on candelized images",
                      action="store_true")
	
	parser.add_option("-o","--output", 
				help="set the filename to write the output summary file",
				default="summary_" + time.strftime("%m-%d-%Y") + ".txt")
	
	# parse the command line using above parameter rules
	# options - list with everthing defined above, 
	# args - anything left over after parsing options
	[options, args] = parser.parse_args()
	pprint.pprint(vars(options))

	if not (len(args) or os.path.isfile(args[0])):
		parser.error("results file must be an existing file with result filenames")
	else:
		inputFilename = args[0]
		
	#this will be the file that will contain the images
	r = open(inputFilename, 'r')
	
	# the lines of file containing results
	resultFilenames = r.readlines()
	
	# close the file now that it has been read
	r.close()
	
	outFile = open(options.output, 'w')
	
	delim = " "
	
	outFile.write("galfit result file run on " + time.strftime("%m-%d-%Y") + "\n" + 
				delim.join(["type","galaxyID","timeStep","age(GYr)","redshift(z)",
							"camera","filter","px(pixels)","py(pixels)","mag",
							"rad(pixels)","rad(kpc)", 
							"sersicIndex(n)","b/a","angle(deg)","sky"]) + "\n")
	
	# this loops through every result, writing summary to output
	for resultFilename in resultFilenames:
	
		models = run_pyfits(resultFilename)	
		
		# get the id of the centermost galaxy
		centerID = getCentermostID(600,600,models)
		nextCenterID = getNextCentermostID(600,600,models,centerID)
	
		# summarize galfit and write to output
		if options.bulge:
			centerIDs = [centerID, nextCenterID]
		else:
			centerIDs = [centerID]
		outFile.write( sum_galfit(resultFilename.strip(), models, delim, centerIDs, options) )
		
	outFile.close()
	
	print ("Output written to " + options.output)
######################### done ################################################