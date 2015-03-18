#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 3/9/2015
Colby College Astrophysics Research
'''

import sys
import os
try:
	import pyfits as fits
	print("using pyfits")
except ImportError:
	from astropy.io import fits
	print("using astropy")

def updateHeader(simFilename):
	'''
	update the given simulation .fits file with the header info needed
	from its candelized counterpart.
	
	parameter simFilename - the filename of the simulation .fits image
	returns - boolean indicating if an error occurs while trying to update
	'''
	
	# attempt to open the sim file using fits and get the header dictionary
	try:
		simHDUList = fits.open(simFilename)
		simHeader = simHDUList[0].header
	except:
		print("failed to open sim file " + simFilename + 
				" with fits")
		return False
		
	# use the sim filename to get the candelized image filename
	candelFilename = "_".join(simFilename.split("_")[:-1]) + "_candelized_nonoise.fits"
	
	# attempt to open the candelized file using fits and get the header dictionary
	try:
		candelHDUList = fits.open(candelFilename)
		candelHeader = candelHDUList[0].header
	except:
		print("failed to open candelized file " + candelFilename + 
				" with fits")
		simHDUList.close()
		return False
	
	try:
		simHeader.set("SCALESIM", candelHeader["SCALESIM"], "pc per pixel of SUNRISE img")
	except:
		print("failed to write to sim header " + simFilename + 
				" from candel header " + candelFilename + 
				" with fits")
		simHDUList.close()
		candelHDUList.close()
		return False
		
	simHDUList.writeto(simFilename, clobber=True)
	simHDUList.close()
	candelHDUList.close()
	return True

if __name__ == "__main__":

	# define the command line interface with simUtility.py
	usage = ("USAGE: python updateHeader.py " + 
			"<file containing list of simulation .fits filenames to be updated>")

	# must have exactly one positional command line argument
	if len(sys.argv) != 2 or not os.path.isfile(sys.argv[1]):
		print(usage)
		exit()
		
	# read all simulation .fits filenames from command line file
	with open(sys.argv[1], 'r') as simListFile:
		simFilenames = simListFile.readlines()
	
	# iterate over filenames, updating the header info
	for simFilename in simFilenames:
		if not updateHeader(simFilename.strip()):
			continue
