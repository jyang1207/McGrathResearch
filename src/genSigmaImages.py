#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 3/9/2015
Colby College Astrophysics Research
'''

import sys
import os
import numpy as np
try:
	import pyfits as fits
	print("using pyfits")
except ImportError:
	from astropy.io import fits
	print("using astropy")

def genSigmaImages(simFilename):
	'''
	update the given simulation .fits file with the header info needed
	from its candelized counterpart.
	
	parameter simFilename - the filename of the simulation .fits image
	returns - boolean indicating if an error occurs while trying to update
	'''
	result = []
	prefix = "_".join(simFilename.split("_")[:-1])
	sigmaFilename = prefix + "_simulation_sigma.fits"
	
	# attempt to open the sim file using fits and get the data
	try:
		simHDUList = fits.open(simFilename)
		simData = simHDUList[0].data
		# for simulations, sigma is sqrt(simData)
		data = np.power(simData, 0.5)
		fits.writeto(sigmaFilename, data, clobber=True)
		simHDUList.close()
		result.append(sigmaFilename + " written")
	except:
		result.append(simFilename + " failed")
		
	# use the sim filename to get the candelized image filename
	candelNoNoiseFilename = prefix + "_candelized_nonoise.fits"
	candelNoiseFilename = prefix + "_candelized_noise.fits"
	sigmaNoNoiseFilename = prefix + "_candelized_nonoise_sigma.fits"
	sigmaNoiseFilename = prefix + "_candelized_noise_sigma.fits"
	try:
		candelNoNoiseHDUList = fits.open(candelNoNoiseFilename)
		candelNoNoiseData = candelNoNoiseHDUList[0].data
		candelNoiseHDUList = fits.open(candelNoiseFilename)
		candelNoiseData = candelNoiseHDUList[0].data
		# for candelized, sigma is sqrt(rms^2 + (sqrt(pixelVal*expTime)/expTime)^2)
		rms = candelNoiseData - candelNoNoiseData
		expTime = 3000.0 # TODO: adjust for your data
		data = np.sqrt(np.square(rms) + 
					np.square(np.sqrt(np.abs(candelNoiseData * expTime)) / expTime))
		fits.writeto(sigmaNoiseFilename, data, clobber=True)
		data = np.sqrt(candelNoNoiseData) / expTime
		fits.writeto(sigmaNoNoiseFilename, data, clobber=True)
		candelNoiseHDUList.close()
		candelNoNoiseHDUList.close()
		result.append(sigmaNoiseFilename + " written")
		result.append(sigmaNoNoiseFilename + " written")
	except:
		result.append(candelNoNoiseFilename + " failed")
		
	return result

if __name__ == "__main__":

	# define the command line interface with simUtility.py
	usage = ("USAGE: python genSigmaImages.py " + 
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
		for result in genSigmaImages(simFilename.strip()):
			print(result)
