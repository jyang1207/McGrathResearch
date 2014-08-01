#!/usr/bin/env python


import os
import time
import sys
from math import sqrt, exp, sin, pi

	
def red_shift_to_gyr(z):
	'''
	TODO: paste url of site where this method came from
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
	
	return '%1.3f' % zage_Gyr


def sum_galfit(resultFilename, delim):
	'''
	returns a string summary of the result specified by the parameter result filename
	
	parameter resultFilename - the result filename from running galfit to be summarized
	return - a string summarizing the results, ending in a new line
	'''
	
	resultFile = open(resultFilename, 'r')
	
	resultLines = resultFile.readlines()
	
	resultFile.close()
		
	skipSky = False
	componentResults = ""
	for resultLine in resultLines:
		
		# use output filename to gather some general info
		if resultLine.strip()[:2] == "B)":
			
			fullPathOutput = resultLine.split()[-1]
			# results/VELA02.....
			
			outputFilename = fullPathOutput.split("/")[-1].strip()
			# VELA02MRP_0.201015_0002949__skipir_CAMERA0-BROADBAND_F160W_simulation_bulge_multi.fits

			#TODO: this is dependent on the particular structure of the filename
			#		the below code works for VELA simulations
			galaxyID = outputFilename.split("_")[0]
			#"VELA01"
			
			filt = outputFilename.split("_")[6]
			#"F125W"	
			
			timeStep = outputFilename.split("_")[1].split(".")[1]
			# "110"
			
			cam_str = outputFilename.split("_")[5].split("-")[0]
			# "CAMERA0"
			
			# these statements accomodate for camera numbers
			camera = ""
			for c in cam_str:
				if c.isdigit():
					camera = camera + c
			
			# a = 1/(1+z)
			# z = 1/a - 1
			age_gyr = str(red_shift_to_gyr(1/float("0."+timeStep) - 1))
	
			componentList = [galaxyID, timeStep, age_gyr, camera, filt]
			
		# toggle flag on the line indicating sky for use by future lines
		if resultLine.strip()[:6] == "A) sky":
			skipSky = True
		elif resultLine.strip()[:6] == "B) ser":
			skipSky = False
			
		# dont record results for sky component
		if not skipSky:
		
			#position
			if resultLine.strip()[:2] == "1)":
			
				#1) 301.6210 299.7872 1 1  #  Position x, y
				componentList.append(resultLine.strip().split()[1])
				componentList.append(resultLine.strip().split()[2])
			
			# magnitude
			elif resultLine.strip()[:2] == "3)":
				componentList.append(resultLine.strip().split()[1])
	
			# radius
			elif resultLine.strip()[:2] == "4)":
				componentList.append(resultLine.strip().split()[1])
				
			# sersic index
			elif resultLine.strip()[:2] == "5)":
				componentList.append(resultLine.strip().split()[1])
				
			# b/a
			elif resultLine.strip()[:2] == "9)":
				componentList.append(resultLine.strip().split()[1])
				
			# position angle
			elif resultLine.strip()[:3] == "10)":
				componentList.append(resultLine.strip().split()[1])
				
				# add this completed component as a line to the string of results
				componentResults = componentResults + errorFlag + delim.join(componentList) + "\n"
				
				# reset component list for any remaining components in this file
				componentList = [galaxyID, timeStep, age_gyr, camera, filt]
		
	return componentResults


if __name__ == "__main__":

	if len(sys.argv) == 2:
		inputFilename = sys.argv[1]
		outputFilename = "summary_" + time.strftime("%m-%d-%Y") + ".txt"
	elif len(sys.argv) == 3:
		inputFilename = sys.argv[1]
		outputFilename = sys.argv[2]
	else:
		print ("Usage: %prog input output\nincorrect number of arguments, " + 
			"must provide an input filename followed by a desired output filename")
		exit()
		
	if not os.path.isfile(inputFilename):
		print ("Usage: %prog input output\ninput file " + sys.argv[0] + 
				" either does not exist or is not accessible")
		
	#this will be the file that will contain the images
	r = open(inputFilename, 'r')
	
	# the lines of file containing results
	resultFilenames = r.readlines()
	
	# close the file now that it has been read
	r.close()
	
	outFile = open(outputFilename, 'w')
	
	delim = " "
	
	outFile.write("galfit result file run on " + time.strftime("%m-%d-%Y") + "\n" + 
				delim.join(["galaxyID", "timeStep", "age(GYr)", "camera", 
							"filter", "px", "py", "mag", "rad", 
							"sersicIndex", "b/a", "angle"]) + "\n")
	
	# this loops through every result, writing summary to output
	for resultFilename in resultFilenames:
	
		# summarize galfit and write to output
		outFile.write( sum_galfit(resultFilename.strip(), delim) )
		
	outFile.close()
	
	print ("Output written to " + outputFilename)
######################### done ################################################