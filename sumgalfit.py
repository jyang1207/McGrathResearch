
#! /usr/bin/python


import os
import time
import sys
from math import sqrt, exp, sin, pi
from optparse import OptionParser


def compute_distance(p0, p1):
	'''
	http://stackoverflow.com/questions/5407969/distance-formula-between-two-points-in-a-list
	'''
	
	return sqrt((float(p0[0]) - float(p1[0]))**2 + (float(p0[1]) - float(p1[1]))**2)
	
	
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


def sum_galfit(resultFilename, maxDist, minSersIndex, maxSersIndex):
	'''
	returns a string summary of the result specified by the parameter result filename
	
	parameter resultFilename - the result filename from running galfit to be summarized
	return - a string summarizing the results, ending in a new line
	'''
	
	delim = " "
	
	resultFile = open(resultFilename, 'r')
	
	resultLines = resultFile.readlines()
	
	resultFile.close()
	
	skipSky = False
	galaxyID = ""
	px1 = ""
	px2 = ""
	mag1 = ""
	mag2 = ""
	rad1 = ""
	rad2 = ""
	sersIndex1 = ""
	sersIndex2 = ""
	ba1 = ""
	ba2 = ""
	pa1 = ""
	pa2 = ""
	
	for resultLine in resultLines:
		
		if not galaxyID and resultLine.strip()[:2] == "B)":
		
			# B) a0.220/VELA01_220_cam0_F160W_multi.fits      # Output data image block
			galaxyID = resultLine.split("/")[-1].split("_")[0]
			timeStep = resultLine.split("/")[-1].split("_")[1]
			camera = resultLine.split("/")[-1].split("_")[2]
			filt = resultLine.split("/")[-1].split("_")[3]
	
					
		if resultLine.strip()[:6] == "0) sky":
			skipSky = True
		elif resultLine.strip()[:6] == "0) ser":
			skipSky = False
			
	
		if not skipSky:
			if not px1 and resultLine.strip()[:2] == "1)":
			
				#1) 301.6210 299.7872 1 1  #  Position x, y
				px1 = resultLine.strip().split(" ")[1]
				py1 = resultLine.strip().split(" ")[2]
	
			elif not px2 and resultLine.strip()[:2] == "1)":
				px2 = resultLine.strip().split(" ")[1]
				py2 = resultLine.strip().split(" ")[2]
	
			elif not mag1 and resultLine.strip()[:2] == "3)":
				mag1 = resultLine.strip().split(" ")[1]
				
			elif not mag2 and resultLine.strip()[:2] == "3)":
				mag2 = resultLine.strip().split(" ")[1]
	
			elif not rad1 and resultLine.strip()[:2] == "4)":
				rad1 = resultLine.strip().split(" ")[1]
				
			elif not rad2 and resultLine.strip()[:2] == "4)":
				rad2 = resultLine.strip().split(" ")[1]
				
			elif not sersIndex1 and resultLine.strip()[:2] == "5)":
				sersIndex1 = resultLine.strip().split(" ")[1]
				
			elif not sersIndex2 and resultLine.strip()[:2] == "5)":
				sersIndex2 = resultLine.strip().split(" ")[1]
				
			elif not ba1 and resultLine.strip()[:2] == "9)":
				ba1 = resultLine.strip().split(" ")[1]
				
			elif not ba2 and resultLine.strip()[:2] == "9)":
				ba2 = resultLine.strip().split(" ")[1]
				
			elif not pa1 and resultLine.strip()[:3] == "10)":
				pa1 = resultLine.strip().split(" ")[1]
				
			elif not pa2 and resultLine.strip()[:3] == "10)":
				pa2 = resultLine.strip().split(" ")[1]
	
	# a = 1/(1+z)
	# z = 1/a - 1
	age_gyr = str(red_shift_to_gyr(1/(float(timeStep)/1000) - 1))
	
	errorFlag1 = ""
	errorFlag2 = ""
	posFlag = ""
	
	# sersic index of 0.5 < si < 10 is good
	if float(sersIndex1) < minSersIndex or float(sersIndex1) > maxSersIndex:
		sersFlag1 = "*"
		errorFlag1 = "*"
	else:
		sersFlag1 = ""

	
	if not px2:
		if float(sersIndex1) > 2.5:
			type0 = "bulge"
		else:
			type0 = "disk"
			
		return delim.join(errorFlag1 + galaxyID, timeStep, age_gyr, 
				camera, filt, type0, px1, py1, 
				sersIndex1, mag1, rad1, ba1, pa1, "0\n")
				
	
	# component seperation of greater than 5 is bad 
	dist = compute_distance([px1, py1], [px2, py2])
	if dist > maxDist:
		posFlag = "*"
		errorFlag1 = "*"
		errorFlag2 = "*"
		
		
	# sersic index of 0.5 < si < 10 is good
	if float(sersIndex2) < minSersIndex or float(sersIndex2) > maxSersIndex:
		sersFlag2 = "*"
		errorFlag2 = "*"
	else:
		sersFlag2 = ""
		
	# test for type (bulge or disk)
	if float(sersIndex1) > 2.5 and float(sersIndex2) < 2.5:
		type1 = "bulge"
		type2 = "disk"
	elif float(sersIndex1) < 2.5 and float(sersIndex2) > 2.5:
		type2 = "bulge"
		type1 = "disk"
	elif float(sersIndex1) > 2.5 and float(sersIndex2) > 2.5:
		type1 = "bulge"
		type2 = "bulge"
	else:
		if rad1 < rad2:
			type1 = "bulge"
			type2 = "disk"
		else:
			type2 = "bulge"
			type1 = "disk"
			
	result1 = delim.join(errorFlag1 + galaxyID, timeStep, age_gyr, 
				camera, filt, type1, px1, py1, 
				sersIndex1, mag1, rad1, ba1, pa1, str(dist) + "\n")
	result2 = delim.join(errorFlag2 + galaxyID, timeStep, age_gyr, 
				camera, filt, type2, px2, py2, 
				sersIndex2, mag2, rad2, ba2, pa2, str(dist) + "\n")
	
	# just to get rid of unused variable warnings, these flags are not used anymore
	if not (sersFlag1 + sersFlag2 + posFlag):
		continue
	
	return result1 + result2


if __name__ == "__main__":
	
	usage = "usage: %prog [-h for help] [options] input output" 

	# used to parse command line arguments
	parser = OptionParser(usage)
						
	# file specifies the full path filename of the list of images to run
	parser.add_option("--maxDistance", 
						help="set the threshold for distance between multiple components of the same image, beyond which a * will indicate the error",
						type=float, default=5.0)
						
	# file specifies the full path filename of the list of images to run
	parser.add_option("--minSersicIndex", 
						help="set the min threshold for sersic index, below which a * will indicate the error",
						type=float, default=0.5)
						
	# file specifies the full path filename of the list of images to run
	parser.add_option("--maxSersicIndex", 
						help="set the max threshold for sersic index, above which a * will indicate the error",
						type=float, default=10.0)
						
	# Magnitude photometric zeropoint					
	# Plate scale
	# PSF
	
	# parse the command line using above parameter rules
	[options, args] = parser.parse_args()
		
	if len(args) != 2:
		parser.error("incorrect number of arguments, " + 
					"must provide an input file followed by an output file")
	elif not os.path.isfile(args[0]):
		parser.error("input file " + args[0] + " either does not exist or is not accessible")
		
	#this will be the file that will contain the images
	r = open(args[0], 'r')		
	
	# the lines of file containing results
	resultFilenames = r.readlines()
	
	# close the file now that it has been read
	r.close()
		
	outFile = open(args[1], 'w')
	
	outFile.write("galfit result file run on " + time.strftime("%m-%d-%Y") + " with " + 
				"max distance = " + str(options.maxDistance) + ", " + 
				"min sersic index = " + str(options.minSersicIndex) + ", and " + 
				"max sersic index = " + str(options.maxSersicIndex) + "\n" + 
				"galaxyID timeStep age(GYr) camera filter type " + 
				"px py sersicIndex mag rad b/a angle" + 
				", component separation distance\n")
	
	# this loops through every result, writing summary to output
	for resultFilename in resultFilenames:
	
		try:
			# summarize galfit and write to output
			outFile.write( sum_galfit(resultFilename.strip(), 
										options.maxDistance, 
										options.minSersicIndex, 
										options.maxSersicIndex))
		except:
			print (str(sys.exc_info()[0]) + str(sys.exc_info()[1]))
		
	outFile.close()
	
	print ("Done!")
######################### done ################################################