
#! /usr/bin/python


import os
import sys
import re
from pyraf import iraf
import imexam					#works
import dimensions
import readin
import math


def compute_distance(p0, p1):
	'''
	http://stackoverflow.com/questions/5407969/distance-formula-between-two-points-in-a-list
	'''
	return math.sqrt((float(p0[0]) - float(p1[0]))**2 + (float(p0[1]) - float(p1[1]))**2)


def init_galfit_parameter_files():
	'''
	creates and writes header to data lower and upper
	'''
	
	# create a txt file for the result data of the lower half of the picture in this directory
	os.system("touch " + "data_lower.txt")
	
	# open the result stroing file in write mode
	format = open("data_lower.txt", "w")
	
	# write a header line to define what will go in the result file
	format.write("Input Image" + "	" + "Time Step" + "	" + "X Value" + 
					"	" + "X Error" + "	" + "Y Value" + "	" + "Y Error" + 
					"	" + "Magnitude" + "	" + "Magnitude Error" + "	" + "Radius" + 
					"	" + "Radius Error" + "	" + "Sersic" + "	" + "Sersic Error" + 
					"	" + "BA Value" + "	" + "BA Error" + "	" + "Position Angle" + 
					"	" + "Position Angle Error" + "	" + "Chi Squared" + "	" + "ndof" + 
					"	" + "Chi Squared/nu" + "\n")
	
	# close the result file
	format.close()
	
	# create a txt file for the result data of the upper half of the picture in this directory
	os.system("touch " + "data_upper.txt")
	
	# open the result stroing file in write mode
	format = open("data_upper.txt", "w")
	
	# write a header line to define what will go in the result file
	format.write("Input Image" + "	" + "Time Step" + "	" + "X Value" + 
					"	" + "X Error" + "	" + "Y Value" + "	" + "Y Error" + 
					"	" + "Magnitude" + "	" + "Magnitude Error" + "	" + "Radius" + 
					"	" + "Radius Error" + "	" + "Sersic" + "	" + "Sersic Error" + 
					"	" + "BA Value" + "	" + "BA Error" + "	" + "Position Angle" + 
					"	" + "Position Angle Error" + "	" + "Chi Squared" + "	" + "ndof" + 
					"	" + "Chi Squared/nu" + "\n")
	
	# close the result file
	format.close()


def run_imhead(imageFilename):
	'''
	returns a list of the form [cam_number, image_number, frame_height, frame_width]
	'''
	
	imhead_return = iraf.imhead(imageFilename, Stdout=1)[0].strip()
	# "broadband_200_cam0_g.fits[1893,1893][real]:" 
	
	frame_dimensions = imhead_return.split("[")[1].replace("]", "").split(",")
	# ["1893","1893"]

	frame_width = frame_dimensions[0]
	# "1893"
	
	frame_height = frame_dimensions[1]
	# "1893"
	
	image_number = imhead_return.split("_")[1]
	# "200"
	
	cam_str = imhead_return.split("_")[2]
	# "cam0"
	
	#these statements accomodate for camera numbers greater than 9
	conditions = ["1","2","3","4","5","6","7","8","9"]			
	
	#determines if camera number is greater than 10, does not allow greater than 99
	# ex: "cam0" -> M6[-2] = "m"
	# ex: "cam12" -> M6[-2] = "1"
	if cam_str[-2] in conditions:										
		cam_number = cam_str[-2:]
	else:
		cam_number = cam_str[-1]					
	#print cam_number
	
	return [cam_number, image_number, frame_height, frame_width]
	
				
def run_imexam(image):
	'''
	returns list of info in the form [line, col, mag, radius, b_a, PA]
	'''
	
	#writes the output of the minmax function of iraf to a file for later use. 
	os.system('touch'+' coords.tmp')
	write_coords = open('coords.tmp', 'w')
	
	#call iraf's minmax function
	maxCoords = iraf.minmax(image, Stdout=1, update=0
		)[0].strip().split(" ")[3].replace("[", "").replace("]", "").split(",")
	# "broadband_200_cam0_g.fits[1:1893,1:946] [1765,691] -0.2845799624919891 [1283,692] 19.96048927307129"
	# ["1283","692"]
	write_coords.write(maxCoords[0]) 				
	write_coords.write(" " + maxCoords[1])
	write_coords.close()	
	
	#run imexam passing coords.tmp
	imexam_out = str(iraf.imexam(image, use_display=0, imagecur="coords.tmp", Stdout=1))
	# "display frame (1:) (1): ['#   COL    LINE    COORDINATES', '#     " +
	# "R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT', " +
	# "'1283.40  692.16 1283.40 692.16', '  63.22  14.53  15461.  0.1202   19.53 0.52    " +
	# "1 1.75    21.58    18.85  21.39']"	
	
	# delete coords.tmp, no longer needed
	os.system('rm' + ' coords.tmp')	
	
	imexam_array = str.split(imexam_out, "'")
	# ["display frame (1:) (1): [", 
	# "#   COL    LINE    COORDINATES", 
	# ", ", "#     " +
	# "R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT", ", ",
	# "1283.40  692.16 1283.40 692.16", ", ", "  63.22  14.53  15461.  0.1202   19.53 0.52    " +
	# "1 1.75    21.58    18.85  21.39']"			
	
	#print imexam_out
	#print imexam_array	
	
	xy = imexam_array[5].strip()
	data = imexam_array[7].strip()
	data = str.split(data)
	xy = str.split(xy)
	
	col = xy[0] 									#this gives x value	
	line = xy[1]									#this gives y value
	mag = data[1]									#this prints the magnitude
	radius = data[0]								#this prints the radius
	b_a = int(1) - float(data[5])					#this prints the E variable that relates to b/a. note b/a=1-e
	PA = float(data[6]) - int(90)					#this gives the position angle. iraf measures up from x. Galfit down from y

# 	print col
# 	print line
# 	print mag
# 	print radius
# 	print b_a
# 	print PA
	#exit()	
	
	return [line, col, mag, radius, b_a, PA]
	

def read_param(logfile, distance):
	'''
	logfile is the filename of galfits log file (e.g. fit.log)
	distance is an integer
	'''

	lines = logfile.readlines()
	

	line1 = lines[7]
	line2 = lines[8]
	line3 = lines[9]
	line4 = lines[10]
	line5 = lines[13]
	line6 = lines[14]
	line7 = lines[11]
	line8 = lines[12]
	infoline = lines[2]
	#print line1
	#print line2
	#print line3
	#print line4
	#print line5
	#print line6
	#print line7
	#print line8
	
	imname = str.split(infoline, " ")
	imagename = imname[7]
	time = str.split(imagename, "_")
	timestep = time[1]
	
	dist_min = 20	#need to change this value in both rungalfit and readin if you wish to modify it
	#print distance
	
	split = line1.rstrip()
	split1 = " ".join(split.split())
	split2 = split1.replace("(", "")
	split3 = split2.replace(")", "")
	split4 = split3.replace(",", "")
	split5 = str.split(split4, " ")
	split6 = split5[3:]
	#print split6
	[X_fin,Y_fin,mag_fin,rad_fin,sersic_fin,BA_fin, PA_fin] = split6				

	split7 = line2.rstrip()
	split8 = " ".join(split7.split())
	split9 = split8.replace("(", "")
	split10 = split9.replace(")", "")
	split11 = split10.replace(",", "")
	split12 = str.split(split11, " ")
	split13 = split12[1:]
	#print split13
	[X_err,Y_err,mag_err,rad_err,sersic_err,BA_err, PA_err] = split13
	#print X_err
	
	if distance > dist_min: 
		split14 = line3.rstrip()
		split15 = " ".join(split14.split())
		split16 = split15.replace("(", "")
		split17 = split16.replace(")", "")
		split18 = split17.replace(",", "")
		split19 = str.split(split18, " ")
		split20 = split19[3:]
		#print split20
		[X_fin2,Y_fin2,mag_fin2,rad_fin2,sersic_fin2,BA_fin2, PA_fin2] = split20			#note 2 refers to top component
		#print BA_fin2

		split21 = line4.rstrip()
		split22 = " ".join(split21.split())
		split23 = split22.replace("(", "")
		split24 = split23.replace(")", "")
		split25 = split24.replace(",", "")
		split26 = str.split(split25, " ")
		split27 = split26[1:]
		#print split27
		[X_err2,Y_err2,mag_err2,rad_err2,sersic_err2,BA_err2, PA_err2] = split27
		#print X_err 
	
		split28 = line5.rstrip()
		split29 = split28.replace(",", "")
		split30 = str.split(split29, " ")
		chi_square = split30[3]
		ndof = split30[7]
		#print chi_square
		#print ndof

		split31 = line6.rstrip()
		split32 = split31.replace(",", "")
		split33 = str.split(split32, " ")
		chi_square_nu = split33[3]
	else: 
		split34 = line7.rstrip()
		split35 = split34.replace(",", "")
		split36 = str.split(split35, " ")
		chi_square = split36[3]
		ndof = split36[7]
		print chi_square
		print ndof

		split37 = line8.rstrip()
		split38 = split37.replace(",", "")
		split39 = str.split(split38, " ")
		chi_square_nu = split39[3]
		print chi_square_nu
	
	
	data_lower = open("data_lower.txt", "a")
	data_upper = open("data_upper.txt", "a")
	#open_dist = open("dist.tmp", "r")
	#d = open_dist.readlines()
	#open_dist.close()
	#distance = int(d[0])
	#print distance
	
	if distance <= dist_min: 
		data_lower.write(imagename + "*" + "	" + timestep + "	" + X_fin + "	" + X_err + "	" + Y_fin + "	" + Y_err + "	" + mag_fin + "	" + mag_err + "	" + rad_fin + "	" + rad_err + "	" + sersic_fin + "	" + sersic_err + "	" + BA_fin + "	" + BA_err + "	" + PA_fin + "	" + PA_err + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_lower.close()
		data_upper.write(imagename + "*" + "	" + timestep + "	" + X_fin + "	" + X_err + "	" + Y_fin + "	" + Y_err + "	" + mag_fin + "	" + mag_err + "	" + rad_fin + "	" + rad_err + "	" + sersic_fin + "	" + sersic_err + "	" + BA_fin + "	" + BA_err + "	" + PA_fin + "	" + PA_err + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_upper.close()
	else: 
		data_lower.write(imagename + "	" + timestep + "	" + X_fin + "	" + X_err + "	" + Y_fin + "	" + Y_err + "	" + mag_fin + "	" + mag_err + "	" + rad_fin + "	" + rad_err + "	" + sersic_fin + "	" + sersic_err + "	" + BA_fin + "	" + BA_err + "	" + PA_fin + "	" + PA_err + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_lower.close()
		data_upper.write(imagename + "	" + timestep + "	" + X_fin2 + "	" + X_err2 + "	" + Y_fin2 + "	" + Y_err2 + "	" + mag_fin2 + "	" + mag_err2 + "	" + rad_fin2 + "	" + rad_err2 + "	" + sersic_fin2 + "	" + sersic_err2 + "	" + BA_fin2 + "	" + BA_err2 + "	" + PA_fin2 + "	" + PA_err2 + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_upper.close()
	logfile.close()
	
	
def run_galfit(imageListFilename):
	'''
	opens file (parameter) containing list of images and loops over every image, running galfit
	'''

	#this will be the file that will contain the images
	f = open(imageListFilename)		
	
	# the first line of file containing images
	imageFilenames = f.readlines()
	
	# close the file now that it has been read
	f.close()
	
	# this loops through every image in images file and
	for imageFilename in imageFilenames:
		
		# remove the "\n" on the end of filename
		imageFilename = imageFilename.strip()
		
		# run the method of dimensions.py
		[camera, image, height, width] = run_imhead(imageFilename)
		
		# calls a method of the imexam.py file, passing filename (with width and height) as a parameter
		# imexam gets the coordinates of the max pixel and uses iraf.imexam to set global
		# variables detailing the results of iraf's examination of the max coordinate
		[Y, X, magnitude, rad, BA, angle] = run_imexam(
			imageFilename + '[1:' + width + ',1:' + str(int(height)/2) + ']')
		
		# same as above, but for the 'top' of the image
		[Y_top, X_top, magnitude_top, rad_top, BA_top, angle_top] = run_imexam(
			imageFilename + '[1:' + width + ',' + str((int(height)/2)+1) + ':' + height + ']')
		
		distance = int(compute_distance([X, Y], [X_top, Y_top]))
		print distance
		print str(distance)
		
		#adjustable parameter for when to run galfit on second component	
		if distance <= 20:																							
			Z = "# "
		else:
			Z = " "
	
		#write the galfit parameter file
		galfit_parameter_filename = 'bb' + image + '_c' + str(camera) + '.txt'
		galfit_output_filename = 'bb' + image + '_c' + str(camera) + "_multi.fits"
		galfit_constraint_filename = 'bb' + image + '_c' + str(camera) + "_constraint.txt"
		os.system('touch ' + galfit_parameter_filename)								
	
	
######### writes gathered galfit parameters to file ###########################
		WP = open(galfit_parameter_filename,'w')
		WP.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
		WP.write("A) " + imageFilename + 
				"						#Input data image block\n")
		WP.write("B) " + galfit_output_filename + 
				"						#Output data image block\n")
		WP.write("C)" + " none" + 
				"						#Sigma image name (made from data if blank or 'none')\n")
		WP.write("D)" + " ps4.fits" + 
				"						#Input PSF image and (optional) diffusion kernel\n")
		WP.write("E)" + " 1" + 
				"						#PSF fine sampling factor relative to data\n")
		WP.write("F)" + " none" + 							
				"						#Bad pixel mask (FITS file or ASCIIcoord list)\n")
		WP.write("#G) " + galfit_constraint_filename + 		
				"						#File with parameter constraints (ASCII file)\n")
		WP.write("H)" + " 1	" + width + " 1	" + height + 
				"						#Image region to fit (xmin xmax ymin ymax)\n")
		WP.write("I)" + " 200 200" + 
				"						#Size of the concolution box (x y)\n")
		WP.write("J)" + " 25.0" + 
				"						#Magnitude photometric zeropoint\n")			#this will be user inputed 
		WP.write("K)" + " 0.038" + "  0.038" + 
				"						#Plate scale (dx dy)  [arcsec per pixel]\n")
		WP.write("O)" + " regular" + 
				"						#display type (regular, curses, both\n")
		WP.write("P)" + " 0" + 
				"						#Options: 0=normal run; 1,2=make model/imgblock & quit\n")
		WP.write("S)" + " 0" + 
				"						#Modify/create components interactively?\n")
		WP.write("\n")
		WP.write("# INITIAL FITTING PARAMETERS\n")
		WP.write("#\n")
		WP.write("#	For component type, the allowed functions are:\n")
		WP.write("#		nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n")
		WP.write("#		ferrer, coresersic, sky and isophote.\n")
		WP.write("#\n")
		WP.write("#	Hidden parameters will only appear when they are specified:\n")
		WP.write("#		C0 (diskyness/boxyness),\n")
		WP.write("#		Fn (n=interger, Azimuthal Fourier Modes).\n")
		WP.write("#		R0-R10 (PA rotation, for creating spiral structures).\n")
		WP.write("#\n")
		WP.write("# ------------------------------------------------------------------------------\n")
		WP.write("#		par)	par value(s)	fit toggle(s)	# parameter description\n")
		WP.write("# ------------------------------------------------------------------------------\n")
		WP.write("\n")
		WP.write("# Componenet number: 1\n")
		WP.write(" 0) sersic					#Component type\n")
		WP.write(" 1) " + str(X) + "	" + str(Y) + "	1	1			#Position x,y\n")
		WP.write(" 3) " + str(magnitude) + "	1			#Integrated Magnitude\n")
		WP.write(" 4) " + str(rad) + "			1			#R_e (half-light radius)	[pix]\n")
		WP.write(" 5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
		WP.write(" 6) 0.0000		0			#	-----\n")
		WP.write(" 7) 0.0000		0			#	-----\n")
		WP.write(" 8) 0.0000		0			#	-----\n")
		WP.write(" 9) " + str(BA) + "			1			#Axis ratio (b/a)\n")
		WP.write(" 10) " + str(angle) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
		WP.write(" Z) 0							#Leave in [1] or subtract [0] this comp from data?\n")
		WP.write("\n")
	
		WP.write("# Componenet number: 2\n")
		WP.write(Z + "0) sersic					#Component type\n")
		WP.write(Z + "1) " + str(X_top) + "	" + str(Y_top) + "	1	1			#Position x,y\n")
		WP.write(Z + "3) " + str(magnitude_top) + "	1			#Integrated Magnitude\n")
		WP.write(Z + "4) " + str(rad_top) + "		1			#R_e (half-light radius)	[pix]\n")
		WP.write(Z + "5) " + "1.0000		1			#Sersic index n (de Vaucouleurs n=4)\n")
		WP.write(Z + "6) 0.0000		0			#	-----\n")
		WP.write(Z + "7) 0.0000		0			#	-----\n")
		WP.write(Z + "8) 0.0000		0			#	-----\n")
		WP.write(Z + "9) " + str(BA_top) + "			1			#Axis ratio (b/a)\n")
		WP.write(Z + "10) " + str(angle_top) + "		1			#Position angle (PA) [deg: Up=0, left=90]\n")
		WP.write(Z + "Z) 0						#Leave in [1] or subtract [0] this comp from data?\n")
		WP.write("\n")
	
		WP.write("# Componenet number: 3\n")
		WP.write(" 0) sky						#Component type\n")
		WP.write(" 1) 0.0000		0			#	Sky background at center of fitting region [ADUs]\n")
		WP.write(" 2) 0.0000		0			#	dsky/dx (sky gradient in x) [ADUs/pix]\n")
		WP.write(" 3) 0.0000		0			#	dsky/dy (sky gradient in y) [ADUs/pix]\n")
		WP.write(" Z) 0						#Leave in [1] or subtract [0] this comp from data?\n")
	
		WP.close()
	
		'''
		# run galfit on paramter file
		os.system('galfit ' + galfit_parameter_filename)
		
		# open logfile and use readin.py to append galfit result summaries to data_lower and data_upper
		logfile = open('fit.log') 
		read_param(logfile, int(distance))
	
		# remove temp files
		os.system("rm " + "fit.log")
		'''

if __name__ == "__main__":

	filePatternToMatch = "*_g.fits"
	imageListFilename = "images.txt"

	'''
	usageStr = ("Usage: command line should be 'python <full path>/rungalfit.py " +
				"<full path to directory containing images>")
	if len(sys.argv) != 2:
		print usageStr
	elif not os.path.isdir(sys.argv[1]):
		print "directory given on command line does not exist or is not accessible\n" + usageStr
	else:
		imagesDirectory = sys.argv[1]
	
		print "input files exist in the directory {}".format(imagesDirectory)
		if imagesDirectory[-1] != "/":
			imagesDirectory = imagesDirectory + "/"
			
		os.system("ls " + sys.argv[1] + filePatternToMatch + " > " + imageListFilename)
	'''
		
	os.system("ls " + filePatternToMatch + " > " + imageListFilename)
	#init_galfit_parameter_files()
	#run_galfit(imageListFilename)
	os.system("rm " + imageListFilename)




