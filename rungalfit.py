
#! /usr/bin/python


import os
import sys
import re
from pyraf import iraf
import imexam					#works
import dimensions
import readin

os.system("touch " + "data_lower.txt")
format = open("data_lower.txt", "w")
format.write("Input Image" + "	" + "Time Step" + "	" + "X Value" + "	" + "X Error" + "	" + "Y Value" + "	" + "Y Error" + "	" + "Magnitude" + "	" + "Magnitude Error" + "	" + "Radius" + "	" + "Radius Error" + "	" + "Sersic" + "	" + "Sersic Error" + "	" + "BA Value" + "	" + "BA Error" + "	" + "Position Angle" + "	" + "Position Angle Error" + "	" + "Chi Squared" + "	" + "ndof" + "	" + "Chi Squared/nu" + "\n")
format.close()

os.system("touch " + "data_upper.txt")
format2 = open("data_upper.txt", "w")
format2.write("Input Image" + "	" + "Time Step" + "	" + "X Value" + "	" + "X Error" + "	" + "Y Value" + "	" + "Y Error" + "	" + "Magnitude" + "	" + "Magnitude Error" + "	" + "Radius" + "	" + "Radius Error" + "	" + "Sersic" + "	" + "Sersic Error" + "	" + "BA Value" + "	" + "BA Error" + "	" + "Position Angle" + "	" + "Position Angle Error" + "	" + "Chi Squared" + "	" + "ndof" + "	" + "Chi Squared/nu" + "\n")
format2.close()

f = open("images.txt")		#this will be the file that will contain the images
line = f.readline()

while line:
	filename = line.rstrip()
	line = f.readline()
	
	dimensions.run_imhead([filename])
	
	camera = dimensions.cam_number					
	image = dimensions.image_number
	height = dimensions.frame_height
	width = dimensions.frame_width
	
	imexam.run_imexam([filename + '[1:' + width + ',1:' + str(int(height)/2) + ']'])
	
	angle = imexam.PA																		#redefine variables in this script for convenience
	BA = imexam.b_a
	rad = imexam.radius
	magnitude = imexam.mag
	X = imexam.col
	Y = imexam.line
	
	imexam.run_imexam([filename + '[1:' + width + ',' + str((int(height)/2)+1) + ':' + height + ']'])
	
	angle_top = imexam.PA
	BA_top = imexam.b_a
	rad_top = imexam.radius
	magnitude_top = imexam.mag
	X_top = imexam.col
	Y_top = imexam.line
	
	delta_x = float(X_top) - float(X)														#this block computes the distance between objects
	delta_y = float(Y_top) - float(Y)
	delta_x_square = int(delta_x)**2
	delta_y_square = int(delta_y)**2
	#global distance
	distance = int((delta_x_square + delta_y_square)**(0.5))
	print distance
	print str(distance)
	os.system("touch " + "dist.tmp")
	write_dist = open('dist.tmp', 'w')
	write_dist.write(str(distance))
	write_dist.close()
	
	if distance <= 20:																		#adjustable parameter for when to run galfit on second component						
		Z = "# "
	else:
		Z = " "

	os.system('touch' + ' bb' + image + '_c' + str(camera) + '.txt')								#below this writes the galfit parameter file
	def make_galfit(line):
		WP = open('bb' + image + '_c' + str(camera) + '.txt','w')
		WP.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
		WP.write("A)" + " broadband_" + image + "_cam" + str(camera) +"_g.fits" + " 			#Input data image block\n")
		WP.write("B)" + " bb" + image + "_c" + str(camera) + "_multi.fits" + "					#Output data image block\n")
		WP.write("C)" + " none" + "									#Sigma image name (made from data if blank or 'none')\n")
		WP.write("D)" + " ps4.fits" + "								#Input PSF image and (optional) diffusion kernel\n")
		WP.write("E)" + " 1" + "									#PSF fine sampling factor relative to data\n")
		WP.write("F)" + " none" + "									#Bad pixel mask (FITS file or ASCIIcoord list)\n")
		WP.write("#G)" + " bb" + image + "_c" + str(camera) + "_constraint.txt" + "				#File with parameter constraints (ASCII file)\n")
		WP.write("H)" + " 1	" + width + "  1	" + height + "						#Image region to fit (xmin xmax ymin ymax)\n")
		WP.write("I)" + " 200	" + "200" + "								#Size of the concolution box (x y)\n")
		WP.write("J)" + " 25.0" + "								#Magnitude photometric zeropoint\n")			#this will be user inputed 
		WP.write("K)" + " 0.038" + "  0.038" + "							#Plate scale (dx dy)  [arcsec per pixel]\n")
		WP.write("O)" + " regular" + "								#display type (regular, curses, both\n")
		WP.write("P)" + " 0" + "									#Options: 0=normal run; 1,2=make model/imgblock & quit\n")
		WP.write("S)" + " 0" + "									#Modify/create components interactively?\n")
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
	
		os.system('galfit ' + 'bb' + image + '_c' + str(camera) + '.txt')
		logfile = open('fit.log') 
		readin.read_param(logfile)
	
		os.system("rm " + "dist.tmp")
		os.system("rm " + "fit.log")
	
	if __name__ == "__main__":
		make_galfit(sys.argv[1:])
	

f.close()



