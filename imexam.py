#! /usr/bin python

#this script will run imexam and minmax to get the x,y and magnitude values for use in GALFIT

import os
import sys
import re
from pyraf import iraf

os.system('touch'+' coords.tmp')			
def run_imexam(input):
	print input
	inputstring = str(input)
	input1 = str.split(inputstring, ",")[1]
	input2 = str.split(input1, ":")[0]
	
	if input2 == str(1):													#Note 'offset' is used to correct for only referencing half of the image
		offset = int(0)
	else: 
		offset = int(input2)

	iraf.images()
	for image in input:
		global line, col, mag, radius, b_a, PA								#writes the output of the minmax function of iraf to a file for later use. 
		write_coords = open('coords.tmp', 'a')
		X=str(iraf.minmax(image, Stdout=1, update=0))
		xnew=str.split(X," ")
		x1 = (xnew[7])[1:]							
		x2 = x1[:-1]								
		x3 = str.split(x2,",")						
		x4 = x3[0]									
		x5 = int(x3[1]) + offset													
		write_coords.write(str(x4)) 				
		write_coords.write(" " + str(x5))
		write_coords.close()						
		imexam_out = str(iraf.imexam(image, use_display=0, imagecur="coords.tmp", Stdout=1))				
		imexam_array = str.split(imexam_out, "'")
		
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
		#print col
		#print line
		#print mag
		#print radius
		#print b_a
		#print PA
		
		os.system('rm' + ' coords.tmp')				
	
if __name__ == "__main__":
	run_imexam(sys.argv[1:])

#to initiate program, type 'iraf', cd into correct directory then type statement below
#in terminal import python pyraftest.py "broadband_xxx_cam0_g.fits[1:946,1:473]"


