#! /usr/bin python

#this script will find the image dimensions for use in the imexam script

import os
import sys
import re
from pyraf import iraf

def run_imhead(image):

	#iraf.images() -> I.T. 6/9 - I don't know what this does and commenting it out seems to have no effect
	
	# globals are bad, fix this later TODO
	global frame_width, frame_height, image_number, cam_number
	
	M = iraf.imhead(image, Stdout=1)[0]
	# "broadband_200_cam0_g.fits[1893,1893][real]: " 
	
	M2 = str.split(M,"[")
	# ["broadband_200_cam0_g.fits", "1893,1893]", "real]: "]
	
	M3 = str.split(M2[1],"]")
	# ["1893,1893", ""]
	
	M4 = str.split(M3[0],",")
	# ["1893", "1893"]
	
	frame_width = M4[0]
	# "1893"
	
	frame_height = M4[1]
	# "1893"
	
	M5 = str.split(M,"_")
	# ["broadband", "200", "cam0", "g.fits[1893,1893][real]: "]
	
	image_number = M5[1]
	# "200"
	
	M6 = M5[2]
	# "cam0"
	
	#these statements accomodate for camera numbers greater than 9
	conditions = ["1","2","3","4","5","6","7","8","9"]			
	
	#determines if camera number is greater than 10, does not allow greater than 99
	# ex: "cam0" -> M6[-2] = "m"
	# ex: "cam12" -> M6[-2] = "1"
	if M6[-2] in conditions:										
		cam_number = M6[-2:]
	else:
		cam_number = int(M6[-1])						
	#print cam_number
		
if __name__ == "__main__":
	run_imhead(sys.argv[1:])

#to initiate program, type 'iraf', cd into correct directory then type statement below
#in terminal import python dimensions.py "broadband_xxx_cam0_g.fits"
