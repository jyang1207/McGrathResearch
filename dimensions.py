#! /usr/bin python

#this script will find the image dimensions for use in the imexam script

import os
import sys
import re
from pyraf import iraf

def run_imhead(input):
	iraf.images()
	for image in input:
		global frame_width, frame_height, image_number, cam_number
		M = str(iraf.imhead(image, Stdout=1))
		M2 = str.split(M,"[")
		M3 = M2[2]		
		M3 = str.split(M2[2],"]")				
		M4 = str.split(M3[0],",")
		frame_width = M4[0]
		frame_height = M4[1]
		M5 = str.split(M,"_")
		image_number = M5[1]
		M6 = M5[2]
		conditions = ["1","2","3","4","5","6","7","8","9"]				#these statements accomodate for camera numbers greater than 9
		if M6[-2] in conditions:										#determines if camera number is greater than 10
			cam_number = M6[-2:]
		else:
			cam_number = int(M6[-1])						
		#print cam_number
		
if __name__ == "__main__":
	run_imhead(sys.argv[1:])

#to initiate program, type 'iraf', cd into correct directory then type statement below
#in terminal import python dimensions.py "broadband_xxx_cam0_g.fits"

