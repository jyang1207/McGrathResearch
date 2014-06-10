#! /usr/bin python

#this script will run imexam and minmax to get the x,y and magnitude values for use in GALFIT

# '[1:' + width + ',1:' + str(int(height)/2) + ']'

import os
import sys
import re
from pyraf import iraf

os.system('touch'+' coords.tmp')			
def run_imexam(image):

# We're not splitting the image in half, also this only works if single parameter
# of the form -> [filename + '[1:' + width + ',1:' + str(int(height)/2) + ']']
# 	print input
# 	# [filename + '[1:' + width + ',1:' + str(int(height)/2) + ']']
# 	
# 	inputstring = str(input)
# 	# "[filename + '[1:' + width + ',1:' + str(int(height)/2) + ']']"
# 	
# 	input1 = str.split(inputstring, ",")[1]
# 	# ["[filename + '[1:' + width + '", "1:' + str(int(height)/2) + ']']"]
# 	# "1:' + str(int(height)/2) + ']']"
# 	
# 	input2 = str.split(input1, ":")[0]
# 	"1"
# 	
# 	#Note 'offset' is used to correct for only referencing half of the image
# 	if input2 == "1":													
# 		offset = int(0)
# 	else: 
# 		offset = int(input2)

	#iraf.images() -> I.T. 6/9 - I don't know what this does and commenting it out seems to have no effect
	
	
	#for image in input: -> wrong way of passing parameters
	
	# globals are bad, fix this later TODO
	global line, col, mag, radius, b_a, PA	
	
	#writes the output of the minmax function of iraf to a file for later use. 
	# TODO why append? might be a bad thing
	write_coords = open('coords.tmp', 'w')
	
	#call iraf's minmax function
	maxCoords = iraf.minmax(image, Stdout=1, update=0
		)[0].strip().split(" ")[3].replace("[", "").replace("]", "").split(",")

	# "broadband_200_cam0_g.fits[1:1893,1:946] [1765,691] -0.2845799624919891 [1283,692] 19.96048927307129"
	#xnew=str.split(X," ")
	# ["broadband_200_cam0_g.fits[1:1893,1:946]", "[1765,691]", "-0.2845799624919891", "[1283,692]", "19.96048927307129"]
	#x1 = (xnew[7])[1:]	
	# "1283, 692]"	
	#x2 = x1[:-1]	
	# "1283, 692"	
	#x3 = str.split(x2,",")		
	# ["1283", "692"]	
	#x4 = x3[0]		
	# "1283"
	#x5 = int(x3[1]) # + offset
	# 692
	
	write_coords.write(maxCoords[0]) 				
	write_coords.write(" " + maxCoords[1])
	write_coords.close()						
	imexam_out = str(iraf.imexam(image, use_display=0, imagecur="coords.tmp", Stdout=1))
	os.system('rm' + ' coords.tmp')		
	# "display frame (1:) (1): ['#   COL    LINE    COORDINATES', '#     " +
	# "R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT', " +
	# "'1283.40  692.16 1283.40 692.16', '  63.22  14.53  15461.  0.1202   19.53 0.52    " +
	# "1 1.75    21.58    18.85  21.39']"	
	
	# TODO: this should probably just return the above string, allowing the caller to parse for info
	
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
			
	
if __name__ == "__main__":
	run_imexam(sys.argv[1:])

#to initiate program, type 'iraf', cd into correct directory then type statement below
#in terminal import python pyraftest.py "broadband_xxx_cam0_g.fits[1:946,1:473]"


