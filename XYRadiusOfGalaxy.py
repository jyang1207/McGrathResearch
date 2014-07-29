#!/usr/bin/env python

import os
import sys
# getting the parameters for radius and location of a galaxy from generic.cat

	
def getXYRadius(catalogueFileName):

	content = open(catalogueFileName, 'r')
	lines= content.readlines()
	content.close()
	#print lines
	indexDict = {}
	
	errorStr = ""
	
	# join method is to handle the relative path names
	outPutFileName = ".".join(catalogueFileName.split(".")[:-1])+".reg"
	
	outFile = open(outPutFileName, 'w')
	
	outFile.write(
'''#Region file format: DS9 version 4.1 
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 
physical
''')
	
	for line in lines:
		galaxyOutputList = line.split()
		
		if galaxyOutputList[0] == "#":
			indexDict[galaxyOutputList[2].upper()] = int(galaxyOutputList[1]) - 1
		else:
			if "X_IMAGE" in indexDict:
				galaxyX = galaxyOutputList[indexDict["X_IMAGE"]]
			else:
				errorStr = (errorStr + 
					"X_IMAGE is a required field of the parameter file\n")
			
			if "Y_IMAGE" in indexDict:
				galaxyY = galaxyOutputList[indexDict["Y_IMAGE"]]
			else:
				errorStr = (errorStr + 
					"Y_IMAGE is a required field of the parameter file\n")
			
			if "FLUX_RADIUS" in indexDict:
				galaxyRadius = galaxyOutputList[indexDict["FLUX_RADIUS"]]
			else:
				errorStr = (errorStr + 
					"FLUX_RADIUS is a required field of the parameter file\n")				
						
			if errorStr:
				print(errorStr)
				return False	
				
			circleString = "circle("+galaxyX + "," + galaxyY+ ","+ galaxyRadius+")\n"
			
			outFile.write(circleString)	
	outFile.close() 
	return True
	
if __name__=="__main__":
	if not os.path.isfile(sys.argv[1]):
		print("argument is not a file")
	else:
		filenames = open(sys.argv[1],'r')
		catalogs = filenames.readlines()
		filenames.close()
		
		for catalogFilename in catalogs:
			if not getXYRadius(catalogFilename.strip()):
				print "failed because of parameter file missing parameters."
		
	
	
	
	
	
	
	
	
	




