# classifyGalaxies.py
# Jianing Yang
# McGrath research

import csv
import numpy as np
import sys
from optparse import OptionParser

class Classifier:
	def __init__(self, filenames = None):
		# create and initialize fields for the class

		#list of all headers
		self.headers = []
		
		#list of all types
		self.types = []
		
		#list of lists of all data for three types of components. Each row is a list of strings
		self.central_data = []
		self.disk_data = []
		self.bulge_data = []
		
		#list of lists of all MRP data for three types of components. Each row is a list of strings
		self.central_dataMRP = []
		self.disk_dataMRP = []
		self.bulge_dataMRP = []
		
		#dictionary mapping header string to index of column in raw data
		self.header2col = {}
		
		if filenames != None:
			f = file(filenames, 'rU')
			fread = csv.reader(f)
			filenames = []
			for row in fread:
				filenames.append(row)
			print filenames
			self.readFiles(filenames)
		
	def readFiles(self, filenames):
		#read in all the csv files and put the data into a matrix of strings
		for filename in filenames:
			fp = file(filename[0], 'rU')
			cread = csv.reader(fp)
			
			self.headers = cread.next()
			self.types = cread.next()
			
			i = 0
			j = 0
			#read in data for galaxies in the list
			for row in cread:
				if row[2] in ["VELA02", "VELA03", "VELA04", "VELA05", "VELA12", "VELA14", "VELA15", "VELA26", "VELA27", "VELA28"]:
					if row[1] == "central":
						self.central_data.append(row)
						i += 1
					elif row[1] == "disk":
						self.disk_data.append(row)
						j += 1
					elif row[1] == "bulge":
						self.bulge_data.append(row)
				elif row[2] in ["VELA02MRP", "VELA03MRP", "VELA04MRP", "VELA05MRP", "VELA12MRP", "VELA14MRP", "VELA15MRP", "VELA26MRP", "VELA27MRP", "VELA28MRP"]:
					if row[1] == "central":
						self.central_dataMRP.append(row)
						i += 1
					elif row[1] == "disk":
						self.disk_dataMRP.append(row)
						j += 1
					elif row[1] == "bulge":
						self.bulge_dataMRP.append(row)
			print filename[0] + ' ' + str(i) + ', ' + str(j)
			
			for i in range(len(self.headers)):
				self.header2col[self.headers[i]] = i
				
	#classify the galaxies into one-component and two-component according to (rff1-rff2)/rff1
	def classifyComponent(self, isMRP=False):
		if isMRP:
			central_data = self.central_data
			disk_data = self.disk_data
			bulge_data = self.bulge_data
		else:
			central_data = self.central_dataMRP
			disk_data = self.disk_dataMRP
			bulge_data = self.bulge_dataMRP
		
		if len(central_data) != len(disk_data):
			print "Numbers of galaxies for single component and bulge component don't match."
			print str(len(central_data)) + ', ' + str(len(disk_data))
			return
		
		#column indicating whether one-component fit or two-component fit is better
		self.headers.append("comp_num")
		self.types.append("enum")
		self.header2col["comp_num"] = len(self.headers) - 1
		
		for i in range(len(central_data)):
			rff1 = float(central_data[i][self.header2col["rff"]])
			rff2 = float(disk_data[i][self.header2col["rff"]])
			rff_ratio = (rff1-rff2)/rff1
			
			#the cut-off rff ratio is 0.2
			if rff_ratio < 0.2:
				central_data[i].append("1")
				disk_data[i].append("1")
				bulge_data[i].append("1")
			else:
				central_data[i].append("2")
				disk_data[i].append("2")
				bulge_data[i].append("2")
			
	#classify the galaxies into spheroids and disk galaxies
	def classifyShape(self, isMRP=False):
		if isMRP:
			central_data = self.central_data
			disk_data = self.disk_data
			bulge_data = self.bulge_data
		else:
			central_data = self.central_dataMRP
			disk_data = self.disk_dataMRP
			bulge_data = self.bulge_dataMRP
		
		#column indicating whether the galaxy is a spheroid or a disk
		self.headers.append("shape")
		self.types.append("enum")
		self.header2col["shape"] = len(self.headers) - 1
		
		#for one-component, the cut-off sersic index is 2.5
		for i in range(len(central_data)):
			if central_data[i][self.header2col["comp_num"]] == "1":
				if float(central_data[i][self.header2col["sersicIndex(n)"]]) < 2.5:
					central_data[i].append("disk")
					disk_data[i].append("disk")
					bulge_data[i].append("disk")
				else:
					central_data[i].append("spheroid")
					disk_data[i].append("spheroid")
					bulge_data[i].append("spheroid")
					
			else:
				bMag = float(bulge_data[i][self.header2col["mag"]])
				dMag = float(disk_data[i][self.header2col["mag"]])
				
				bFlux = 10**(-0.4*bMag)
				dFlux = 10**(-0.4*dMag)
				tFlux = bFlux + dFlux
				
				if bFlux/tFlux < 0.5:
					central_data[i].append("disk")
					disk_data[i].append("disk")
					bulge_data[i].append("disk")
				else:
					central_data[i].append("spheroid")
					disk_data[i].append("spheroid")
					bulge_data[i].append("spheroid")
				
	#write the new data with component and shape information into 3 csv files
	def writeFiles(self, outputFilename, isMRP=False):
		if isMRP:
			central_data = self.central_data
			disk_data = self.disk_data
			bulge_data = self.bulge_data
		else:
			central_data = self.central_dataMRP
			disk_data = self.disk_dataMRP
			bulge_data = self.bulge_dataMRP
		
		outputFilename = outputFilename.split('.')[0]
		singleFilename = outputFilename + "_single.csv"
		diskFilename = outputFilename + "_disk.csv"
		bulgeFilename = outputFilename + "_bulge.csv"
		
		with open(singleFilename, 'wb') as file:
			fWriter = csv.writer(file, delimiter = ',')
			fWriter.writerow(self.headers)
			fWriter.writerow(self.types)
			fWriter.writerows(central_data)
			
		with open(diskFilename, 'wb') as file:
			fWriter = csv.writer(file, delimiter = ',')
			fWriter.writerow(self.headers)
			fWriter.writerow(self.types)
			fWriter.writerows(disk_data)
			
		with open(bulgeFilename, 'wb') as file:
			fWriter = csv.writer(file, delimiter = ',')
			fWriter.writerow(self.headers)
			fWriter.writerow(self.types)
			fWriter.writerows(bulge_data)
			
	#write the new data combined into 1 csv file according to component number and shape
	def writeCombined(self, outputFilename, isMRP=False):
		if isMRP:
			central_data = self.central_data
			disk_data = self.disk_data
			bulge_data = self.bulge_data
		else:
			central_data = self.central_dataMRP
			disk_data = self.disk_dataMRP
			bulge_data = self.bulge_dataMRP
		
		data = central_data
		
		for i in range(len(data)):
			if data[i][self.header2col["comp_num"]] == "2":
				if data[i][self.header2col["shape"]] == "disk":
					data[i] = disk_data[i]
				else:
					data[i] = bulge_data[i]
		
		with open(outputFilename, 'wb') as file:
			fWriter = csv.writer(file, delimiter = ',')
			fWriter.writerow(self.headers)
			fWriter.writerow(self.types)
			fWriter.writerows(data)
			
if __name__ == '__main__':
	if len(sys.argv) < 3:
		print "Usage: python %s <csv_filename>" % sys.argv[0]
		print "		  where <csv_filename> specifies a csv file"
		exit()		
	
	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# indicate that the galaxies are MRP
	parser.add_option("-m", "--isMRP",
					  help="to plot MRP counterparts adjacent",
					  action="store_true")

	# indicate whether to write a combined file or seperate files
	parser.add_option("-c", "--combined",
					  help="to write a combined file",
					  action="store_true")
					  
	# parse the command line using above parameter rules
	# options - list with everthing defined above, 
	# args - anything left over after parsing options
	[options, args] = parser.parse_args()
		
	# classify the a0 or MRP galaxies by component number and shape
	classifier = Classifier(sys.argv[1])
	classifier.classifyComponent(isMRP=options.isMRP)
	classifier.classifyShape(isMRP=options.isMRP)
	
	if options.combined:
		classifier.writeCombined(sys.argv[2])
		print "combined files written"
	else:
		classifier.writeFiles(sys.argv[2])
		print "seperate files written"
	
	
		
		
		
		
		
		
		
		
		