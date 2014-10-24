#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 8/13/2104
Colby College Astrophysics Research
'''

import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
import matplotlib.ticker  as ticker
import numpy as np
import os
# import itertools
from optparse import OptionParser
from pprint import pprint
from collections import OrderedDict

	
def vararg_callback(option, opt_str, value, parser):
	'''
	allows for galaxy names to be given as an arbitrary list spearated by spaces
	'''
	assert value is None
	value = []
	
	def floatable(strParam):
		try:
			float(strParam)
			return True
		except ValueError:
			return False
	
	for arg in parser.rargs:
		# stop on --foo like options
		if arg[:2] == "--" and len(arg) > 2:
			break
		# stop on -a, but not on -3 or -3.0
		if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
			break
		value.append(arg)
	
	del parser.rargs[:len(value)]
	setattr(parser.values, option.dest, value)
	

def getFormats():
	'''
	generates a list of formats for plot markers of various shape and color
	'''
	colors = ['b', 'g', 'r', 'c', 'y', 'm']
	markers = ['o', '^', 'd', 'h', '*', '+']
	formats = []
	# use itertools to create product of the list if this isnt enough
	for index, color in enumerate(colors):
		formats.append(color + markers[index])
	
	for index, color in enumerate(colors):
		formats.append(color + markers[(index+1)*-1])
	
	return formats

def plotAllCamera(	plotComponents, fieldDescriptions, xFieldName='age', yFieldName='ser', curSubPlot=plt, includeLegend = True):
	'''
	plot all cameras in plot components and optionally include legend
	plots the components given on the plot given, using fieldDescriptions to set axis scales
	'''
	formats = [fStr + "-" for fStr in getFormats()]
	
	# plot sersic index vs age, separate for each camera angles
	ageSet = np.unique(plotComponents[xFieldName])
	
	camSet = np.unique(plotComponents['cam'])
	fieldValByAge = dict([[age, []] for age in ageSet])
	for componentNumber, fieldVal in enumerate(plotComponents[yFieldName]):
		fieldValByAge[plotComponents[xFieldName][componentNumber]].append(
								[fieldVal,plotComponents['cam'][componentNumber]])

	for ind, cam in enumerate(camSet):
		xVals = []
		yVals = []
		for age in ageSet:
			found = False
			for [fieldVal,curCam] in fieldValByAge[age]:
				if curCam == cam and not found:
					yVals.append(fieldVal)
					found = True
			if found:
				xVals.append(age)
		curSubPlot.plot(xVals, yVals, formats[ind], label='Camera '+str(cam))
		curSubPlot.xlim([fieldDescriptions[xFieldName][2],fieldDescriptions[xFieldName][3]])
		curSubPlot.ylim([fieldDescriptions[yFieldName][2],fieldDescriptions[yFieldName][3]])
		#if yFieldName == "rad":
		#	 curSubPlot.set_yscale("log")
		if includeLegend:
			curSubPlot.legend(loc='upper left', prop={'size':10})
	
def plotAvgCamera(plotComponents, fieldDescriptions, xFieldName='age', yFieldName='ser', curSubPlot=plt):
	'''
	average cameras in plot components and plot errorbars
	plots the components given on the plot given, using fieldDescriptions to set axis scales
	'''
	# plot sersic index vs age, averaging all camera angles
	ageSet = np.unique(plotComponents[xFieldName])
	fieldValByAge = dict([[age, []] for age in ageSet])
	for componentNumber, fieldVal in enumerate(plotComponents[yFieldName]):
		fieldValByAge[plotComponents[xFieldName][componentNumber]].append(fieldVal)
	
	yVals = []
	yErr = []
	for age in ageSet:
		yVals.append(np.mean(fieldValByAge[age]))
		yErr.append(np.std(fieldValByAge[age]))
	
	xVals = ageSet
	
	# compute a linear fit to overplot
	#fit = pyl.polyfit(xVals, yVals, 1)
	#print(fit)
	#fitFn = pyl.poly1d(fit)
	#curSubPlot.plot(xVals, yVals, 'bo', xVals, fitFn(xVals), 'b-')
	curSubPlot.errorbar(xVals, yVals, yerr=yErr, ecolor='r', linewidth=0.1, fmt="ro")
	curSubPlot.xlim([fieldDescriptions[xFieldName][2],fieldDescriptions[xFieldName][3]])
	curSubPlot.ylim([fieldDescriptions[yFieldName][2],fieldDescriptions[yFieldName][3]])
	
	return

# TODO not implemented after revision
def calcBulgeToTotal(bulgeMag, diskMag):
	'''
	method comment
	
	bulgeToTotal = 1 / (np.power([10], [0.4 * (bulgeMag - diskMag)]) + 1)
	return bulgeToTotal
	'''
	return

# TODO not implemented after revision
def getBulgeToTotalRatios(allFieldNames, data, galaxyName=""):
	'''
	method comment
	
	exclusionList = getBulgeErrors(data, galaxyName)
	bulgeGalaxies = getGalaxies(allFieldNames, data, "bulge", galaxyName, exclusionList)
	diskGalaxies = getGalaxies(allFieldNames, data, "disk", galaxyName, exclusionList)
	for componentNumber, bulgeMag in enumerate(bulgeGalaxies['mag']):
		diskMag = diskGalaxies['mag'][componentNumber]
		bulgeToTotal = calcBulgeToTotal(bulgeMag, diskMag)
		diskGalaxies['mag'][componentNumber] = bulgeToTotal
	
	return diskGalaxies
	'''
	return
	
def getBulgeErrors(data, galaxyName):
	'''
	used by getGalaxies to only include bulge with disks and visa versa, prints excluded
	'''
	magByAge = {}
	for componentNumber, curName in enumerate(data["id"]):
		curAge =	data['ts'][componentNumber]
		curType =	data['typ'][componentNumber]
		curMag =	data['mag'][componentNumber]
		curCam =	data['cam'][componentNumber]
		if ((not galaxyName) or (curName == galaxyName)) and (curType in ["bulge","disk"]):
			curKey = str(curName)+"_"+str(curAge)+"cam"+str(curCam)
			if not (curKey in magByAge):
				magByAge[curKey] = [curMag, curType]
			else:
				magByAge[curKey] = [[curMag, curType], magByAge[curKey]]
				
	excludeList = []
	for curKey in magByAge:
		curMags = magByAge[curKey]
		try:
			if curMags[0][1] == curMags[1][1]:
				print(" ".join(["excluding",str(curKey),str(curMags)]))
				excludeList.append([str(curKey).split("cam")[0].split("_")[1], 
									str(curKey).split("cam")[1]])
		except:
			print(" ".join(["excluding",str(curKey),str(curMags)]))
			excludeList.append([str(curKey).split("cam")[0].split("_")[1], 
								str(curKey).split("cam")[1]])
	return excludeList

def getGalaxies(allFieldNames, data, galaxyType="central", galaxyName=""):
	'''
	filters data by type and name of galaxy, returning matching galaxies
	'''
	excluded=[]
	exclusionList = getBulgeErrors(data, galaxyName)
	
	# filter records down to just the galaxies of given type
	typedGalaxies = dict([[fieldName, []] for fieldName in allFieldNames])
	for componentNumber, curType in enumerate(data['typ']):
		curName =	data['id'][componentNumber]
		curTS =		data['ts'][componentNumber]
		curCam =	data['cam'][componentNumber]
		
		# if not excluded and of the correct type
		if (curType == galaxyType) and not([str(curTS), str(curCam)] in exclusionList):
		
			# if the the galaxy name is not restricted or if it matches
			if (not galaxyName) or (curName == galaxyName):
			
				# append all the data for the current component to the galaxies to be returned
				for fieldName in allFieldNames:
					typedGalaxies[fieldName].append(data[fieldName][componentNumber])
		
		# if not added because of exclusion list, add to list of excluded galaxies
		if ((curType == galaxyType) and (curName == galaxyName) and 
			([str(curTS), str(curCam)] in exclusionList)):
			excluded.append([str(curTS), str(curCam)])

	return typedGalaxies
	
def getData(summaryFilename, fieldDescriptions):
	'''
	given a summary file, return a formatted dictionary of field name to values
	
	parameter summaryFilenames - 
		summary filename containing the data in space delimited columns
	parameter fieldDescriptions - 
		a list of tuples containing field names and format strings
	returns - 
		dictionary with field names as keys and array of field values as value
	'''

	# read the columns into a 2D array with names and formats as above
	rawData = pyl.loadtxt(summaryFilename, 
						  dtype={'names':list(fieldDescriptions.keys()),
								 'formats':[fieldDescriptions[name][0] for name in fieldDescriptions]},
						  skiprows=2)
	
	# dictionary with field names as keys and array of field values as value
	data = {}
	for colIndex, name in enumerate(fieldDescriptions):
		
		# loadtxt does an annoying byte array thing, decode undoes it so strings are strings
		if fieldDescriptions[name][0][0] in ['a','S']:
			data[name] = [component[colIndex].decode('utf-8') for component in rawData]
		else:
			data[name] = [component[colIndex] for component in rawData]

	return data
	
	
if __name__ == "__main__":
	
	# master list of all available plot types
	plotTypes = ["default", "allGalaxies", "allFields", "bulgeToTotal", "special"]
	
	#define the command line interface with simUtility.py
	usage = ("\n%prog summaryFile [-h help] [options (with '-'|'--' prefix)]")
			
	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-p","--plotType",
					  help="the type of plot, available options are: " + str(plotTypes),
					  default=plotTypes[0])
	
	# indicate that a there are MRP counterparts to all galaxy names
	parser.add_option("-m","--includeMRP", 
					  help="if all galaxy names have corresponding MRP, plot them too",
					  action="store_true")

	# indicate that you want all cameras plotted separately
	parser.add_option("-c","--allCameras", 
                      help="to show all cameras",
                      action="store_true")

	# indicate that you want all cameras plotted separately
	parser.add_option("-r","--candelized", 
                      help="to indicate candelized results are being plotted",
                      action="store_true")
                      
	# indicate that you have errors in the summary file
	parser.add_option("-e","--errors", 
                      help="to indicate summary file has GALFIT errors for each field",
                      action="store_true")
	
	# pass the list of galaxies
	parser.add_option("-n","--galaxyNames",
					  help="the space separated list of galaxy names to be plotted (must exist in summary file)",
					  dest="galaxyNames",
					  action="callback", callback=vararg_callback, default=[""])
	
	# pass the list of y field names
	parser.add_option("-y","--yFields",
					  help="the list of y field names to be plotted (see field descriptions variable)",
					  dest="yFields",
					  action="callback", callback=vararg_callback, default=[])
	
	# pass the field name to be plotted on the x axis
	parser.add_option("-x","--xFieldName",
					  help="the field name of the x values (see field descriptions variable)",
					  default="red")
	
	# pass the type of component to plot
	parser.add_option("-t","--componentType",
					  help="the type of component to be plotted (central, bulge, disk)",
					  default="central")
	
	# parse the command line using above parameter rules
	# options - list with everthing defined above, 
	# args - anything left over after parsing options
	[options, args] = parser.parse_args()
	
	if not len(args) or not os.path.isfile(args[0]):
		parser.error(	"summary file must be an existing file with space " +
						"delimited summary data")
		
	error = 5
	if options.candelized:
		poslow = 0
		poshigh = 100
	else:
		poslow = 285
		poshigh = 315
	# dictionary of lists [format, text, lower, upper], one for each field in summary file
	# careful changing this list, ordered to match the order of columns in summary file
	fieldDescriptions = OrderedDict()
	fieldDescriptions['typ'] =	 	['a10','Galaxy Type']
	fieldDescriptions['id'] =	 	['a10','Galaxy ID']
	fieldDescriptions['ts'] =	 	['a10','Time Step (a)']
	fieldDescriptions['age'] =	 	['f4','Age (Gyr)',0,8]
	fieldDescriptions['red'] =	 	['f4','Redshift (z)',5,0]
	fieldDescriptions['cam'] =	 	['i4','Camera Number']
	fieldDescriptions['fil'] =	 	['a10','Filter']
	fieldDescriptions['px'] =	 	['f4','X Position (pixels)',poslow,poshigh]
	if options.errors:
		fieldDescriptions['epx'] =	 	['f4','Error in x position (pixels)', -error, error]
	fieldDescriptions['py'] =	 	['f4','Y Position (pixels)',poslow,poshigh]
	if options.errors:
		fieldDescriptions['epy'] =	 	['f4','Error in y position (pixels)', -error, error]
	fieldDescriptions['mag'] =	 	['f4','Magnitude',30,15]
	if options.errors:
		fieldDescriptions['emag'] =	 	['f4','Error in magnitude', -error, error]
	fieldDescriptions['rpix'] =	 	['f4',r"$R_{eff}$ (pixels)",0.5,50]
	if options.errors:
		fieldDescriptions['erpix'] =	['f4','Error in radius (pixels)', -error, error]
	fieldDescriptions['rad'] =	 	['f4',r"$R_{eff}$ (kpc)",0.5,15]
	if options.errors:
		fieldDescriptions['erad'] =	 	['f4','Error in radius (kpc)', -error, error]
	fieldDescriptions['ser'] =	 	['f4','Sersic Index',0.05,8.5]
	if options.errors:
		fieldDescriptions['eser'] =	 	['f4','Error in sersic index', -error, error]
	fieldDescriptions['ba'] =	 	['f4','Axis Ratio',0.05,1.2]
	if options.errors:
		fieldDescriptions['eba'] =	 	['f4','Error in axis ratio', -error, error]
	fieldDescriptions['pa'] =	 	['f4','Position Angle (deg)',-180,180]
	if options.errors:
		fieldDescriptions['epa'] =	 	['f4','Error in position angle (deg)', -error, error]
	fieldDescriptions['sky'] =	 	['f4','sky value']
	
	###
	print("\nCommand line options...")
	###
	pprint(vars(options))
	print("summary file used: " + args[0])
	
	###
	print("\nValidating field names...")
	###
	# verify given y fields
	yFields = []
	for yFieldName in options.yFields:
		if (yFieldName in fieldDescriptions) and (len(fieldDescriptions[yFieldName]) > 3):
			yFields.append(yFieldName)
	if not yFields:
		print("\tno valid y fields specified, using sersic index as the default")
		yFields = ['ser']
		
	# verify given x field
	if options.xFieldName in fieldDescriptions:
		xFieldName = options.xFieldName
	else:
		print("\tno valid x field specified, using redshift as the default")
		xFieldName = 'red'
		
	###
	print("\nValidating plot type...")
	###
	if options.plotType in plotTypes:
		plotType = options.plotType
	else:
		print("\tno valid plot type specified, using '" + plotTypes[0] + "' as the default")
		plotType = plotTypes[0]
	
	###
	print("\nReading data from summary file...")
	###
	# dictionary with field names as keys (above) and array of field values as value
	data = getData(args[0], fieldDescriptions)
	print("\tData read successfully, using y fields:")
	for yFieldName in yFields:
		print("\t\t" + fieldDescriptions[yFieldName][1] + 
				": num elements = " + str(len(data[yFieldName])))
	print("\tagainst x field:")
	print("\t\t" + fieldDescriptions[xFieldName][1] + 
			": num elements = " + str(len(data[xFieldName])))
				
	###
	print("\nPlotting...")
	###
	# plot each galaxy name and y field in its own figure
	if plotType == "default":
		for galaxyName in options.galaxyNames:
		
			curData = getGalaxies(	fieldDescriptions.keys(), data, 
									options.componentType, galaxyName)
			if not curData[xFieldName]:
				print("no data for galaxy " + galaxyName)
				continue	
			numCols = 1
			if options.includeMRP:
				mrpData = getGalaxies(	fieldDescriptions.keys(), data, 
										options.componentType, galaxyName+"MRP")
				if not mrpData[xFieldName]:
					print("no data for galaxy " + galaxyName + "MRP")
				else:
					numCols = 2
					
			for yFieldName in yFields:
				xlabel = fieldDescriptions[xFieldName][1]
				ylabel = fieldDescriptions[yFieldName][1]
				plt.figure().suptitle(galaxyName + " " + ylabel + " vs " + xlabel)
						
				if options.allCameras:
					numRows = 2
					plt.subplot(numRows,numCols,1)
					plt.title("all cameras")
					plotAllCamera(curData, fieldDescriptions, xFieldName, yFieldName)
					plt.xlabel(xlabel)
					plt.ylabel(ylabel)
					if numCols == 2:
						plt.subplot(numRows,numCols,2)
						plt.title("all cameras")
						plotAllCamera(mrpData, fieldDescriptions, xFieldName, yFieldName)
						plt.xlabel(xlabel)
						plt.ylabel(ylabel + " (with MRP)")
						plt.subplot(numRows,numCols,4)
						plt.title("avg cameras")
						plotAvgCamera(mrpData, fieldDescriptions, xFieldName, yFieldName)
						plt.xlabel(xlabel)
						plt.ylabel(ylabel + " (with MRP)")
						plt.subplot(numRows,numCols,3)
					else:
						plt.subplot(numRows,numCols,2)
					plt.title("avg cameras")
					plotAvgCamera(curData, fieldDescriptions, xFieldName, yFieldName)
					plt.xlabel(xlabel)
					plt.ylabel(ylabel)
				else:
					numRows = 1
					plt.subplot(numRows,numCols,1)
					plotAvgCamera(curData, fieldDescriptions, xFieldName, yFieldName)
					plt.xlabel(xlabel)
					plt.ylabel(ylabel)
					if numCols == 2:
						plt.subplot(numRows,numCols,2)
						plotAvgCamera(mrpData, fieldDescriptions, xFieldName, yFieldName)
						plt.xlabel(xlabel)
						plt.ylabel(ylabel + " (with MRP)")
						
	# single component fit results, one plot per field, all galaxies
	elif plotType == "allGalaxies":
	
		if options.includeMRP:
			numRows = 2
		else:
			numRows = 1
			
		for yFieldName in yFields:		
			xlabel = fieldDescriptions[xFieldName][1]
			ylabel = fieldDescriptions[yFieldName][1]
			plt.figure().suptitle(ylabel + " vs " + xlabel)
			
			numCols = len(options.galaxyNames)
			for galaxyIndex, galaxyName in enumerate(options.galaxyNames, start=1):	
				curData = getGalaxies(	fieldDescriptions.keys(), data, 
										options.componentType, galaxyName)
				if not curData[xFieldName]:
					print("no data for galaxy " + galaxyName + " y field " + yFieldName)
					continue	
				if options.includeMRP:
					mrpData = getGalaxies(	fieldDescriptions.keys(), data, 
											options.componentType, galaxyName+"MRP")
					if not mrpData[xFieldName]:
						print("no data for galaxy " + galaxyName + "MRP")
				
				plt.subplot(numRows,numCols,galaxyIndex)
				plt.title(galaxyName)
				plotAvgCamera(curData, fieldDescriptions, xFieldName, yFieldName)
				plt.xlabel(xlabel)
				plt.ylabel(ylabel)
				
				if options.includeMRP and mrpData[xFieldName]:
					plt.subplot(numRows,numCols,galaxyIndex+numCols)
					plotAvgCamera(mrpData, fieldDescriptions, xFieldName, yFieldName)
					plt.xlabel(xlabel)
					plt.ylabel(ylabel + " (with MRP)")
				
	elif plotType == "allFields":
		if options.includeMRP:
			numRows = 2
		else:
			numRows = 1
			
		for galaxyName in options.galaxyNames:
		
			curData = getGalaxies(	fieldDescriptions.keys(), data, 
									options.componentType, galaxyName)
			if not curData[xFieldName]:
				print("no data for galaxy " + galaxyName)
				continue	
			if options.includeMRP:
				mrpData = getGalaxies(	fieldDescriptions.keys(), data, 
										options.componentType, galaxyName+"MRP")
				if not mrpData[xFieldName]:
					print("no data for galaxy " + galaxyName + "MRP")
					
			plt.figure().suptitle(galaxyName)
			numCols = len(yFields)
			
			for yFieldIndex, yFieldName in enumerate(yFields, start=1):		
				xlabel = fieldDescriptions[xFieldName][1]
				ylabel = fieldDescriptions[yFieldName][1]	
				
				plt.subplot(numRows,numCols,yFieldIndex)
				plotAvgCamera(curData, fieldDescriptions, xFieldName, yFieldName)
				plt.xlabel(xlabel)
				plt.ylabel(ylabel)
				
				if options.includeMRP and mrpData[xFieldName]:
					plt.subplot(numRows,numCols,yFieldIndex+numCols)
					plotAvgCamera(mrpData, fieldDescriptions, xFieldName, yFieldName)
					plt.xlabel(xlabel)
					plt.ylabel(ylabel + " (with MRP)")
					
	elif plotType == "bulgeToTotal":
		print("plot type '" + plotType + "' not yet implemented")
		exit()
		
	# this is for manual use, limited/inconsistent use of other command line flags
	elif plotType == "special":
		print("plot type '" + plotType + "' not yet implemented")
		exit()
			
	else:
		print("plot type '" + plotType + "' not yet implemented")
		exit()
		
	# display all figures resulting from above calls
	plt.show()


	
