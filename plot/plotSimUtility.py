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
import pprint

	
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


def plotAllCamera(plotComponents, xFieldName='age', fieldName='ser', curSubPlot=plt):
	'''
	
	'''
	formats = [fStr + "-" for fStr in getFormats()]
	
	# plot sersic index vs age, separate for each camera angles
	ageSet = np.unique(plotComponents[xFieldName])
	
	camSet = np.unique(plotComponents['cam'])
	fieldValByAge = dict([[age, []] for age in ageSet])
	for componentNumber, fieldVal in enumerate(plotComponents[fieldName]):
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
		#if fieldName == "rad":
		#	 curSubPlot.set_yscale("log")
		curSubPlot.legend(loc='upper left', prop={'size':10})
	

def plotAvgCamera(plotComponents, xFieldName='age', fieldName='ser', curSubPlot=plt):
	'''
	
	'''
	# plot sersic index vs age, averaging all camera angles
	ageSet = np.unique(plotComponents[xFieldName])
	fieldValByAge = dict([[age, []] for age in ageSet])
	for componentNumber, fieldVal in enumerate(plotComponents[fieldName]):
		fieldValByAge[plotComponents[xFieldName][componentNumber]].append(fieldVal)
	
	yVals = []
	yErr = []
	for age in ageSet:
		yVals.append(np.mean(fieldValByAge[age]))
		yErr.append(np.std(fieldValByAge[age]))
	
	xVals = ageSet
	'''
	# compute a linear fit to overplot
	fit = pyl.polyfit(xVals, yVals, 1)
	print(fit)
	fitFn = pyl.poly1d(fit)'''
	#curSubPlot.plot(xVals, yVals, 'bo', xVals, fitFn(xVals), 'b-')
	curSubPlot.errorbar(xVals, yVals, yerr=yErr, ecolor='r', linewidth=0.1, fmt="ro")
		

def getBulgeErrors(data, galaxyName):
	'''
	
	'''
	magByAge = {}
	for componentNumber, curName in enumerate(data["id"]):
		curAge =	data['ts'][componentNumber]
		curType =	data['typ'][componentNumber]
		curMag =	data['mag'][componentNumber]
		curCam =	data['cam'][componentNumber]
		if (curName == galaxyName) and (curType in ["bulge","disk"]):
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


def getGalaxies(allFieldNames, data, galaxyType="central", galaxyName="", exclusionList=[]):
	'''
	
	'''
	excluded=[]
	# filter records down to just the galaxies of given type
	typedGalaxies = dict([[fieldName, []] for fieldName in allFieldNames])
	for componentNumber, curType in enumerate(data['typ']):
		curName =	data['id'][componentNumber]
		curTS =		data['ts'][componentNumber]
		curCam =	data['cam'][componentNumber]
		if (curType == galaxyType) and not ([str(curTS), str(curCam)] in exclusionList):
			if (not galaxyName) or (curName == galaxyName):
				for fieldName in allFieldNames:
					typedGalaxies[fieldName].append(data[fieldName][componentNumber])
		if curType == galaxyType and curName == galaxyName and ([str(curTS), str(curCam)] in exclusionList):
			excluded.append([str(curTS), str(curCam)])

	return typedGalaxies
	
	
def getData(summaryFilename, fieldDescriptions):
	'''
	parameter summaryFilenames - 
		summary filename containing the data in space delimited columns
	parameter fieldDescriptions - 
		a list of tuples containing field names and format strings
	returns - 
		dictionary with field names as keys and array of field values as value
	'''

	# read the columns into a 2D array with names and formats as above
	rawData = pyl.loadtxt(summaryFilename, 
						  dtype={'names':[row[0] for row in fieldDescriptions],
								 'formats':[row[1] for row in fieldDescriptions]},
						  skiprows=2)
	
	# dictionary with field names as keys and array of field values as value
	data = {}
	for colIndex, desc in enumerate(fieldDescriptions):
		
		# loadtxt dies an annoying byte array thing, decode undoes it so strings are strings
		if desc[1][0] in ['a','S']:
			data[desc[0]] = [component[colIndex].decode('utf-8') for component in rawData]
		else:
			data[desc[0]] = [component[colIndex] for component in rawData]

	return data


def calcBulgeToTotal(bulgeMag, diskMag):
	'''
	
	'''
	bulgeToTotal = 1 / (np.power([10], [0.4 * (bulgeMag - diskMag)]) + 1)
	return bulgeToTotal


def getBulgeToTotalRatios(allFieldNames, data, galaxyName=""):
	'''
	
	'''
	exclusionList = getBulgeErrors(data, galaxyName)
	bulgeGalaxies = getGalaxies(allFieldNames, data, "bulge", galaxyName, exclusionList)
	diskGalaxies = getGalaxies(allFieldNames, data, "disk", galaxyName, exclusionList)
	for componentNumber, bulgeMag in enumerate(bulgeGalaxies['mag']):
		diskMag = diskGalaxies['mag'][componentNumber]
		bulgeToTotal = calcBulgeToTotal(bulgeMag, diskMag)
		diskGalaxies['mag'][componentNumber] = bulgeToTotal
	
	return diskGalaxies


def plotBulgeTotalRatio(fieldDescriptions, allFieldNames, data, 
						galaxyName="", figureIndex=1, xFieldIndex=3, MRP=False):
	'''
	
	'''
	xFieldName = fieldDescriptions[xFieldIndex][0]
	xLabel = fieldDescriptions[xFieldIndex][2]
	yLabel = 'Bulge/Total'
	axisLimits = [fieldDescriptions[xFieldIndex][3],
				  fieldDescriptions[xFieldIndex][4],
				  0.05, 1.2]
	fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, num=galaxyName+" bt")#figureIndex)
	axes = [[axes[0]],[axes[1]]]
	'''if not MRP:
		fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, num=figureIndex)
		axes = [[axes[0]],[axes[1]]]
	else:
		fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, num=figureIndex)
		axes[1][1].set_xlabel(xLabel)'''
		
	axes[0][0].set_ylabel(yLabel)#+"\n(without MRP)")#(by camera)")
	axes[0][0].text(0.1,0.8,"No RP",transform=axes[0][0].transAxes,fontsize=30, 
					bbox=dict(facecolor='cyan', alpha=0.5))
	axes[1][0].set_ylabel(yLabel)#+"\n(with MRP)")#(average)")
	axes[1][0].text(0.1,0.8,"With RP",transform=axes[1][0].transAxes,fontsize=30, 
					bbox=dict(facecolor='cyan', alpha=0.5))
	axes[1][0].set_xlabel(xLabel)
	axes[0][0].axis(axisLimits)
	
	# without MRP
	galaxiesBulgeToTotal = getBulgeToTotalRatios(allFieldNames, data, galaxyName)
	axes[0][0].set_title(galaxyName + " High Resolution", y=1.04)#galaxyName)#High Resolution
	plotAvgCamera(galaxiesBulgeToTotal, xFieldName, 'mag', axes[0][0])
	galaxiesBulgeToTotal = getBulgeToTotalRatios(allFieldNames, data, galaxyName + "MRP")
	plotAvgCamera(galaxiesBulgeToTotal, xFieldName, 'mag', axes[1][0])
	'''
	# with MRP
	if MRP:
		galaxiesBulgeToTotal = getBulgeToTotalRatios(allFieldNames, data, galaxyName + "MRP")
		axes[0][1].set_title(galaxyName + "MRP")
		plotAllCamera(galaxiesBulgeToTotal, xFieldName, 'mag', axes[0][1])
		plotAvgCamera(galaxiesBulgeToTotal, xFieldName, 'mag', axes[1][1])
	'''
	fig.subplots_adjust(hspace=0.1, wspace=0.1)


def plotAllAndAvgCamera(fieldDescriptions, allFieldNames, data, 
						galaxyName="", figureIndex=1,
						xFieldIndex=3, yFieldIndex=12, MRP=False,
						compType="central"):
	'''
	
	'''
	xFieldName = fieldDescriptions[xFieldIndex][0]
	yFieldName = fieldDescriptions[yFieldIndex][0]
	xLabel = fieldDescriptions[xFieldIndex][2]
	yLabel = fieldDescriptions[yFieldIndex][2]
	axisLimits = [fieldDescriptions[xFieldIndex][3],
				  fieldDescriptions[xFieldIndex][4],
				  fieldDescriptions[yFieldIndex][3],
				  fieldDescriptions[yFieldIndex][4]]
	figureIndex = str(figureIndex) + " " + yFieldName
	if not MRP:
		fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, num=figureIndex)
		axes = [[axes[0]],[axes[1]]]
	else:
		fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, num=figureIndex)
		axes[1][1].set_xlabel(xLabel)
		
	axes[0][0].set_ylabel(yLabel+" (by camera)")
	axes[1][0].set_ylabel(yLabel+" (average)")
	axes[1][0].set_xlabel(xLabel)
	axes[0][0].axis(axisLimits)
	#axes[0][0].yaxis.set_major_locator(MaxNLocator(prune='upper'))
	
	exclusionList = getBulgeErrors(data, galaxyName)
	centralGalaxies = getGalaxies(allFieldNames, data, compType, galaxyName, exclusionList)
	axes[0][0].set_title(galaxyName)
	plotAllCamera(centralGalaxies, xFieldName, yFieldName, axes[0][0])
	plotAvgCamera(centralGalaxies, xFieldName, yFieldName, axes[1][0])
	if yFieldName == "rad":
		axes[0][0].set_yscale("log")
		axes[0][0].get_yaxis().labelpad = -20
		axes[1][0].set_yscale("log")
		axes[1][0].get_yaxis().labelpad = -20
	
	if MRP:
		centralGalaxies = getGalaxies(allFieldNames, data, compType, galaxyName + "MRP", exclusionList)
		axes[0][1].set_title(galaxyName + "MRP")
		plotAllCamera(centralGalaxies, xFieldName, yFieldName, axes[0][1])
		plotAvgCamera(centralGalaxies, xFieldName, yFieldName, axes[1][1])
		if yFieldName == "rad":
			axes[0][1].set_yscale("log")
			axes[0][1].get_yaxis().labelpad = -20
			axes[1][1].set_yscale("log")
			axes[1][1].get_yaxis().labelpad = -20
	
	if yFieldName == "rad":
		axes[0][0].set_yticks([1.0,3.0,10.0])
		axes[0][0].get_yaxis().set_major_formatter(ticker.FormatStrFormatter("%10.1f"))
		
	fig.subplots_adjust(hspace=0.1, wspace=0.1)


def plotAllGalaxiesAvgCamera(fieldDescriptions, allFieldNames, data,
							 galaxyNames=[""], figureIndex=1, 
							 xFieldIndex=3, yFieldIndex=12, MRP=False,
							 compType="central"):
	'''
	
	'''

	
	xFieldName = fieldDescriptions[xFieldIndex][0]
	yFieldName = fieldDescriptions[yFieldIndex][0]
	xLabel = fieldDescriptions[xFieldIndex][2]
	yLabel = fieldDescriptions[yFieldIndex][2]
	axisLimits = [fieldDescriptions[xFieldIndex][3],
				  fieldDescriptions[xFieldIndex][4],
				  fieldDescriptions[yFieldIndex][3],
				  fieldDescriptions[yFieldIndex][4]]
	figureIndex = str(figureIndex) + " " + yFieldName
	
	if not MRP:
		fig, axes = plt.subplots(1, len(galaxyNames), sharex=True, sharey=True, num=figureIndex)
		if len(galaxyNames) == 1:
			axes = [[axes],[]]
		else:
			axes = [axes,[]]
		axes[0][0].set_ylabel(yLabel)
		axes[0][0].axis(axisLimits)
	else:
		fig, axes = plt.subplots(2, len(galaxyNames), sharex=True, sharey=True, num=figureIndex)
		if len(galaxyNames) == 1:
			axes = [[axes[0]],[axes[1]]]
		axes[0][0].set_ylabel(yLabel)#+"\n(without MRP)")#(by camera)")
		axes[0][0].text(0.1,0.8,"No RP",transform=axes[0][0].transAxes,fontsize=30, 
						bbox=dict(facecolor='cyan', alpha=0.5))
		axes[1][0].set_ylabel(yLabel)#+"\n(with MRP)")#(average)")
		axes[1][0].text(0.1,0.8,"With RP",transform=axes[1][0].transAxes,fontsize=30, 
						bbox=dict(facecolor='cyan', alpha=0.5))
		axes[0][0].axis(axisLimits)
	
	for index, galaxyName in enumerate(galaxyNames):
		exclusionList = getBulgeErrors(data, galaxyName)
		centralGalaxies = getGalaxies(allFieldNames, data, compType, galaxyName, exclusionList)
		axes[0][index].set_title(galaxyName + " High Resolution", y=1.04)#galaxyName)#
		plotAvgCamera(centralGalaxies, xFieldName, yFieldName, axes[0][index])
		
		if not MRP:
			axes[0][index].set_xlabel(xLabel)
			if yFieldName == "rad":
				axes[0][index].set_yscale("log")
				axes[0][index].get_yaxis().labelpad = -50
		if MRP:
			centralGalaxies = getGalaxies(allFieldNames, data, compType, galaxyName + "MRP", exclusionList)
			axes[1][index].set_xlabel(xLabel)
			plotAvgCamera(centralGalaxies, xFieldName, yFieldName, axes[1][index])
			if yFieldName == "rad":
				axes[0][index].set_yscale("log")
				axes[0][index].get_yaxis().labelpad = -50
				axes[1][index].set_yscale("log")
				axes[1][index].get_yaxis().labelpad = -50
		
	if yFieldName == "rad":
		axes[0][0].set_yticks([1.0,3.0,10.0])
		axes[0][0].get_yaxis().set_major_formatter(ticker.FormatStrFormatter("%10.1f"))
		
	fig.subplots_adjust(hspace=0.1, wspace=0.1)


def plotAllFieldsAvgCamera(fieldDescriptions, allFieldNames, data, 
						   galaxyName="", figureIndex=1,
						   xFieldIndex=3,yFields=[12], MRP=False,
						   compType="central"):
	'''
	
	'''
	
	xFieldName = fieldDescriptions[xFieldIndex][0]
	xLabel = fieldDescriptions[xFieldIndex][2]
	
	yFieldNames = []
	for yFieldIndex in yFields:
		yFieldNames.append(fieldDescriptions[yFieldIndex][0])
	figureIndex = str(figureIndex) + " " + " ".join(yFieldNames)
	
	if not MRP:
		fig, axes = plt.subplots(1, len(yFields), sharex=True, sharey=False, num=figureIndex)
		if len(yFields) == 1:
			axes = [[axes],[]]
		else:
			axes = [axes,[]]
	else:
		fig, axes = plt.subplots(2, len(yFields), sharex=True, sharey=False, num=figureIndex)
		if len(yFields) == 1:
			axes = [[axes[0]],[axes[1]]]
	
	for index, yFieldIndex in enumerate(yFields):
		
		yFieldName = fieldDescriptions[yFieldIndex][0]
		yLabel = fieldDescriptions[yFieldIndex][2]
		axisLimits = [fieldDescriptions[xFieldIndex][3],
					  fieldDescriptions[xFieldIndex][4],
					  fieldDescriptions[yFieldIndex][3],
					  fieldDescriptions[yFieldIndex][4]]
		
		
		exclusionList = getBulgeErrors(data, galaxyName)
		centralGalaxies = getGalaxies(allFieldNames, data, compType, galaxyName, exclusionList)
		axes[0][index].set_title(galaxyName)
		axes[0][index].axis(axisLimits)
		plotAvgCamera(centralGalaxies, xFieldName, yFieldName, axes[0][index])
		
			
		if not MRP:
			axes[0][index].set_xlabel(xLabel)
			axes[0][index].set_ylabel(yLabel)
			if yFieldName == "rad":
				axes[0][index].set_yscale("log")
				axes[0][index].set_yticks([1.0,3.0,10.0])
				axes[0][index].get_yaxis().set_major_formatter(ticker.FormatStrFormatter("%10.1f"))
				axes[0][index].get_yaxis().labelpad = -20
		else:
			centralGalaxies = getGalaxies(allFieldNames, data, compType, galaxyName + "MRP", exclusionList)
			axes[0][index].set_ylabel(yLabel+" (without MRP)")
			axes[1][index].set_xlabel(xLabel)
			axes[1][index].set_ylabel(yLabel+" (with MRP)")
			axes[1][index].axis(axisLimits)
			plotAvgCamera(centralGalaxies, xFieldName, yFieldName, axes[1][index])
			
			if yFieldName == "rad":
				axes[0][index].set_yscale("log")
				axes[0][index].set_yticks([1.0,3.0,10.0])
				axes[0][index].get_yaxis().set_major_formatter(ticker.FormatStrFormatter("%10.1f"))
				axes[0][index].get_yaxis().labelpad = -20
				axes[1][index].set_yscale("log")
				axes[1][index].set_yticks([1.0,3.0,10.0])
				axes[1][index].get_yaxis().set_major_formatter(ticker.FormatStrFormatter("%10.1f"))
				axes[1][index].get_yaxis().labelpad = -20
	
	fig.subplots_adjust(left=0.05, right=0.95)
	
if __name__ == "__main__":
	
	#define the command line interface with simUtility.py
	usage = ("\n%prog summaryFile [-h help] [options (with '-'|'--' prefix)]")
			
	# used to parse command line arguments
	parser = OptionParser(usage)
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-s","--special", 
					  help="to run special plot, ignore many other options",
					  action="store_true")
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-c","--allCameras", 
					  help="to show all cameras",
					  action="store_true")
	
	# indicate that you want all galaxies average on one plot
	parser.add_option("-g","--allGalaxies", 
					  help="if you want to plot all galaxies avg field value in one figure",
					  action="store_true")
	
	# indicate that you want all galaxies average on one plot
	parser.add_option("-f","--allFields", 
					  help="if you want to plot all fields avg value in one figure",
					  action="store_true")
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-b","--bulge", 
					  help="if running on results of GALFIT bulge (two component) fit",
					  action="store_true")
	
	# indicate that a bulge component was run to produce results
	parser.add_option("-m","--includeMRP", 
					  help="if all galaxy names have corresponding MRP, plot them too",
					  action="store_true")
	
	# pass the list of galaxies
	parser.add_option("-n","--galaxyNames",
					  help="the space separated list of galaxy names to be plotted (must exist in summary file)",
					  dest="galaxyNames",
					  action="callback", callback=vararg_callback, default=[""])
	
	# pass the list of galaxies
	parser.add_option("-y","--yFields",
					  help="the list of y field names to be plotted (see field descriptions variable)",
					  dest="yFields",
					  action="callback", callback=vararg_callback, default=[])
	
	# pass the field name to be plotted on the y axis
	parser.add_option("-x","--xFieldName",
					  help="the field name of the x values (see field descriptions variable)",
					  default="age")
	
	# pass the field name to be plotted on the y axis
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
	  
	pprint.pprint(vars(options))
	print("summary file used: " + args[0])
	
	'''pyl.rc("axes", linewidth=2)
	pyl.rc("axes", labelsize=20)
	pyl.rc("axes", titlesize=30)
	pyl.rc("xtick.major", size=10)
	pyl.rc("xtick.major", width=2)
	pyl.rc("xtick", labelsize=20)
	pyl.rc("ytick.major", size=10)
	pyl.rc("ytick.major", width=2)
	pyl.rc("ytick", labelsize=20)
	pyl.rc("figure.subplot", right=0.95)'''
	
	# list of lists [name, format, text, lower, upper], one for each field in summary file
	# careful changing this list, there are places in the code where name and
	# field index are referred to as they appear now, check before changing
	fieldDescriptions = [('typ','a10','Galaxy Type'),		   # Field index 0
						 ('id','a10','Galaxy ID'),			   # Field index 1
						 ('ts','a10','Time Step (a)'),			# Field index 2
						 ('age','f4','Age (Gyr)',0,8),			   # Field index 3
						 ('red','f4','Redshift (z)',5,0),		   # Field index 4
						 ('cam','i4','Camera Number'),		   # Field index 5
						 ('fil','a10','Filter'),			   # Field index 6
						 ('px','f4','X Position (pixels)',285,315),	   # Field index 7
						 ('py','f4','Y Position (pixels)',285,315),	   # Field index 8
						 ('mag','f4','Magnitude',30,15),			 # Field index 9
						 ('rpix','f4',r"$R_{eff}$ (pixels)",0.5,50),	  # Field index 10
						 ('rad','f4',r"$R_{eff}$ (kpc)",0.5,15),	  # Field index 11
						 ('ser','f4','Sersic Index',0.05,8.5),			# Field index 12
						 ('ba','f4','Axis Ratio',0.05,1.2),# Field index 13
						 ('pa','f4','Position Angle (deg)',-180,180),	# Field index 14
						 ('sky','f4','sky value')]	 # Field index 14
	
	xFieldIndex = None
	yFields = []
	for index, fieldDescription in enumerate(fieldDescriptions):
		if fieldDescription[0] == options.xFieldName:
			xFieldIndex = index
		for yField in options.yFields:
			if fieldDescription[0] == yField:
				yFields.append(index)
	if xFieldIndex == None:
		xFieldIndex = 4
		
	if not yFields:
		yFields = [12]
		
	# convenience array of just the field names
	allFieldNames = [desc[0] for desc in fieldDescriptions]
	
	# dictionary with field names as keys (above) and array of field values as value
	data = getData(args[0], fieldDescriptions)

	# this is for manual use, limited/inconsistent use of other command line flags
	if options.special:
		figureName = "special"
		xFieldIndex = 4
		xFieldName = "red"
		redLow = 1.4
		redHigh = 2.6
		
		[allSimID, allSimTS, allSimSSFR, allSimMassStars] = pyl.loadtxt(
							"moody_sim_parameters_ssfr.txt", delimiter="\t",
							skiprows=1, usecols=[0,1,54,55], unpack=True,
							dtype={"names":["id","timeStep","ssfr","mass"],
							"formats":["a10","a10","f4","f4"]})
		allSimID = [curID.decode("utf-8") for curID in allSimID]
		allSimTS = [curTS.decode("utf-8").split(".")[1] for curTS in allSimTS]
		
			
		#print(simTS)
		#exit()
		'''
		# read and format the visual classifications
		[visID, visTS, visCAM, visClump] = pyl.loadtxt(
							"mozena_visclass_results.cat.txt", delimiter=",",
							skiprows=1, usecols=[0,1,3,58], unpack=True,
							dtype={"names":["id","timeStep","camera","clump"],
							"formats":["a10","a10","i4","a10"]})
		visID = [curID.decode("utf-8") for curID in visID]
		visTS = [curTS.decode("utf-8").split(".")[1] for curTS in visTS]
		visClumps = []
		for curClump in visClump:
			try:
				visClumps.append(float(curClump))
			except:
				visClumps.append(0.0)
		'''
		textstr = "Simulations"
		# the MRP flag here means either do MRP of all galaxies or not
		if options.includeMRP:
			galaxyNames = [galaxyName+"MRP" for galaxyName in options.galaxyNames]
			textstr = textstr + " with MRP"
			# only allowed because no visual classifications for any MRP	
			visID = ["speed_optimization"] 
		else:
			galaxyNames = options.galaxyNames
			textstr = textstr + " without MRP"
			
		textstr = textstr + "\n" + args[0]
		if options.allCameras:
			textstr = textstr + "\n(all cameras)"
		else:
			textstr = textstr + "\n(avg cameras)"
		textstr = textstr + "\nz = " + str(redLow) + " - " + str(redHigh)  

		# define subplots
		fig = plt.figure(figureName + " " + " ".join(galaxyNames))
		'''ax00 = plt.subplot(221, ylabel="Sersic (n)")
		ax01 = plt.subplot(222)
		ax10 = plt.subplot(223, sharex=ax00, xlabel=r"$R_{eff}$ (kpc)", ylabel="Axis Ratio (q)")
		ax11 = plt.subplot(224, sharey=ax10, xlabel="Sersic (n)")
		totalClumps = 0
		totalOther = 0
		'''
		simID = []
		simTS = []
		simSSFR = []
		simMassStars = []
		for curIndex, curID in enumerate(allSimID):
			if curID in galaxyNames:
				simID.append(allSimID[curIndex])
				simTS.append(allSimTS[curIndex])
				simSSFR.append(allSimSSFR[curIndex])
				simMassStars.append(allSimMassStars[curIndex])
				
		plots = []
		formats = getFormats()
		for gIndex, galaxyName in enumerate(galaxyNames):
				
			# get the data for the current galaxy
			curGalaxyComponents = getGalaxies(allFieldNames, data, galaxyName=galaxyName)
			
			xValsByAge = {}
			yValsByAge = {}
			# subdivide plot components
			for componentNumber, curRed in enumerate(curGalaxyComponents[xFieldName]):
				
				# xVals -> simMassStars/rad^1.5
				# yVals -> simSSFR
				# restrict points to the current galaxy 
				# restrict points to the specified redshift range
				for curIndex, curID in enumerate(simID):
				
					# if the galaxy and time step match, populate the plot
					if ((curID == galaxyName) and 
						(simTS[curIndex].strip("0") == curGalaxyComponents["ts"][componentNumber].strip("0"))):
						if curRed in xValsByAge:
							xValsByAge[curRed].append(np.log10(simMassStars[curIndex]/np.power(curGalaxyComponents["rad"][componentNumber],1.5)))
						else:
							xValsByAge[curRed] = [np.log10(simMassStars[curIndex]/np.power(curGalaxyComponents["rad"][componentNumber],1.5))]
							yValsByAge[curRed] = simSSFR[curIndex]
				
			plots.append([galaxyName, xValsByAge, yValsByAge, formats[gIndex]])
						  
		axes = []
		redshifts = []
		axes.append(plt.subplot(2,3,1))
		redshifts.append([0.5,1.0])
		for i in range(5):
			axes.append(plt.subplot(2,3,i+2))
			redshifts.append([redshifts[-1][-1],redshifts[-1][-1]+0.4])

		for aIndex, curAxes in enumerate(axes):
			redLow = redshifts[aIndex][0]
			redHigh = redshifts[aIndex][1]
					

			for [galaxyName, xValsByAge, yValsByAge, fmtStr] in plots:
				xVals = []
				xErr = []
				yVals = []
				for curRed in xValsByAge:
					# only consider galaxies in the defined redshift range
					if (curRed > redLow) and (curRed < redHigh):
						xVals.append(np.mean(xValsByAge[curRed]))
						xErr.append(np.std(xValsByAge[curRed]))
						yVals.append(yValsByAge[curRed])

				curAxes.errorbar(xVals, yVals, xerr=xErr, ecolor='r', linewidth=0.1, fmt=fmtStr, label=galaxyName)
				curAxes.legend(loc="best", prop={"size":10})
				curAxes.text(0.1,0.8,"%.1f" % redLow + "<z<%.1f" % redHigh,transform=curAxes.transAxes, 
					bbox=dict(facecolor='cyan', alpha=0.5))
				curAxes.set_ylim([2,-3])
				curAxes.set_xlim(8.1,11.9)
			
		axes[0].set_ylabel(r"$log sSFR[Gyr^{-1}]$")
		axes[3].set_ylabel(r"$log sSFR[Gyr^{-1}]$")
		axes[3].set_xlabel(r"$log \sum_{1.5}[M_{solar}kpc^{-1.5}]$")
		axes[4].set_xlabel(r"$log \sum_{1.5}[M_{solar}kpc^{-1.5}]$")
		axes[5].set_xlabel(r"$log \sum_{1.5}[M_{solar}kpc^{-1.5}]$")
		'''
			for galaxyName in galaxyNames:
				
				# get the data for the current galaxy
				curGalaxyComponents = getGalaxies(allFieldNames, data, galaxyName=galaxyName)
				print(np.unique(curGalaxyComponents["ts"]))
				
				# define dictionaries to hold field values with age as their key
				serByAge = {}
				radByAge = {}
				baByAge = {}
				serClumpByAge = {}
				radClumpByAge = {}
				baClumpByAge = {}
				
				xValsByAge = {}
				yValsByAge = {}
				# subdivide plot components
				for componentNumber, curRed in enumerate(curGalaxyComponents[xFieldName]):
					
					# xVals -> simMassStars/rad^1.5
					# yVals -> simSSFR
					# restrict points to the current galaxy 
					# restrict points to the specified redshift range
					for curIndex, curID in enumerate(simID):
						
						# only consider galaxies in the defined redshift range
						if (curRed > redLow) and (curRed < redHigh):
							
							# if the galaxy and time step match, populate the plot
							if ((curID == galaxyName) and 
								(simTS[curIndex].strip("0") == curGalaxyComponents["ts"][componentNumber].strip("0"))):
								if curRed in xValsByAge:
									xValsByAge[curRed].append(np.log10(simMassStars[curIndex]/np.power(curGalaxyComponents["rad"][componentNumber],1.5)))
								else:
									xValsByAge[curRed] = [np.log10(simMassStars[curIndex]/np.power(curGalaxyComponents["rad"][componentNumber],1.5))]
									yValsByAge[curRed] = simSSFR[curIndex]
				
				# different formats for different galaxies
				xVals = []
				xErr = []
				yVals = []
				for curRed in xValsByAge:
					xVals.append(np.mean(xValsByAge[curRed]))
					xErr.append(np.std(xValsByAge[curRed]))
					yVals.append(yValsByAge[curRed])

				# use visual classifications of clumpiness
				for curIndex, curID in enumerate(visID):
					
					# only consider galaxies in the defined redshift range
					if (curRed > redLow) and (curRed < redHigh):
						
						# clumpy
						if ((not options.includeMRP) and (curID == galaxyName) and 
							(int(visTS[curIndex]) == int(curGalaxyComponents["ts"][componentNumber])) and 
							(visCAM[curIndex] == curGalaxyComponents["cam"][componentNumber]) and
							(visClumps[curIndex] >= 0.5)):
							if curRed in serClumpByAge:
								serClumpByAge[curRed].append(curGalaxyComponents["ser"][componentNumber])
								radClumpByAge[curRed].append(curGalaxyComponents["rad"][componentNumber])
								baClumpByAge[curRed].append(curGalaxyComponents["ba"][componentNumber])
							else:
								serClumpByAge[curRed] = [curGalaxyComponents["ser"][componentNumber]]
								radClumpByAge[curRed] = [curGalaxyComponents["rad"][componentNumber]]
								baClumpByAge[curRed] = [curGalaxyComponents["ba"][componentNumber]]
						
						# not clumpy
						else:
							if curRed in serClumpByAge:
								serByAge[curRed].append(curGalaxyComponents["ser"][componentNumber])
								radByAge[curRed].append(curGalaxyComponents["rad"][componentNumber])
								baByAge[curRed].append(curGalaxyComponents["ba"][componentNumber])
							else:
								serByAge[curRed] = [curGalaxyComponents["ser"][componentNumber]]
								radByAge[curRed] = [curGalaxyComponents["rad"][componentNumber]]
								baByAge[curRed] = [curGalaxyComponents["ba"][componentNumber]]
				
			# plots all cameras
			if options.allCameras:
				for i in range(len(np.unique(curGalaxyComponents["cam"]))):
					i=1
					print(i)
					serVals = []
					radVals = []
					baVals = []
					serClumpVals = []
					radClumpVals = []
					baClumpVals = []
					for curRed in serByAge:
						try:
							serVals.append(serByAge[curRed][i])
							radVals.append(radByAge[curRed][i])
							baVals.append(baByAge[curRed][i])
						except:
							continue
					for curRed in serClumpByAge:
						try:
							serClumpVals.append(serClumpByAge[curRed][i])
							radClumpVals.append(radClumpByAge[curRed][i])
							baClumpVals.append(baClumpByAge[curRed][i])
						except:
							continue
						  
					# plot all value for current galaxy
					ax00.semilogx(radVals, serVals, 'g.')
					ax00.semilogx(radClumpVals, serClumpVals, 'b^')
					ax10.semilogx(radVals, baVals, 'g.')
					ax10.semilogx(radClumpVals, baClumpVals, 'b^')
					ax11.plot(serVals, baVals, 'g.')
					ax11.plot(serClumpVals, baClumpVals, 'b^')
					#plt.Axes.
					# track number of clumps vs not
					totalClumps = totalClumps + len(serClumpVals)
					totalOther = totalOther + len(serVals)
					break
			else:
				# compute averages
				avgSerByAge = {}
				avgRadByAge = {}
				avgBAByAge = {}
				for age in serByAge:
					avgSerByAge[age] = np.mean(serByAge[age])
					avgRadByAge[age] = np.mean(radByAge[age])
					avgBAByAge[age] = np.mean(baByAge[age])
				avgSerClumpByAge = {}
				avgRadClumpByAge = {}
				avgBAClumpByAge = {}
				for age in serClumpByAge:
					avgSerClumpByAge[age] = np.mean(serClumpByAge[age])
					avgRadClumpByAge[age] = np.mean(radClumpByAge[age])
					avgBAClumpByAge[age] = np.mean(baClumpByAge[age])
				
				# create lists of averages for plots
				ageVals = list(avgSerByAge.keys())
				serVals = list(avgSerByAge.values())
				radVals = list(avgRadByAge.values())
				baVals = list(avgBAByAge.values())
				ageClumpVals = list(avgSerClumpByAge.keys())
				serClumpVals = list(avgSerClumpByAge.values())
				radClumpVals = list(avgRadClumpByAge.values())
				baClumpVals = list(avgBAClumpByAge.values())
				
				# plot all value for current galaxy
				ax00.semilogx(radVals, serVals, 'g.')
				ax00.semilogx(radClumpVals, serClumpVals, 'b^')
				ax10.semilogx(radVals, baVals, 'g.')
				ax10.semilogx(radClumpVals, baClumpVals, 'b^')
				ax11.plot(serVals, baVals, 'g.')
				ax11.plot(serClumpVals, baClumpVals, 'b^')
				#plt.Axes.
				# track number of clumps vs not
				totalClumps = totalClumps + len(serClumpVals)
				totalOther = totalOther + len(serVals)
				
		# adjust top left subplot
		ax00.set_ylim(fieldDescriptions[12][3],
					  fieldDescriptions[12][4])
		ax00.get_yaxis().set_minor_locator(ticker.AutoMinorLocator(2))
		plt.setp(ax00.get_xticklabels(), visible=False)
		
		# adjust top right subplot
		ax01.set_axis_off()
		textstr = textstr + "\n" + "N = " + str(totalClumps+totalOther)
		textstr = textstr + "\n" + r"$\vartriangle \Longrightarrow$ " + str(totalClumps) + " clumpy"
		textstr = textstr + "\n" + r"$\circ \Longrightarrow$ " + str(totalOther) + " not clumpy"
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)# these are matplotlib.patch.Patch properties
		ax01.text(0.05, 0.95, textstr, transform=ax01.transAxes, fontsize=14,
				verticalalignment='top', bbox=props)
		
		# adjust bottom left subplot
		ax10.set_xlim(fieldDescriptions[11][3],
					  fieldDescriptions[11][4])
		ax10.set_xticks([1.0,3.0,10.0])
		ax10.get_xaxis().set_major_formatter(ticker.FormatStrFormatter("%10.1f"))
		
		# adjust bottom right subplot
		ax11.set_xlim(fieldDescriptions[12][3],
					  fieldDescriptions[12][4])
		ax11.get_xaxis().set_minor_locator(ticker.AutoMinorLocator(2))
		ax11.set_ylim(fieldDescriptions[13][3],
					  fieldDescriptions[13][4])
		ax11.set_yticks([0.2,0.4,0.6,0.8,1.0])
		ax11.get_yaxis().set_minor_locator(ticker.AutoMinorLocator(2))
		plt.setp(ax11.get_yticklabels(), visible=False)
		
		# adjust figure and show
		fig.subplots_adjust(hspace=0, wspace=0)
		'''
		plt.show()
		exit()
	
	
	# single component fit results, one plot per field, all galaxies
	if options.allGalaxies and not options.bulge:
		for figureID, yFieldIndex in enumerate(yFields):				
			plotAllGalaxiesAvgCamera(fieldDescriptions, allFieldNames, data, 
									 galaxyNames=options.galaxyNames,
									 figureIndex=" ".join(options.galaxyNames),
									 xFieldIndex=xFieldIndex,
									 yFieldIndex=yFieldIndex,
									 MRP=options.includeMRP,
									 compType=options.componentType)
			# so that they are maximized
			plt.get_current_fig_manager().window.state('zoomed')

	else:		 
		for galaxyNum, galaxyName in enumerate(options.galaxyNames):
			
			# bulge to total ratios, not easily modified
			if options.bulge:
				plotBulgeTotalRatio(fieldDescriptions, allFieldNames, data, 
									galaxyName=galaxyName, 
									figureIndex=galaxyName+" bulge/total", 
									xFieldIndex=xFieldIndex,
									MRP=options.includeMRP)
				# so that they are maximized
				plt.get_current_fig_manager().window.state('zoomed')
				
			# single component fit results, one plot per galaxy, all fields
			elif options.allFields:			  
				plotAllFieldsAvgCamera(fieldDescriptions, allFieldNames, data, 
									   galaxyName=galaxyName,
									   figureIndex=galaxyName,
									   xFieldIndex=xFieldIndex,
									   yFields=yFields,
									   MRP=options.includeMRP)
				
				# so that they are maximized
				plt.get_current_fig_manager().window.state('zoomed')
				
			# single component fit results, one plot per galaxy per field
			else:	
				for figureID, yFieldIndex in enumerate(yFields):  
					plotAllAndAvgCamera(fieldDescriptions, allFieldNames, data, 
										galaxyName=galaxyName, 
										figureIndex=galaxyName, 
										xFieldIndex=xFieldIndex,
										yFieldIndex=yFieldIndex,
										MRP=options.includeMRP,
										compType=options.componentType)
					# so that they are maximized
					plt.get_current_fig_manager().window.state('zoomed')
			
	# display all figures resulting from above calls
	plt.show()


	
