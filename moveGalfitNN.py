import os
import sys

if __name__ == "__main__":
	if not os.path.isfile(sys.argv[1]):
		print("argument is not a file")
	else:
		filenames = open(sys.argv[1],'r')
		catalogs = filenames.readlines()
		filenames.close()
		
		for gNNFilename in catalogs:
			gNNFile = open(gNNFilename.strip(), 'r')
			gNNLines = gNNFile.readlines()
			gNNFile.close()
			
			for gNNLine in gNNLines:
				if gNNLine.strip() and gNNLine.strip().split()[0] == "B)":
					template = gNNLine.strip().split()[1]
					dest = "_".join(template.split("_")[:-1])+"_result.txt"
					os.system("mv " + gNNFilename.strip() + " " + dest)
					break