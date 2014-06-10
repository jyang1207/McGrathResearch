import os
import sys
import re


def read_param(logfile):
	#logfile = f
	#logfile = open('fit.log')
	lines = logfile.readlines()
	

	line1 = lines[7]
	line2 = lines[8]
	line3 = lines[9]
	line4 = lines[10]
	line5 = lines[13]
	line6 = lines[14]
	line7 = lines[11]
	line8 = lines[12]
	infoline = lines[2]
	#print line1
	#print line2
	#print line3
	#print line4
	#print line5
	#print line6
	#print line7
	#print line8
	
	imname = str.split(infoline, " ")
	imagename = imname[7]
	time = str.split(imagename, "_")
	timestep = time[1]
	
	open_dist = open("dist.tmp", "r")
	d = open_dist.readlines()
	open_dist.close()
	distance = int(d[0])
	dist_min = 20							#need to change this value in both rungalfit and readin if you wish to modify it
	#print distance
	
	split = line1.rstrip()
	split1 = " ".join(split.split())
	split2 = split1.replace("(", "")
	split3 = split2.replace(")", "")
	split4 = split3.replace(",", "")
	split5 = str.split(split4, " ")
	split6 = split5[3:]
	#print split6
	[X_fin,Y_fin,mag_fin,rad_fin,sersic_fin,BA_fin, PA_fin] = split6				

	split7 = line2.rstrip()
	split8 = " ".join(split7.split())
	split9 = split8.replace("(", "")
	split10 = split9.replace(")", "")
	split11 = split10.replace(",", "")
	split12 = str.split(split11, " ")
	split13 = split12[1:]
	#print split13
	[X_err,Y_err,mag_err,rad_err,sersic_err,BA_err, PA_err] = split13
	#print X_err
	
	if distance > dist_min: 
		split14 = line3.rstrip()
		split15 = " ".join(split14.split())
		split16 = split15.replace("(", "")
		split17 = split16.replace(")", "")
		split18 = split17.replace(",", "")
		split19 = str.split(split18, " ")
		split20 = split19[3:]
		#print split20
		[X_fin2,Y_fin2,mag_fin2,rad_fin2,sersic_fin2,BA_fin2, PA_fin2] = split20			#note 2 refers to top component
		#print BA_fin2

		split21 = line4.rstrip()
		split22 = " ".join(split21.split())
		split23 = split22.replace("(", "")
		split24 = split23.replace(")", "")
		split25 = split24.replace(",", "")
		split26 = str.split(split25, " ")
		split27 = split26[1:]
		#print split27
		[X_err2,Y_err2,mag_err2,rad_err2,sersic_err2,BA_err2, PA_err2] = split27
		#print X_err 
	
		split28 = line5.rstrip()
		split29 = split28.replace(",", "")
		split30 = str.split(split29, " ")
		chi_square = split30[3]
		ndof = split30[7]
		#print chi_square
		#print ndof

		split31 = line6.rstrip()
		split32 = split31.replace(",", "")
		split33 = str.split(split32, " ")
		chi_square_nu = split33[3]
	else: 
		split34 = line7.rstrip()
		split35 = split34.replace(",", "")
		split36 = str.split(split35, " ")
		chi_square = split36[3]
		ndof = split36[7]
		print chi_square
		print ndof

		split37 = line8.rstrip()
		split38 = split37.replace(",", "")
		split39 = str.split(split38, " ")
		chi_square_nu = split39[3]
		print chi_square_nu
	
	
	data_lower = open("data_lower.txt", "a")
	data_upper = open("data_upper.txt", "a")
	#open_dist = open("dist.tmp", "r")
	#d = open_dist.readlines()
	#open_dist.close()
	#distance = int(d[0])
	#print distance
	
	if distance <= dist_min: 
		data_lower.write(imagename + "*" + "	" + timestep + "	" + X_fin + "	" + X_err + "	" + Y_fin + "	" + Y_err + "	" + mag_fin + "	" + mag_err + "	" + rad_fin + "	" + rad_err + "	" + sersic_fin + "	" + sersic_err + "	" + BA_fin + "	" + BA_err + "	" + PA_fin + "	" + PA_err + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_lower.close()
		data_upper.write(imagename + "*" + "	" + timestep + "	" + X_fin + "	" + X_err + "	" + Y_fin + "	" + Y_err + "	" + mag_fin + "	" + mag_err + "	" + rad_fin + "	" + rad_err + "	" + sersic_fin + "	" + sersic_err + "	" + BA_fin + "	" + BA_err + "	" + PA_fin + "	" + PA_err + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_upper.close()
	else: 
		data_lower.write(imagename + "	" + timestep + "	" + X_fin + "	" + X_err + "	" + Y_fin + "	" + Y_err + "	" + mag_fin + "	" + mag_err + "	" + rad_fin + "	" + rad_err + "	" + sersic_fin + "	" + sersic_err + "	" + BA_fin + "	" + BA_err + "	" + PA_fin + "	" + PA_err + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_lower.close()
		data_upper.write(imagename + "	" + timestep + "	" + X_fin2 + "	" + X_err2 + "	" + Y_fin2 + "	" + Y_err2 + "	" + mag_fin2 + "	" + mag_err2 + "	" + rad_fin2 + "	" + rad_err2 + "	" + sersic_fin2 + "	" + sersic_err2 + "	" + BA_fin2 + "	" + BA_err2 + "	" + PA_fin2 + "	" + PA_err2 + "	" + chi_square + "	" + ndof + "	" + chi_square_nu + "\n")
		data_upper.close()
	logfile.close()
	
if __name__ == "__main__":
		read_param(sys.argv[1:])
		

#NOTE: if the distance is less than the set value we print the same data to data_upper and data_lower. This is when we have run only one component. If the distance is greater than the set value,
#		we seperate the data into its corresponding data file
