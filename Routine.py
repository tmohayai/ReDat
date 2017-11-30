import ReDat
#from pylab import *
import os
import math
import numpy as np 

"""
Author: Tanaz A. Mohayai
Copyright (C) 2016-present by Tanaz A. Mohayai, Illinois Institute of Technology. All rights reserved.
This routine takes the data file from the ReDat module and turns it into the for009 formatted file. The 
root_dummy is a dummy ROOT file name and is used for illustrative purposes, only. 
"""


# sample root file 
root_file_name = "root_dummy.root"
# sample data quality and PID cuts 
tof_min, tof_max, plane_number, usp_min, usp_max = 28, 31, 2, 130, 150

ReDat.redat(root_file_name, tof_min, tof_max, plane_number, usp_min, usp_max)

# below are a list of necessary steps for converting the ReDat output file to for009 formatted file
arr_us = np.loadtxt('for009_US.dat')
arr_us =arr_us[np.lexsort((arr_us[:,0],arr_us[:,4]))]
np.savetxt('for009_US_Sort.dat',arr_us, fmt='%i %i %i %i %i %g %g %g %g %g %g %g %i %i %i %i %i %i %i %i %i %i %i')

arr_ds = np.loadtxt('for009_DS.dat')
arr_ds =arr_ds[np.lexsort((arr_ds[:,0],arr_ds[:,4]))]
np.savetxt('for009_DS_Sort.dat',arr_ds, fmt='%i %i %i %i %i %g %g %g %g %g %g %g %i %i %i %i %i %i %i %i %i %i %i')

infile_us = file('for009_US_Sort.dat')
newopen_us = open('for009_US_Sort_No_nan.dat', 'w')
for line_us in infile_us:
	if 'nan' in line_us:
		continue
	newopen_us.write(line_us)
newopen_us.close()

infile_ds = file('for009_DS_Sort.dat')
newopen_ds = open('for009_DS_Sort_No_nan.dat', 'w')
for line_ds in infile_ds:
    if 'nan' in line_ds:
        continue
    newopen_ds.write(line_ds)
newopen_ds.close()

filenames = [ "for009_US_Sort_No_nan.dat" , "for009_DS_Sort_No_nan.dat"]
with open('for009.dat', 'w') as outfile:
	for fname in filenames:
		with open(fname) as infile_all:
			content = infile_all.read()
			outfile.write(content)
		
