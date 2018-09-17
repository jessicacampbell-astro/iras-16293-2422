from spectral_cube import SpectralCube
from astropy.io.fits import getheader
import numpy as np 

filename_list = ["0a","0b","1a","1b","2a","2b","3a","3b"]

data_dir = "/Users/jessica/Documents/LEAPS2015/Data/"

rest_freq = []
for filename in filename_list:
	file = "IRAS16293_Band9.fixed.rebin.ms.spw"+str(filename)+".image.pbcor.subim.fits"
	data = SpectralCube.read(str(data_dir)+str(file))
	header = getheader(str(data_dir)+str(file))
	print filename
	rest_freq.append(header["restfrq"]/1E9)
	#print abs(data.spectral_axis.value[0] - data.spectral_axis.value[1])/1E3
	#print len(data.spectral_axis.value)

print np.min(rest_freq), np.max(rest_freq)