from os import listdir
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import astropy.units as units
from astropy.io.fits import getheader
from scipy import constants
import numpy as np

#bright_x,bright_y = 141, 146
bright_x,bright_y = 146, 150

spw_files = ["0a","0b","1a","1b","2a","2b","3a","3b"]

# data directory
data_dir = "/Users/jessica/Documents/LEAPS2015/Data/"

fig = plt.figure(figsize=(11,8.5))

rms_all = []
for i in range(len(spw_files)):
	spw = spw_files[i]
	data_cube = SpectralCube.read(str(data_dir)+"IRAS16293_Band9.fixed.rebin.ms.spw"+str(spw)+".image.pbcor.subim.fits")
	data = data_cube.with_spectral_unit(units.GHz)
	spectrum_amp = data[:,bright_y,bright_x].value
	spectrum_freq = data.spectral_axis.value

	#header = getheader(str(data_dir)+"IRAS16293_Band9.fixed.rebin.ms.spw"+str(spw)+".image.pbcor.subim.fits")
	#noise_temp = data[:,70:120,80:130]
	#rms_temp = np.sqrt((noise_temp.flattened()**2).mean()).value
	#rms_all.append(rms_temp)

	# force position
	if spw=="0a":
		j = 1
	if spw=="0b":
		j = 2
	if spw=="1a":
		j = 4
	if spw=="1b":
		j = 3
	if spw=="2a":
		j = 6
	if spw=="2b":
		j = 5
	if spw=="3a":
		j = 0
	if spw=="3b":
		j = 7

	# add data to plot
	ax = fig.add_subplot(4,2,j)
	plt.plot(spectrum_freq,spectrum_amp,color="black")

plt.text(-0.25,2.6,r"$F_{\nu}{\,}(\mathrm{Jy\,beam^{-1}})$",transform=ax.transAxes,rotation="vertical", fontsize=18)
plt.text(1.0,-0.4,r"${\nu}{\,}(\mathrm{GHz})$",transform=ax.transAxes, fontsize=18)
plt.subplots_adjust(hspace=0.25,wspace=0.2,bottom=0.1, top=0.99, left=0.11, right=0.95)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Spectra.pdf", bbox_inches="tight")
plt.show()



