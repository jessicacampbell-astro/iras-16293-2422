from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
from astroquery.splatalogue import Splatalogue
import astropy.units as units
import pywcsgrid2 as pwg
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes
from astropy.io.fits import getheader
import pyfits
import numpy as np 
from scipy import constants

def fnu_LSR(v_LSR,nu_rest):
	nu_LSR = (1.-v_LSR/c_kms) * nu_rest
	return nu_LSR

c_kms = constants.c/1E3

#spw="0a"
spw="0b"
#spw="1a"
#spw="1b"
#spw="2a"
#spw="2b"
#spw="3a"
#spw="3b"

#good_y = 135;  good_x = 140
#bad_y = 135;  bad_x = 142

# Acetaldehyde in spw 0a
filename = "IRAS16293_Band9.fixed.rebin.ms.spw"+str(spw)+".image.pbcor.subim.fits"

if spw=="0a":
	nu_lower = 703.249883899
	nu_upper = 704.176545691
	v_LSR = 1.0
	#y = 145;  x = 147
elif spw=="0b":
	nu_lower = 704.176545691
	nu_upper = 705.103207484
	v_LSR = 2.7
	#y = 145;  x = 147
elif spw=="1a":
	nu_lower = 691.228863879
	nu_upper = 692.155525672
	v_LSR = 2.7
	#y = 145;  x = 147
elif spw=="1b":
	nu_lower = 690.302202086
	nu_upper = 691.228863879
	v_LSR = 2.7
	#y = 145;  x = 147
elif spw=="2a":
	nu_lower = 689.429050373
	nu_upper = 690.355712166
	v_LSR = 2.7
	#y = 145;  x = 147
elif spw=="2b":
	nu_lower = 688.50238858
	nu_upper = 689.429050373
	v_LSR = 2.7
	#y = 145;  x = 147
elif spw=="3a":
	nu_lower = 687.42925759
	nu_upper = 688.355919383
	v_LSR = 2.7
	#y = 145;  x = 147
elif spw=="3b":
	nu_lower = 686.502595797
	nu_upper = 687.42925759
	v_LSR = 2.7
	#y = 145;  x = 147

# Splatalogue query for brightest lines
# upper-level energy limit
energy_max = 1000
# lower-limit on Einstein coefficient
Aij_min = -5
lines = Splatalogue.query_lines(nu_lower*units.GHz,nu_upper*units.GHz,chemical_name="CH3OH",energy_max=energy_max,energy_type="eu_k",version="v2.0")[("Species","Chemical Name","Freq-GHz","Freq Err","Resolved QNs","Log<sub>10</sub> (A<sub>ij</sub>)","E_L (cm^-1)","E_L (K)","E_U (cm^-1)","E_U (K)","Linelist")]
lines_LSR = []
for line in lines:
	#if line[0]=="C":
	if "13C" in line[0]:
		continue
	else:
		nu_rest = line[2]
		EU = line[9]
		nu_LSR = fnu_LSR(v_LSR,nu_rest)
		lines_LSR.append(nu_LSR)
		print line


data_dir = "/Users/jessica/Documents/LEAPS2015/Data/"
data = SpectralCube.read(str(data_dir)+str(filename))
#header = pyfits.open(str(data_dir)+str(filename))[0].header
header = getheader(str(data_dir)+str(filename))
image = data.mean(axis=0).value


# plot average of spectra across entire frequency range
fig = plt.figure()
#plt.step(data.spectral_axis.value/1E9, data[:,y,x].value,color="blue")
#plt.step(data.spectral_axis.value/1E9, data[:,good_y,good_x].value,color="purple")
plt.step(data.spectral_axis.value/1E9, data[:,138,142].value,color="blue")

# plot strongest lines
for nu_LSR in lines_LSR:
	plt.axvline(nu_LSR,linestyle="dashed",color="red")

plt.xlabel("Frequency [GHz]")
plt.ylabel("Flux [Jy/beam]")
plt.show()




