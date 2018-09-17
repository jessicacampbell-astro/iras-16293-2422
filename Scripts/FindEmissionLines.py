###########################################################################
#                              import modules
#
from astropy.coordinates import SkyCoord
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import pyfits
import pywcsgrid2 as pwg
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes
from scipy import constants
import astropy.units as units
from matplotlib import rc
import matplotlib
#
###########################################################################

# calculate local-standard-of-rest frequency
def fnu_LSR(v_LSR,nu_rest):
	nu_LSR = (1. - v_LSR/constants.c)*nu_rest
	return nu_LSR

data = SpectralCube.read("/Users/jessica/Documents/LEAPS2015/Data/IRAS16293_Band9.fixed.rebin.ms.spw0a.image.pbcor.subim.fits")
# convert data in Hz to GHz
data = data.with_spectral_unit(u.GHz)
freq_axis = data.spectral_axis.value

noise_lower = [687.325,687.19,686.81,686.60]
noise_upper = [687.345,687.21,686.84,686.64]

fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(111)
plt.step(freq_axis,data[:,142,146].value,color="black")
plt.subplots_adjust(top=0.99,bottom=0.08,left=0.05,right=0.98)
plt.xlabel(r"${\nu}{\,}[\mathrm{GHz}]$",fontsize=15)
plt.ylabel(r"$\mathrm{F}_{\nu}{\,}[\mathrm{Jybeam^{-1}{\cdot}kms^{-1}}]$",fontsize=15)
#plt.axvline(fnu_LSR(2.7*1E3,690.3326*1E9)/1E9,color="red")

i=0
while i<len(noise_lower):
	plt.axvline(noise_lower[i],linestyle="dashed",color="red")
	plt.axvline(noise_upper[i],linestyle="dashed",color="red")
	i+=1

plt.xlim(freq_axis[0],freq_axis[-1])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Spectra/FindAcetaldehyde.pdf", bbox_inches="tight")
plt.show()




