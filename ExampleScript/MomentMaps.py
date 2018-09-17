###########################################################################
#                              import modules                             #
from astropy.coordinates import SkyCoord
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import pyfits
import pywcsgrid2 as pwg
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes
#                                                                         #
###########################################################################

# coordinates of the sources A1 and A2 in
# IRAS 16293-2422 (from Persson et al. 2013)
ra1 = "16h32m22.8733s"
dec1 = "-24d28m36.5800s"
ra2 = "16h32m22.8616s"
dec2 = "-24d28m36.6600s"
a1 = SkyCoord(ra1,dec1,frame='fk5')
a2 = SkyCoord(ra2,dec2,frame='fk5')

# coordinates of source B
# IRAS 16293-2422 (from Favre et al. 2014)
ra = "16h32m22.61s"
dec = "-24d28m32.4s"
b = SkyCoord(ra,dec,frame='fk5')

# set data directory
data_dir = "/Users/jessica/Documents/LEAPS2015/Data/"
data = SpectralCube.read(str(data_dir)+"IRAS16293_Band9.fixed.rebin.ms.spw0a.image.pbcor.subim.fits")
header = pyfits.open('./Data/IRAS16293_Band9.fixed.rebin.ms.spw0a.image.pbcor.subim.fits')[0].header

# convert frequency axis to velocity axis
data = data.with_spectral_unit(u.km/u.s ,velocity_convention='radio')
# get slab between -175 and -162
slab = data.spectral_slab(-175*u.km/u.s, -162*u.km/u.s)#[:,20:65,30:80]

# calculate the per-channel RMS.
noise = data[:,70:120,80:130]
rms = np.sqrt((noise.flattened()**2).mean())

# calculate moment0 map -- integrated flux density
mom0 = slab.moment0()

# following example adapted from
# http://leejjoon.github.io/matplotlib_astronomy_gallery/ic443/ic443.html

plt.close(1)
fig = plt.figure(1)
ax = pwg.subplot(111, header=header)

divider = make_axes_locatable(ax)
cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
fig.add_axes(cax)
# draw figure
im = ax.imshow(mom0.value, cmap=plt.cm.gray_r, origin="lower", interpolation="nearest")
# calculate one-sigma uncertainty
sigma = abs(np.diff(slab.spectral_axis.value))[0]*rms.value*np.sqrt(slab.shape[0])
nsig = 3; nsjump = 1
im.set_clim(nsig*sigma, mom0.max().value)
# draw contour
cont = ax.contour(mom0.value, np.arange(nsig*sigma, mom0.max().value, nsjump*sigma),
                  colors='#880000', alpha=0.5)

#for col in cont.collections:
#    col.set_linewidth(0.5)
cbar = plt.colorbar(im, cax=cax)
# adjust cbar ticks and and add levels for contour lines
cbar.set_ticks(np.arange(nsig*sigma, mom0.max().value, nsjump*sigma))
cbar.add_lines(cont)

cax.set_ylabel(mom0.unit.to_string())

# plot out A1 and A2
ax["fk5"].plot([a1.ra.value], [a1.dec.value], "+", color='g', ms=9, mew=2, zorder=3)
ax["fk5"].plot([a2.ra.value], [a2.dec.value], "+", color='g', ms=9, mew=2, zorder=3)

##### source A1 & A2 #####
#
#ax.set_xlim(33,73)
#ax.set_ylim(27,57)
#ax.add_size_bar(np.int(1./abs(header["CDELT1"])/3600.), r"$1^{\prime}$", loc=8)
#
#########################

####### source B #######
#
ax.set_xlim(130,155)
ax.set_ylim(140,165)
ax.add_size_bar(np.int(0.5/abs(header["CDELT1"])/3600.), r"$0.5^{\prime}$", loc=8)
#
#########################

ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
ax.add_compass(loc=1)
#ax.add_inner_title("Figure 1", loc=2) # not working???

plt.show()
