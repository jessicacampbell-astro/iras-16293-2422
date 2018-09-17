
#############################
from astropy.coordinates import SkyCoord
# coordinates of the sources A1 and A2 in
# IRAS 16293-2422 (from Persson et al. 2013)
ra1 = "16h32m22.8733s"
dec1 = "-24d28m36.5800s"
ra2 = "16h32m22.8616s"
dec2 = "-24d28m36.6600s"

a1 = SkyCoord(ra1,dec1,frame='fk5')
a2 = SkyCoord(ra2,dec2,frame='fk5')

##############################




# to read in the data we use the SpectralCube package. The package reads many different cube formats
# i.e. gildas lmv, casa image, 'normal' fits etc
from spectral_cube import SpectralCube
import numpy as np
# to get the unit right, BUNIT needs to be in 'Jy/beam 'and not 'JY/BEAM '(parsing error)
# sed -i 's/JY//BEAM/Jy//beam/g' data/IRAS16293_Band9.fixed.rebin.ms.spw0a.image.pbcor.subim.fits
# for all files... or 'edhead filename.fits' and manually change it.
data_spw0a = SpectralCube.read('data/IRAS16293_Band9.fixed.rebin.ms.spw0a.image.pbcor.subim.fits')

# now we plot the just the sum of the whole array, just a first quick look
import matplotlib.pyplot as plt; plt.ion()
image= data_spw0a.sum(axis=0).value
plt.imshow(image, interpolation='nearest', origin='lower')
# ok, let's extract the spectrum from pixel err say 147, 145
# test to zoom with the controls to see where this is

#plt.step(data_spw0a.spectral_axis, data_spw0a[:,145,147], where='mid') # B source
plt.step(data_spw0a.spectral_axis, data_spw0a[:,30:55,40:70].mean(axis=1).mean(axis=1), where='mid') # A source
# x axis in frequency... I want km/s
import astropy.units as u
data_spw0a = data_spw0a.with_spectral_unit(u.km/u.s ,velocity_convention='radio')
# now we can plot with km/s
plt.step(data_spw0a.spectral_axis.to(u.km/u.s), data_spw0a[:,30:55,40:70].mean(axis=1).mean(axis=1), where='mid')
#  line between -175 and -162, lets cut out a box around it!
slab = data_spw0a.spectral_slab(-175*u.km/u.s, -162*u.km/u.s)[:,20:65,30:80]

# we can also cut out a box where there is not emission and calculate
# the per-channel RMS.
noise = data_spw0a[:,70:120,80:130]
rms = np.sqrt((noise.flattened()**2).mean())

# now we can calculate the moment 0 map for this line! (integrated flux density map)
# this can be buggy at times I think, so be weary
mom0 = slab.moment0()
# now we want to plot this whole thing, with coordinates and stuff

# One way is to use pywcsgrid2 (WCS:world coordinate system)
# following example adapted from
# http://leejjoon.github.io/matplotlib_astronomy_gallery/ic443/ic443.html
import pywcsgrid2 as pwg
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes

# prepare figure & axes
plt.close(1)
fig = plt.figure(1)
ax = pwg.subplot(111, header=slab.header)

divider = make_axes_locatable(ax)
cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
fig.add_axes(cax)

# draw image
im = ax.imshow(mom0.value, cmap=plt.cm.gray_r, origin="lower", interpolation="nearest")

sigma = abs(np.diff(slab.spectral_axis.value)[0]*rms.value*np.sqrt(slab.shape[0]))
nsig = 3; nsjump = 1
im.set_clim(nsig*sigma, mom0.max().value)

# draw contour
# have to double check the sigma calculation here...
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


ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
