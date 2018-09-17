
#############################
from astropy.coordinates import SkyCoord
# coordinates of the sources A1 and A2 in
# IRAS 16293-2422 (from Persson et al. 2013)
# A1 source
ra1 = "16h32m22.8733s"
dec1 = "-24d28m36.5800s"
a1 = SkyCoord(ra1, dec1, frame='fk5')
# A2 source
ra2 = "16h32m22.8616s"
dec2 = "-24d28m36.6600s"
a2 = SkyCoord(ra2, dec2, frame='fk5')
# B source
ra3 = "16h32m22.61597s"
dec3 = "-24d28m32.496462s"
b = SkyCoord(ra3, dec3, frame='fk5')


##############################




# to read in the data we use the SpectralCube package. The package reads many different cube formats
# i.e. gildas lmv, casa image, 'normal' fits etc
from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt; plt.ion()
import pywcsgrid2 as pwg
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes

data_cont = SpectralCube.read('data/IRAS16293_Band9.fixed.CONTIN_AP.image.pbcor.subim.fits')

noise = data_cont[0,70:120,80:130]
rms = np.sqrt((noise.flatten()**2).mean())



plt.close(1)
fig = plt.figure(1)
ax = pwg.subplot(111, header=data_cont.header)

divider = make_axes_locatable(ax)
cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
fig.add_axes(cax)

# draw image
im = ax.imshow(data_cont[0].value, cmap=plt.cm.gray_r, origin="lower", interpolation="nearest")

nsig = 3; nsjump = 10
im.set_clim(nsig*rms.value, data_cont[0].max().value)

# draw contour
# have to double check the sigma calculation here...
cont = ax.contour(data_cont[0].value, np.arange(nsig*rms.value, data_cont[0].max().value, nsjump*rms.value),
                  colors='#880000', alpha=0.5)

#for col in cont.collections:
#    col.set_linewidth(0.5)

cbar = plt.colorbar(im, cax=cax)


# adjust cbar ticks and and add levels for contour lines
cbar.set_ticks(np.arange(nsig*rms.value, data_cont[0].max().value, nsjump*rms.value))
cbar.add_lines(cont)

cax.set_ylabel(data_cont[0].unit.to_string())

# plot out A1 and A2
ax["fk5"].plot([a1.ra.value], [a1.dec.value], "+", color='g', ms=9, mew=2, zorder=3)
ax["fk5"].plot([a2.ra.value], [a2.dec.value], "+", color='g', ms=9, mew=2, zorder=3)
ax["fk5"].plot([b.ra.value], [b.dec.value], "+", color='g', ms=9, mew=2, zorder=3)

ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")