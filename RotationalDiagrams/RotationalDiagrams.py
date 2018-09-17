#####################################################################################
#
#                              About this program:
#
#  This program fits rotational diagrams for Acetaldehyde, Methanol, and Methyl
#  Formate, to obtain rotational temperatures and total column densities of each.
#  This program also outputs histograms of the rotational temperatures and total
#  column densities separated by molecule. 
#
#           ******** currently restricted to 0 K < T < 2500 K *********
#     (uncertainties too high to get decent temperature/column density plots)
#
#
#  Written by: Jessica Campbell
#
#####################################################################################

#################################################
#                 import modules
#
from os import listdir
import numpy as np
from scipy import constants
from astropy.io.fits import getdata, getheader
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq
import pywcsgrid2 as pwg
from astropy.coordinates import SkyCoord
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes
import pylab
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=15)
#
#################################################

### functions ###
# partition function
def fQ(a,b,T):
	Q = a*T**b
	return Q
def fZ(T,a,b):
	return a*T**b
### for curve_fit ###
def fline(m,x,b):
	# fit negative slope so temperature should be positive (assuming plot decreases)
	y = -m*x + b
	return y
### for leastsq ###
params = [1,1]
def peval(x,params):
	return -params[0]*x + params[1]
def residuals(params,y,x):
	return y - peval(x,params)

# coordinates of source B
ra_B = "16h32m22.61597s"
dec_B = "-24d28m32.496462s"
B = SkyCoord(ra_B,dec_B,frame='fk5')

k = constants.k

### partition function constants from FitPartitionFunction.py
# Acetaldehyde
A_a = 1.63744357165
A_b = 1.90712717084
# Methanol
M_a = 0.130508059349
M_b = 1.96152117672
# Methyl Formate
MF_a = 9.34139806078
MF_b = 1.74797657064

# "B" source location [127:162,123:159]
B_ymin = 127                                     # lower-limit of y pixel for "B" source
B_ymax = 162                                     # upper-limit of y pixel for "B" source
B_xmin = 123                                     # lower-limit of x pixel for "B" source
B_xmax = 159                                     # upper-limit of x pixel for "B" source

# S/N requirement in integrated line intensity (mom0) map
SNR_threshold = 5.

maps_dir = "/Users/jessica/Documents/LEAPS2015/Maps/Data/"
molecules = ["Acetaldehyde","Methanol","MethylFormate"]
#molecules = ["Acetaldehyde"]

# rotational temperatures for each molecule
Trot_B_A = []
Trot_B_M = []
Trot_B_MF = []
# total column densities for each molecule
Ntot_B_A = []
Ntot_B_M = []
Ntot_B_MF = []

# iterate through molecules
for molecule in molecules:
	# initialize a map to contain column densities of transitions found at each pixel
	Ng_map = []
	EU_map = []
	Trot_map = np.zeros(shape=(196,196))
	Ntot_map = np.zeros(shape=(196,196))
	dir = str(maps_dir)+str(molecule)+"/"
	for file in listdir(dir):
		if file!=".DS_Store":
			# for each transition, find high S/N pixels and obtain column density, appending it N_thin
			# each pixel will later be fit with the function that provides the rotational temperature and total column density of the molecule
			QN = file.split(".")[0]
			# open data
			data = getdata(str(dir)+str(file))
			header = getheader(str(dir)+str(file))
			mom0 = data[0]                                     # integrated line intensity in Jy/beam * km/s
			N = data[3]                                        # column density of transition in cm^-2
			# obtain header information
			B_ymin = header["ymin"]                            # lower-limit of y pixel for "B" source
			B_ymax = header["ymax"]                            # upper-limit of y pixel for "B" source
			B_xmin = header["xmin"]                            # lower-limit of x pixel for "B" source
			B_xmax = header["xmax"]                            # upper-limit of x pixel for "B" source
			delta_mom0 = float(header["deltmom0"])             # uncertainty in integrated intensity
			g = float(header["G"])                             # statistical weight (i.e., upper level degeneracy) of tranition
			EU_K = float(header["EU_K"].split(" ")[0])         # upper level energy in Kelvin# get beam sizes
			beam_min_deg = float(header["BMIN"])               # deg
			beam_maj_deg = float(header["BMAJ"]) 
			# beam position
			beam_pos = float(header["BPA"])                    # deg


			# mask low S/N
			N[abs(mom0)<(SNR_threshold*delta_mom0)] = 0.0
			# mask absorption
			N[mom0<0.0] = 0.0

			if molecule=="Methanol":
				# check K quantum number only for Methanol (apparently not straightforward for other molecules)
				if QN.split("-")[1]=="3":
					Ng_map.append(N/g)
					EU_map.append(EU_K*np.ones(shape=(196,196)))
				else:
					#pass
					Ng_map.append(N/g)
					EU_map.append(EU_K*np.ones(shape=(196,196)))
			elif molecule=="Acetaldehyde" or molecule=="MethylFormate":
				Ng_map.append(N/g)
				EU_map.append(EU_K*np.ones(shape=(196,196)))


			print molecule, QN, g

	# define partition function parameters
	if molecule=="Acetaldehyde":
		a = A_a
		b = A_b
	elif molecule=="Methanol":
		a = M_a
		b = M_b
	elif molecule=="MethylFormate":
		a = MF_a
		b = MF_b

	Ng_map = np.array(Ng_map)
	EU_map = np.array(EU_map)

	j = 0
	# iterate through each pixel, remove zero elements, and print if there is more than one non-zero element (to which we can fit a line)
	for y in range(Ng_map.shape[1]):
		for x in range(Ng_map.shape[2]):
			# column densitites across pixel (x,y)
			Ng_temp = Ng_map[:,y,x]
			EU_temp = EU_map[:,y,x]
			# mask zeros in column densities & mask corresponding statistical weights
			Ng_temp_masked = Ng_temp[Ng_temp!=0.0]
			EU_temp_masked = EU_temp[Ng_temp!=0.0]

			if len(Ng_temp_masked)>1: #and molecule=="Acetaldehyde" and j<10:

				if abs(np.diff(EU_temp_masked)).max()<50.: #len(Ng_temp_masked)==2 and np.diff(EU_temp_masked)[0]<50.:
					pass
				else:
					EU_fit = EU_temp_masked
					Ng_fit = Ng_temp_masked

					# fit
					if len(Ng_fit)>2:
						# take natural logarithm of N/g
						Ng_fit = np.log(Ng_fit)

						######## curve_fit ########
						popt, pcov = curve_fit(fline, EU_fit, Ng_fit)
						errors = np.sqrt(np.diag(pcov))
						slope = popt[0]
						slope_error = errors[0]
						offset = popt[1]
						offset_error = errors[1]
						Trot = 1./slope
						Trot_error = slope_error/slope**2.
						Ntot = np.exp(offset) * fZ(Trot,a,b)
						Ntot_error = offset_error * offset * np.exp(offset-1.) * fZ(Trot,a,b)

						if Trot>0.0 and Trot<2500.:
						#if Trot/Trot_error>=5.0:
							Trot_map[y,x] = Trot
							Ntot_map[y,x] = Ntot


						#########################################################################
						#                                                                       #
						#                     plot rotational diagrams                          #
						#                       in "B" source region                            #
						#                                                                       #
						#########################################################################

						'''
						#if x>B_xmin and x<B_xmax and y>B_ymin and y<B_ymax and (x,y)==(11+B_xmin,7+B_ymin):
						if (x-B_xmin,y-B_ymin) == (9,6):

							#print molecule, (x,y), len(Ng_fit), EU_fit

							#print Trot, Trot_error, abs(Trot/Trot_error)

							# plot
							fit_EU = np.linspace(0,1000,100)
							fit_Ng = fline(slope,fit_EU,offset)

							fig = plt.figure()
							ax = fig.add_subplot(111)
							ax.plot(EU_fit,Ng_fit,"ro",color="red")
							ax.plot(fit_EU,fit_Ng,color="black")
							#ax.set_yscale('log')
							#plt.title(str(molecule))
							# latex
							plt.xlabel("$E_{ul}{\,}(\mathrm{K})$",fontsize=20)
							plt.ylabel(r"$\mathrm{log_{10}} \left( N_{ul}^\mathrm{thin} / g \right)$",fontsize=20)
							plt.text(0.68,0.9,"${0} = {1} \pm {2}$".format("T_\mathrm{rot}",int(Trot),int(Trot_error)),transform=ax.transAxes,fontsize=20)
							#plt.text(0.7,0.9,r"$T_\mathrm{rot} = $"+str(int(Trot))+" $\pm$ "+str(int(Trot_error))+"$\mathrm{K}$",transform=ax.transAxes)
							#plt.text(0.7,0.7,r"$\mathrm{log_{10}}N_{tot} = $"+str(int(np.log10(Ntot)))+"$\mathrm{cm^{-2}}$",transform=ax.transAxes)
							# non-latex
							#plt.xlabel("EU (K)")
							#plt.ylabel("Nthin / g")
							#plt.text(0.7,0.9,"Trot = "+str(int(Trot))+" +/- "+str(int(Trot_error)),transform=ax.transAxes)
							#plt.savefig("/Users/jessica/Documents/LEAPS2015/RotationalDiagrams/"+str(molecule)+"/"+str(x-B_xmin)+"_"+str(y-B_ymin)+".pdf", bbox_inches="tight")
							#plt.savefig("/Users/jessica/Desktop/"+str(x-B_xmin)+"_"+str(y-B_ymin)+".pdf", bbox_inches="tight")
							#plt.savefig("/Users/jessica/Desktop/PopDiagram.pdf", bbox_inches="tight")
							plt.show()
							plt.close()

							print molecule, (x,y), Trot, Trot_error, Ntot, Ntot_error

							j+=1
						'''

	'''
	if molecule=="Acetaldehyde":
		Ntot_A = Ntot_map
	elif molecule=="Methanol":
		Ntot_M = Ntot_map
	elif molecule=="MethylFormate":
		Ntot_MF = Ntot_map
	'''

	########################################################################
	#                                                                      #
	#                          FOR HISTOGRAMS                              #
	#                                                                      #
	########################################################################
	#'''
	Trot_array = Trot_map[B_ymin:B_ymax,B_xmin:B_xmax][Trot_map[B_ymin:B_ymax,B_xmin:B_xmax]!=0.0]
	Ntot_array = Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax][Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax]!=0.0]
	if molecule=="Acetaldehyde":
		for i in range(len(Trot_array)):
			Trot_B_A.append(Trot_array[i])
		for i in range(len(Ntot_array)):
			Ntot_B_A.append(Ntot_array[i])
	elif molecule=="Methanol":
		for i in range(len(Trot_array)):
			Trot_B_M.append(Trot_array[i])
		for i in range(len(Ntot_array)):
			Ntot_B_M.append(Ntot_array[i])
	elif molecule=="MethylFormate":
		for i in range(len(Trot_array)):
			Trot_B_MF.append(Trot_array[i])
		for i in range(len(Ntot_array)):
			Ntot_B_MF.append(Ntot_array[i])

	#'''

	#print molecule, "\n", "Trot:", Trot_array.mean(), Trot_array.min(), Trot_array.max(), "\n", "Ntot:", Ntot_array.mean(), Ntot_array.min(), Ntot_array.max()
	#print molecule, "\n", "Trot:", np.median(Trot_array), "\n", "Ntot:", np.median(Ntot_array)
	#'''








	#########################################################################################
	#                                                                                       #
	#                          ROTATIONAL TEMPERATURE MAPS                                  #
	#                                                                                       #
	#########################################################################################
	'''
	fig = plt.figure()
	ax = pwg.subplot(111, header=header)
	# plot YSO position
	ax["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
	divider = make_axes_locatable(ax)
	cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
	fig.add_axes(cax)

	#im = ax.imshow(Trot_map,cmap=plt.cm.OrRd,origin="lower",interpolation="nearest",vmin=((Trot_map[B_ymin:B_ymax,B_xmin:B_xmax])[~np.isnan(Trot_map[B_ymin:B_ymax,B_xmin:B_xmax])]).min(), vmax=((Trot_map[B_ymin:B_ymax,B_xmin:B_xmax])[~np.isnan(Trot_map[B_ymin:B_ymax,B_xmin:B_xmax])]).max())
	im = ax.imshow(Trot_map,cmap=plt.cm.OrRd,origin="lower",interpolation="nearest",vmin=0.0, vmax=2500.)

	if molecule=="Acetaldehyde":
		# Beam size
		ax.add_beam_size(beam_maj_deg/abs(header["CDELT1"]), beam_min_deg/abs(header["CDELT2"]), beam_pos, loc=3)
		# add reference scale bar
		# size in AU: 0.5/3600.*(np.pi/180.)*160.*constants.parsec/constants.au
		d_parsec = 160.
		size_arcsec = 0.0625
		size_pixels = size_arcsec / (abs(header["CDELT2"])*3600.)
		size_AU = size_arcsec/3600.*(np.pi/180.)*d_parsec*constants.parsec/constants.au
		ax.add_size_bar(size_pixels, r"${0}$".format(int(size_AU))+"${\,}\mathrm{AU}$", loc=8)
		plt.text(0.08,0.85,r"$\mathrm{a)}$",transform=ax.transAxes, fontsize=20)
	elif molecule=="Methanol":
		plt.text(0.08,0.85,r"$\mathrm{b)}$",transform=ax.transAxes, fontsize=20)
	elif molecule=="MethylFormate":
		plt.text(0.08,0.85,r"$\mathrm{c)}$",transform=ax.transAxes, fontsize=20)

	ax.set_ticklabel_type("delta", center_pixel=(140.839,146.088))

	cbar = plt.colorbar(im, cax=cax)
	# set colourbar units
	cax.set_ylabel("$T_\mathrm{rot}\,(\mathrm{K})$", fontsize=17)

	ax.set_xlim(B_xmin,B_xmax)
	ax.set_ylim(B_ymin,B_ymax)
	#plt.savefig("/Users/jessica/Documents/LEAPS2015/RotationalDiagrams/Trot/"+str(molecule)+"/"+(molecule)+"_Tmap.pdf", bbox_inches="tight")
	plt.show()
	plt.close()
	'''


	#########################################################################################
	#                                                                                       #
	#                                COLUMN DENSITY MAPS                                    #
	#                                                                                       #
	#########################################################################################
	'''
	fig = plt.figure()
	ax = pwg.subplot(111, header=header)
	# plot YSO position
	ax["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
	divider = make_axes_locatable(ax)
	cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
	fig.add_axes(cax)

	#im = ax.imshow(np.log10(Ntot_map),cmap=plt.cm.OrRd,origin="lower",interpolation="nearest",vmin=15.0,vmax=np.log10(Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax]).max())#,vmin=((Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax])[~np.isnan(Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax])]).min(), vmax=((Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax])[~np.isnan(Ntot_map[B_ymin:B_ymax,B_xmin:B_xmax])]).max())
	im = ax.imshow(np.log10(Ntot_map),cmap=plt.cm.OrRd,origin="lower",interpolation="nearest",vmin=15.0,vmax=22.0)

	if molecule=="Acetaldehyde":
		# Beam size
		ax.add_beam_size(beam_maj_deg/abs(header["CDELT1"]), beam_min_deg/abs(header["CDELT2"]), beam_pos, loc=3)
		# add reference scale bar
		# size in AU: 0.5/3600.*(np.pi/180.)*160.*constants.parsec/constants.au
		d_parsec = 160.
		size_arcsec = 0.0625
		size_pixels = size_arcsec / (abs(header["CDELT2"])*3600.)
		size_AU = size_arcsec/3600.*(np.pi/180.)*d_parsec*constants.parsec/constants.au
		ax.add_size_bar(size_pixels, r"${0}$".format(int(size_AU))+"${\,}\mathrm{AU}$", loc=8)
		plt.text(0.08,0.85,r"$\mathrm{a)}$",transform=ax.transAxes, fontsize=20)
	elif molecule=="Methanol":
		plt.text(0.08,0.85,r"$\mathrm{b)}$",transform=ax.transAxes, fontsize=20)
	elif molecule=="MethylFormate":
		plt.text(0.08,0.85,r"$\mathrm{c)}$",transform=ax.transAxes, fontsize=20)

	ax.set_ticklabel_type("delta", center_pixel=(140.839,146.088))

	cbar = plt.colorbar(im, cax=cax)
	# set colourbar units
	cax.set_ylabel("$\mathrm{log}_{10}\,N_\mathrm{tot}\,(\mathrm{cm^{-2}})$", fontsize=17)

	ax.set_xlim(B_xmin,B_xmax)
	ax.set_ylim(B_ymin,B_ymax)
	#plt.savefig("/Users/jessica/Documents/LEAPS2015/RotationalDiagrams/Ntot/"+str(molecule)+"/"+(molecule)+"_Nmap.pdf", bbox_inches="tight")
	plt.show()
	plt.close()
	'''


#########################################################################
#                                                                       #
#                    ratio of column densities                          #
#                                                                       #
#########################################################################
'''
# Acetaldehyde to Methanol
Ntot_A_M = np.copy(Ntot_A)
Ntot_A_M[Ntot_M==0.0] = 0.0
Ntot_M_A = np.copy(Ntot_M)
Ntot_M_A[Ntot_A==0.0] = 0.0
# Methyl formate to methanol
Ntot_MF_M = np.copy(Ntot_MF)
Ntot_MF_M[Ntot_M==0.0] = 0.0
Ntot_M_MF = np.copy(Ntot_M)
Ntot_M_MF[Ntot_MF==0.0] = 0.0

Ntot_A_to_M = Ntot_A_M / Ntot_M_A
Ntot_MF_to_M = Ntot_MF_M / Ntot_M_MF


###### Acetaldehyde:Methanol ######

fig = plt.figure()
ax = pwg.subplot(111, header=header)
# plot YSO position
ax["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
divider = make_axes_locatable(ax)
cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
fig.add_axes(cax)

im = ax.imshow(np.log10(Ntot_A_to_M),cmap=plt.cm.OrRd,origin="lower",interpolation="nearest",vmin=np.log10(Ntot_A_to_M[B_ymin:B_ymax,B_xmin:B_xmax][~np.isnan(Ntot_A_to_M[B_ymin:B_ymax,B_xmin:B_xmax])]).min(), vmax=np.log10(Ntot_A_to_M[B_ymin:B_ymax,B_xmin:B_xmax][~np.isnan(Ntot_A_to_M[B_ymin:B_ymax,B_xmin:B_xmax])]).max())

# Beam size
ax.add_beam_size(beam_maj_deg/abs(header["CDELT1"]), beam_min_deg/abs(header["CDELT2"]), beam_pos, loc=3)
# add reference scale bar
# size in AU: 0.5/3600.*(np.pi/180.)*160.*constants.parsec/constants.au
d_parsec = 160.
size_arcsec = 0.0625
size_pixels = size_arcsec / (abs(header["CDELT2"])*3600.)
size_AU = size_arcsec/3600.*(np.pi/180.)*d_parsec*constants.parsec/constants.au
ax.add_size_bar(size_pixels, r"${0}$".format(int(size_AU))+"${\,}\mathrm{AU}$", loc=8)
plt.text(0.08,0.85,r"$\mathrm{a)}$",transform=ax.transAxes, fontsize=20)
ax.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
cbar = plt.colorbar(im, cax=cax)
# set colourbar units
cax.set_ylabel("$\mathrm{log}_{10}(N_A / N_M)$", fontsize=17)
ax.set_xlim(B_xmin,B_xmax)
ax.set_ylim(B_ymin,B_ymax)
plt.savefig("/Users/jessica/Documents/LEAPS2015/RotationalDiagrams/Ntot/Ntot_A_to_M.pdf", bbox_inches="tight")
plt.show()
plt.close()

###### MethylFormate:Methanol ######
fig = plt.figure()
ax = pwg.subplot(111, header=header)
# plot YSO position
ax["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
divider = make_axes_locatable(ax)
cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
fig.add_axes(cax)

im = ax.imshow(np.log10(Ntot_MF_to_M),cmap=plt.cm.OrRd,origin="lower",interpolation="nearest",vmin=np.log10(Ntot_MF_to_M[B_ymin:B_ymax,B_xmin:B_xmax][~np.isnan(Ntot_MF_to_M[B_ymin:B_ymax,B_xmin:B_xmax])]).min(), vmax=np.log10(Ntot_MF_to_M[B_ymin:B_ymax,B_xmin:B_xmax][~np.isnan(Ntot_MF_to_M[B_ymin:B_ymax,B_xmin:B_xmax])]).max())

plt.text(0.08,0.85,r"$\mathrm{b)}$",transform=ax.transAxes, fontsize=20)
ax.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
cbar = plt.colorbar(im, cax=cax)
# set colourbar units
cax.set_ylabel("$\mathrm{log}_{10}(N_{MF} / N_M)$", fontsize=17)
ax.set_xlim(B_xmin,B_xmax)
ax.set_ylim(B_ymin,B_ymax)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/RotationalDiagrams/Ntot/Ntot_MF_to_M.pdf", bbox_inches="tight")
plt.show()
plt.close()
'''


########################################################################################################################
#                                                                                                                      #
#                                                     PLOT                                                             #
#                                                  HISTOGRAMS                                                          #
#                                                                                                                      #
########################################################################################################################

'''
### rotational temperature histograms ###

# take log of rotational temperatures
#Trot_B_A = np.log10(Trot_B_A)
#Trot_B_M = np.log10(Trot_B_M)
#Trot_B_MF = np.log10(Trot_B_MF)

fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)

pylab.hist(Trot_B_A,bins=np.linspace(np.min([np.min(Trot_B_A),np.min(Trot_B_M),np.min(Trot_B_MF)]),np.max([np.max(Trot_B_A),np.max(Trot_B_M),np.max(Trot_B_MF)]),20),histtype="step", color="red", edgecolor="red", alpha=0.3, linewidth=3, label="$\mathrm{CH_3CHO}$")
pylab.hist(Trot_B_M,bins=np.linspace(np.min([np.min(Trot_B_A),np.min(Trot_B_M),np.min(Trot_B_MF)]),np.max([np.max(Trot_B_A),np.max(Trot_B_M),np.max(Trot_B_MF)]),20),histtype="step", color="blue", edgecolor="blue", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OH}$")
pylab.hist(Trot_B_MF,bins=np.linspace(np.min([np.min(Trot_B_A),np.min(Trot_B_M),np.min(Trot_B_MF)]),np.max([np.max(Trot_B_A),np.max(Trot_B_M),np.max(Trot_B_MF)]),20),histtype="step", color="purple", edgecolor="purple", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OCHO}$")

#plt.legend(prop={"size":20})
#ax.tick_params(axis='both', which='major', labelsize=13)
#plt.xlabel("$T_\mathrm{rot}{\;}(\mathrm{K})$", fontsize=17)
#plt.ylabel("$N$", fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel("$T_\mathrm{rot}{\;}(\mathrm{K})$", fontsize=20)
plt.ylabel("$N$", fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Histograms/Trot_hist.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/Trot_hist.pdf", bbox_inches="tight")
plt.show()


### column density histograms ###
# take log of column densities

Ntot_B_A = np.array(Ntot_B_A)
Ntot_B_M = np.array(Ntot_B_M)
Ntot_B_MF = np.array(Ntot_B_MF)

Ntot_B_A = np.log10(Ntot_B_A[Ntot_B_A>0.01])
Ntot_B_M = np.log10(Ntot_B_M[Ntot_B_M>0.01])
Ntot_B_MF = np.log10(Ntot_B_MF[Ntot_B_MF>0.01])

fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)

pylab.hist(Ntot_B_A,bins=np.linspace(np.min([np.min(Ntot_B_A),np.min(Ntot_B_M),np.min(Ntot_B_MF)]),np.max([np.max(Ntot_B_A),np.max(Ntot_B_M),np.max(Ntot_B_MF)]),30),histtype="step", color="red", edgecolor="red", alpha=0.3, linewidth=3, label="$\mathrm{CH_3CHO}$")
pylab.hist(Ntot_B_M,bins=np.linspace(np.min([np.min(Ntot_B_A),np.min(Ntot_B_M),np.min(Ntot_B_MF)]),np.max([np.max(Ntot_B_A),np.max(Ntot_B_M),np.max(Ntot_B_MF)]),30),histtype="step", color="blue", edgecolor="blue", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OH}$")
pylab.hist(Ntot_B_MF,bins=np.linspace(np.min([np.min(Ntot_B_A),np.min(Ntot_B_M),np.min(Ntot_B_MF)]),np.max([np.max(Ntot_B_A),np.max(Ntot_B_M),np.max(Ntot_B_MF)]),30),histtype="step", color="purple", edgecolor="purple", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OCHO}$")

#plt.legend(prop={"size":20})
#ax.tick_params(axis='both', which='major', labelsize=13)
#plt.xlabel("$\mathrm{log}_{10}\,N_u^{\mathrm{tot}}{\;}(\mathrm{cm^{-2}})$", fontsize=17)
#plt.ylabel("$N$", fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel("$\mathrm{log}_{10}\,N_u^{\mathrm{tot}}{\;}(\mathrm{cm^{-2}})$", fontsize=20)
plt.ylabel("$N$", fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Histograms/Ntot_hist.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/Ntot_hist.pdf", bbox_inches="tight")
plt.show()
'''








