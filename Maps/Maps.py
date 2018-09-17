##############################################################################
# 
#                           About this code:
#
#  This program iterates through all expected emission lines of Acetaldehyde,
#  Methanol, and Methyl Formate, to calculate moment 0, 1, and 2 maps, along
#  with noise & column density maps. The resulting images are then saves to FITS 
#  files, where moment 0, 1, and 2 maps are found in the FITS data indices of 0,
#  1, and 2, respectively, the noise map in the data index 3, and the column 
#  density map in index 4. The program also finds the absolute  brightest pixel
#  in the moment 0 map in the region of the "B" source and obtains the spectrum 
#  at that pixel, which is also saved to a separate FITS file. The velocity axis
#  is the zeroth index of this FITS data, and the line intensity in the 1st.
#  The program also outputs subplots of all maps + bright pixel spectrum, as well
#  as histograms of all moment maps, column densities of all transitions, total
#  column densities, and rotational temperatures maps, separated by molecule. 
#
#  ********* several emission lines are corrected for line blending *********
#
#
#  Written by: Jessica Campbell 
#
##############################################################################


###########################################################################
#                              import modules
#
from spectral_cube import SpectralCube
import matplotlib
import numpy as np
import pyfits
from scipy import constants
import astropy.units as units
from astropy.io.fits import update
from astropy.io import fits
from astropy.io.fits import getheader
import matplotlib.pyplot as plt
import pywcsgrid2 as pwg
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.axes import Axes
import matplotlib.gridspec as gridspec
from astropy.coordinates import SkyCoord
from decimal import * 
import pylab
from matplotlib.colors import LogNorm
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=15)
#
###########################################################################

# local-standard-of-rest frequency
def fnu_LSR(v_LSR,nu_rest):
	nu_LSR = (1.-v_LSR/c_kms) * nu_rest
	return nu_LSR
# local-standard-of-rest velocity
def fv_LSR(nu_LSR,nu_rest):
	v_LSR = (1. - nu_LSR/nu_rest) * c_kms
	return v_LSR
# convert linewidth in velocity to frequency
def fsigma_nu(sigma_v,nu_LSR):
	sigma_nu = (sigma_v/c_kms)*nu_LSR
	return sigma_nu
# solid beam angle in steradians
def fBeam(beam_min,beam_maj):
	# Beam in steradians: beam major and minor axes must be in radians
	# Beam in cm^2: beam major and minor axes must be in cm
	Beam = np.pi*beam_maj*beam_min/(4.*np.log(2.))
	return Beam
# conversion factor between flux (Jy) and effective temperature (K)
def fFactor(Beam,nu):
	# rest frequency in Hertz and c in m/s
	# use Boltzmann factor multiplied by 1E26 (J --> Jy)
	factor = (c_ms**2.)/(2.*k_B_Jy*Beam*nu**2.)
	return factor
# column density
def fN(W,d_parsec,Aij,Beam_cm):
	# W: moment 0 in Jy/beam * km/s
	# d_parsec: distance to source in parsecs
	# d_cm: distance to source in cm
	# c_kms: speed of light in km/s
	# Aij: Einstein coefficient in 1/s
	# Beam_cm : area of the beam in cm^2
	d_cm = d_parsec * parsec * 1E2
	num = W * 4.*np.pi*(d_cm**2.)
	den = 1E26 * h*c_kms*Aij*Beam_cm
	N = num / den
	return N
# writing numbers in scientific notation
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

# coordinates of source B
ra_B = "16h32m22.61597s"
dec_B = "-24d28m32.496462s"
B = SkyCoord(ra_B,dec_B,frame='fk5')

# constants
SNR_threshold = 5.0                              # Signal-to-noise cutoff for integrated intensity
linewidth_threshold = 1.5                        # how many standard deviations to integrate over
FWHM_v = 1.9                                     # FWHM of linewidth km/s
sigma_v = FWHM_v /(2.*np.sqrt(2.*np.log(2.)))    # standard deviation of linewidth km/s
d_parsec = 160.                                  # distance in parsecs
parsec = constants.parsec                        # parsec in m
c_ms = constants.c                               # m/s
c_kms = c_ms/1E3                                 # km/s
k_B = constants.k                                # m^2 kg s^-2 K^-1
#k_B_Jy = k_B * 1E26                             # Jy m^2 K^-1
h = constants.h                                  # m^2 kg s^-1
cm_K = (h*c_ms/k_B)*1E2 # cm^(-1) --> K
K_cm = 1./cm_K

# "B" source location [127:162,123:159]
B_ymin = 127                                     # lower-limit of y pixel for "B" source
B_ymax = 162                                     # upper-limit of y pixel for "B" source
B_xmin = 123                                     # lower-limit of x pixel for "B" source
B_xmax = 159                                     # upper-limit of x pixel for "B" source

# frequency ranges for noise maps
noise_0a = np.array([[703.58,703.60]])
noise_0b = np.array([[704.92,704.98]])
noise_1a = np.array([[691.65,691.69]])
noise_1b = np.array([[691.01,691.03]])
noise_2a = np.array([[690.22,690.27],[689.89,689.92],[689.59,689.65],[689.44,689.465]])
noise_2b = np.array([[688.84,688.87]])
noise_3a = np.array([[688.30,688.35],[687.68,687.71],[687.92,687.95],[687.56,687.585]])
noise_3b = np.array([[687.325,687.345],[687.19,687.21],[686.81,686.84],[686.60,686.64]])

# molecular names
molecules = ["Acetaldehyde","Methanol","MethylFormate"]
# chemical formula names
species_names = ["CH3CHO","CH3OH","CH3OCHO"]
# moment maps to calculate
moments = [0,1,2]
# data directory
data_dir = "/Users/jessica/Documents/LEAPS2015/Data/"

### maps for histograms
# moment 0 maps for each molecule
mom0_B_A = []
mom0_B_M = []
mom0_B_MF = []
# moment 1 maps for each molecule
mom1_B_A = []
mom1_B_M =[]
mom1_B_MF = []
# moment 2 maps for each molecule
mom2_B_A = []
mom2_B_M = []
mom2_B_MF = []
# Nthin maps for each molecule
Nthin_B_A = []
Nthin_B_M = []
Nthin_B_MF = []

# EU values for each detected transitions
EU_A = []
EU_M = []
EU_MF = []

### for energy versus moment 1/2 plots
# acetaldehyde
EU_all_A = []
mom1_all_A = []
mom2_all_A = []
# methanol
EU_all_M = []
mom1_all_M = []
mom2_all_M = []
# methyl formate
EU_all_MF = []
mom1_all_MF = []
mom2_all_MF = []

# iterate through molecules 
for i in range(len(molecules)):
	residuals_all = []
	molecule = molecules[i]
	species = species_names[i]
	f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/"+str(molecule)+".txt", "r")
	lines = f.readlines()
	f.close()

	# iterate through spw files
	files = ["0a","0b","1a","1b","2a","2b","3a","3b"]
	'''
	detectionfile = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/Detected"+str(molecule)+".txt", "w=")
	detectionfile.writelines("Chemical Name | Species | Rest Frequency (GHz) | Resolved Quantum Numbers | Aij | E_U (K) | Statistical Weight | SNRmax | Linelist")
	'''
	for filename in files:
		'''
		detectionfile.writelines("\n\n"+str(filename)+"\n")
		'''
		# open data and header
		data_cube = SpectralCube.read(str(data_dir)+"IRAS16293_Band9.fixed.rebin.ms.spw"+str(filename)+".image.pbcor.subim.fits")
		header = pyfits.open(str(data_dir)+"IRAS16293_Band9.fixed.rebin.ms.spw"+str(filename)+".image.pbcor.subim.fits")[0].header
		header_astropy = getheader(str(data_dir)+"IRAS16293_Band9.fixed.rebin.ms.spw"+str(filename)+".image.pbcor.subim.fits")
		data = data_cube.with_spectral_unit(units.GHz)
		if molecule == "Acetaldehyde":
			if filename == "0a":
				int_lower = 3
				int_upper = 48
			elif filename == "0b":
				int_lower = 49
				int_upper = 73
			elif filename =="1a":
				int_lower = 74
				int_upper = 108
			elif filename =="1b":
				int_lower = 109
				int_upper = 141
			elif filename =="2a":
				int_lower = 142
				int_upper = 158
			elif filename =="2b":
				int_lower = 159
				int_upper = 180
			elif filename =="3a":
				int_lower = 181
				int_upper = 208
			elif filename == "3b":
				int_lower = 209
				int_upper = 230
		if molecule == "Methanol":
			if filename == "0a":
				int_lower = 3
				int_upper = 7
			elif filename == "0b":
				int_lower = 8
				int_upper = 14
			elif filename =="1a":
				int_lower = 15
				int_upper = 17
			elif filename =="1b":
				int_lower = 18
				int_upper = 19
			elif filename =="2a":
				int_lower = 20
				int_upper = 23
			elif filename =="2b":
				int_lower = 24
				int_upper = 27
			elif filename =="3a":
				int_lower = 28
				int_upper = 35
			elif filename == "3b":
				int_lower = 36
				int_upper = 42
		if molecule == "MethylFormate":
			if filename == "0a":
				int_lower = 3
				int_upper = 16
			elif filename == "0b":
				int_lower = 17
				int_upper = 37
			elif filename =="1a":
				int_lower = 38
				int_upper = 54
			elif filename =="1b":
				int_lower = 55
				int_upper = 64
			elif filename =="2a":
				int_lower = 65
				int_upper = 78
			elif filename =="2b":
				int_lower = 79
				int_upper = 102
			elif filename =="3a":
				int_lower = 103
				int_upper = 116
			elif filename == "3b":
				int_lower = 117
				int_upper = 126
		# iterate through each rotational transition
		for line in lines[int_lower:int_upper]:
			linesplit = line[:-1].split(" ")
			nu_rest = float(linesplit[2])
			QN = linesplit[3]
			Aij = float(linesplit[4])
			EU_K = float(linesplit[5])
			g = float(linesplit[6])
			linelist = linesplit[-1]

			# calculate local-standard-of-rest velocity and frequency
			if molecule=="Acetaldehyde" or molecule=="Methanol":
				# found via LTE modeling (Mihkel)
				v_LSR = 2.7 # km/s
			elif molecule=="MethylFormate":
				# found via LTE modeling (Mihkel)
				v_LSR = 4.7 # km/s
			else:
				v_LSR = 0.0
				print "-- v_LSR unknown, assuming %.1f km/s" % ( v_LSR )
			nu_LSR = fnu_LSR(v_LSR,nu_rest) # GHz
			# calculate linewidth in frequency
			sigma_nu = fsigma_nu(sigma_v,nu_LSR) # GHz

			# get slab data in GHz
			slab_GHz = data.spectral_slab((nu_LSR-linewidth_threshold*sigma_nu)*units.GHz,(nu_LSR+linewidth_threshold*sigma_nu)*units.GHz)
			# convert from GHz --> km/s by specifying the rest frequency of the transition
			slab_kms = slab_GHz.with_spectral_unit(units.km/units.s, velocity_convention="radio", rest_value=nu_rest*units.GHz)

			########################################################################
			#                                                                      #
			#                            MOMENT 0 MAP                              #
			#                                                                      #
			########################################################################

			# get velocity axis
			slab_kms_axis = slab_kms.spectral_axis.value        # spectral axis in km/s
			# get channel width in velocity units
			if len(slab_kms_axis)==1:
				delta_v = slab_kms_axis[0]                      # km/s
			else:
				delta_v = abs(np.diff(slab_kms_axis))[0]        # km/s
			mom0 = (slab_kms.sum(axis=0).value) * delta_v       # Jy/beam * km/s

			#print delta_v
			
			mom0_spectralcube = slab_kms.moment(order=0).value  # Jy/beam * km/s
			residuals = mom0 - mom0_spectralcube
			residuals_rms = np.sqrt((residuals.mean())**2.+(residuals.std())**2.)
			residuals_all.append(residuals_rms)

			########################################################################
			#                                                                      #
			#                            MOMENT 1 MAP                              #
			#                                                                      #
			########################################################################

			mom1 = np.zeros(shape=(196,196))
			for j in range(len(slab_kms_axis)):
				v = slab_kms_axis[j]
				mom1 += (slab_kms[j,:,:].value)*v
			# normalize moment 1 map
			mom1 *= delta_v/mom0

			#mom1_spectralcube = slab_kms.moment(order=1).value      # km/s

			########################################################################
			#                                                                      #
			#                            MOMENT 2 MAP                              #
			#                                                                      #
			########################################################################

			mom2_temp = np.zeros(shape=(196,196))
			for j in range(len(slab_kms_axis)):
				v = slab_kms_axis[j]
				mom2_temp += abs(slab_kms[j,:,:].value)*(v-mom1)**2.
			mom2 = np.sqrt(delta_v*mom2_temp/abs(mom0))

			#mom2_spectralcube = abs(slab_kms.moment(order=2).value) # km^2/s^2

			########################################################################
			#                                                                      #
			#                         COLUMN DENSITY                               #
			#                                                                      #
			########################################################################

			# get beam sizes
			beam_min_deg = float(header["BMIN"])                # deg
			beam_maj_deg = float(header["BMAJ"])                # deg
			# convert beam sizes from deg --> rad
			beam_maj_rad = beam_maj_deg * (np.pi/180.)          # rad
			beam_min_rad = beam_min_deg * (np.pi/180.)          # rad
			# convert beam sizes from rad --> cm
			beam_maj_cm = beam_maj_rad * (d_parsec*parsec*1E2)  # cm
			beam_min_cm = beam_min_rad * (d_parsec*parsec*1E2)  # cm
			# beam size in arcseconds
			beam_min_arcsec = beam_min_deg * 3600
			beam_maj_arcsec = beam_maj_deg * 3600.
			beam_arcsec = fBeam(beam_min_arcsec,beam_maj_arcsec)
			# beam position
			beam_pos = float(header["BPA"])                     # deg

			# calculate column density
			W = mom0
			# calculate beam solid angle (area) in units of cm
			beam_cm = fBeam(beam_min_cm,beam_maj_cm)            # cm^2
			N = abs(fN(W,d_parsec,Aij,beam_cm))
			# convert from m^(-2) --> cm^(-2)
			N *= 1E-4


			#########################################################################
			#                                                                       #
			#                      BRIGHTEST PIXEL SPECTRUM                         #
			#                                                                       #
			#########################################################################

			# if there are non-zero moment 0 values in the "B" source region, find absolute brightest pixel location
			bright_pixel = False
			if abs(mom0[B_ymin:B_ymax,B_xmin:B_xmax]).max()>0.0:
				bright_pixel = True
				c_mom0 = max(abs(mom0[B_ymin:B_ymax,B_xmin:B_xmax].min()),mom0[B_ymin:B_ymax,B_xmin:B_xmax].max())
				bright_y,bright_x = np.where(abs(mom0[B_ymin:B_ymax,B_xmin:B_xmax]) - c_mom0==0.0)
				bright_x = (bright_x+B_xmin)[0]
				bright_y = (bright_y+B_ymin)[0]
				#if len(bright_x)>1:
				#	bright_x = bright_x[0]
				#	bright_y = bright_y[0]
				#	print "multiple brightest pixels found -- choosing first one"
				#bright_x = 25 + B_xmin
				#bright_y = 22 + B_ymin

				# extract spectral slab -- want more data than what was used in the analysis for visualization
				bright_cube = data.spectral_slab((nu_LSR-linewidth_threshold*5.*sigma_nu)*units.GHz,(nu_LSR+linewidth_threshold*5.*sigma_nu)*units.GHz)
				# convert from frequency to velocity
				bright_cube_kms = bright_cube.with_spectral_unit(units.km/units.s, velocity_convention="radio", rest_value=nu_rest*units.GHz)
				# extract spectrum at known bright pixel location
				bright_spectrum = bright_cube_kms[:,bright_y,bright_x].value
				# extract velocity spectral axis
				bright_spectrum_axis = bright_cube_kms.spectral_axis.value
				spectra_data = np.array([bright_spectrum_axis,bright_spectrum])
				# flip spectrum if backwards in velocity
				if len(bright_spectrum_axis)>1 and np.diff(bright_spectrum_axis)[0]<0.0:
					bright_spectrum_axis = bright_spectrum_axis[::-1]
					bright_spectrum = bright_spectrum[::-1]


			########################################################################
			#                                                                      #
			#                           NOISE MAPS                                 #
			#                                                                      #
			########################################################################

			# for noise maps
			if filename=="0a":
				noise=noise_0a
			elif filename=="0b":
				noise=noise_0b
			elif filename=="1a":
				noise=noise_1a
			elif filename=="1b":
				noise=noise_1b
			elif filename=="2a":
				noise=noise_2a
			elif filename=="2b":
				noise=noise_2b
			elif filename=="3a":
				noise=noise_3a
			elif filename=="3b":
				noise=noise_3b

			slab_rms_all = []
			# iterate through frequency ranges to calculate average baseline rms across each spw file
			for m in range(len(noise)):
				noise_values=noise[m]
				noise_lower=noise_values[0]
				noise_upper=noise_values[1]
				# get rms of spectral slab between two given frequencies
				slab_noise = data.spectral_slab(noise_lower*units.GHz,noise_upper*units.GHz)
				# rms is square-root(mean**2 + standard-deviation**2)
				slab_rms = np.sqrt((slab_noise.mean(axis=0).value)**2.+(slab_noise.std(axis=0).value)**2.)
				# append rms of slab to array and take average of these rms values for each pixel
				slab_rms_all.append(slab_rms)

			slab_rms_all = np.array(slab_rms_all)
			# per-channel rms noise in flux density
			noise = np.mean(slab_rms_all, axis=0)
			# 
			noise_temp = data[:,70:120,80:130]
			rms_temp = np.sqrt((noise_temp.flattened()**2).mean()).value
			# uncertainty in the integrated line intensity
			delta_mom0 = rms_temp * np.sqrt(slab_GHz.shape[0]) * delta_v


			########################################################################
			#                                                                      #
			#                               SAVE                                   #
			#                               MAPS                                   #
			#                                                                      #
			########################################################################

			maps = np.array([mom0,mom1,mom2,N,noise])

			# add some things to the header

			header_astropy.set("molecule",molecule)                      # molecular species name
			header_astropy.set("species",species)                        # chemical formula for molecular species
			header_astropy.set("deltmom0",delta_mom0)                    # uncertainty in integrated line intensity (Jy/beam km/s)
			header_astropy.set("Nchannel",slab_GHz.shape[0])                # number of spectral channels in spectral slab 
			header_astropy.set("QN",QN)                                  # quantum number
			header_astropy.set("nurest",str(nu_rest)+" (GHz)")           # rest frequency
			#header_astropy.set("nurester",str(nu_rest_err) + " (GHz)")  # rest frequency uncertainty (GHz)
			header_astropy.set("nuLSR",str(nu_LSR) + " (GHz)")           # local-standard-of-rest frequency (GHz)
			header_astropy.set("vLSR",str(v_LSR) + " (km/s)")            # local-standard-of-rest velocity (km/s)
			header_astropy.set("veldisp",str(sigma_v) + " (km/s)")       # linewidth in velocity (km/s)
			header_astropy.set("freqdisp",str(sigma_nu) + " (GHz)")      # linewidth in frequency (GHz)
			header_astropy.set("Aij",Aij)                                # Einstein coefficient
			header_astropy.set("g",g)                                    # upper energy degeneracy (i.e., statistical weight)
			#header_astropy.set("EU_cm",str(EU_cm) + " (1/cm)")          # upper energy level (inverse cm)
			header_astropy.set("EU_K",str(EU_K) + " (K)")                # upper energy level (K)
			#header_astropy.set("EL_cm",str(EL_cm) + " (1/cm)")          # lower energy level (inverse cm)
			#header_astropy.set("EL_K",str(EL_K) + " (K)")               # lower energy level (K)
			header_astropy.set("linelist",linelist)                      # linelist from Splatalogue query
			header_astropy.set("xmin",B_xmin)                            # lower-limit pixel along x for "B" source visualization
			header_astropy.set("xmax",B_xmax)                            # upper-limit pixel along x for "B" source visualization
			header_astropy.set("ymin",B_ymin)                            # lower-limit pixel along y for "B" source visualization
			header_astropy.set("ymax",B_ymax)                            # upper-limit pixel along y for "B" source visualization
			# specify which axis contains which image/data
			header_astropy.set("axis0","mom0")
			header_astropy.set("axis1","mom1")
			header_astropy.set("axis2","mom2")
			header_astropy.set("axis3","N_thin")
			header_astropy.set("axis4","rms")
			if bright_pixel==True:
				header_astropy.set("brightx",bright_x)
				header_astropy.set("brighty",bright_y)
				header_spectrum = header_astropy
				header_spectrum.set("axis0","Velocity (km/s)")
				header_spectrum.set("axis1","Flux (Jy/beam)")
			# write maps to FITS files
			#fits.writeto("/Users/jessica/Documents/LEAPS2015/Maps/Data/"+str(molecule)+"/"+str(QN)+".fits",maps,header_astropy,clobber=True)
			# write bright pixel spectrum to FITS file
			#fits.writeto("/Users/jessica/Documents/LEAPS2015/Spectra/BrightPixelSpectra/"+str(molecule)+"/"+str(QN)+".fits",spectra_data,header_spectrum,clobber=True)



			#print molecule, "\n", QN, "\n", EU_K, " (K) \n", mom0[bright_y,bright_x], "Jy/beam km/s \n", N[bright_y,bright_x], " cm^(-2) \n", g, " (g)\n", Aij, " (Aij)\n", (bright_x,bright_y), "\n\n"
			print molecule, filename, nu_rest
			#temp_delta_v.append(delta_v)
		#print filename, np.mean(noise)

	########################################################################
	#                                                                      #
	#                             RESIDUALS                                #
	#                                                                      #
	########################################################################

	#residuals_all = np.array(residuals_all)
	#print molecule, "moment0 residuals:", residuals_all.mean(), "\n"



			########################################################################
			#                                                                      #
			#                              PLOTS                                   #
			#                                                                      #
			########################################################################
#'''
			# set S/N requirement
			mom0[abs(mom0)<(SNR_threshold*delta_mom0)] = 0.0
			#mom0[abs(mom0)<(SNR_threshold*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]))] = 0.0
			mom1[abs(mom0)<(SNR_threshold*delta_mom0)] = 0.0
			#mom1[abs(mom0)<(SNR_threshold*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]))] = 0.0
			mom2[abs(mom0)<(SNR_threshold*delta_mom0)] = 0.0
			#mom2[abs(mom0)<(SNR_threshold*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]))] = 0.0
			N[abs(mom0)<(SNR_threshold*delta_mom0)] = 0.0
			#N[abs(mom0)<(SNR_threshold*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]))] = 0.0
			# only want N where moment map is positive (i.e., ignoring absorption)
			N[mom0<0.0] = 1E-10
			N[mom0==0.0] = 1E-10


			#if abs(mom0[bright_y,bright_x])>0.0 and EU_K==500.735812468 or EU_K==524.452259817: # energies where N ratios change w/ pixel
			# only plot maps for those that have at least one non-zero pixel
			if np.max(abs(mom0[B_ymin:B_ymax,B_xmin:B_xmax]))>0.0:# and molecule=="Acetaldehyde" and QN=="(7-3-5-6)": #and nu_rest==703.2674413:

				#print filename, molecule, species, np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax])
				print molecule, filename, nu_rest

				# correct blended emission lines
				if molecule=="Acetaldehyde":
					if QN=="(221310-3)" or QN=="(2213-9-3)" or QN=="(1912-8-0)" or QN=="(1912-7-0)" or QN=="(1010-0-6)" or QN=="(1010-1-6)" or QN=="(9-6-3-0)" or QN=="(9-6-4-0)" or QN=="(1311-2-3)" or QN=="(1311-3-3)":
						mom0*=0.5
					elif QN=="(9-6-3-2)":
						mom0*=0.03632730598147303
					elif QN=="(1912-7-2)":
						mom0*=0.963672694018527
				elif molecule=="Methanol":
					if QN=="(8-3-1)" or QN=="(8-3-6-1)" or QN=="(9-3-0)" or QN=="(9-3-6-0)" or QN=="(20-2-0)" or QN=="(20-2-19-0)" or QN=="(23-2-0)" or QN=="(23-2-22-0)":
						mom0*=0.5
					elif QN=="(15-8-1)":
						mom0*=0.5024815008538652
					elif QN=="(15-8-7-1)" or QN=="(15-8-8-1)":
						mom0*=0.24875924957306733

				########################################################################
				#                                                                      #
				#                              HISTOGRAMS                              #
				#                                                                      #
				########################################################################

				mom0_array = mom0[B_ymin:B_ymax,B_xmin:B_xmax][mom0[B_ymin:B_ymax,B_xmin:B_xmax]!=0.0]
				mom1_array = mom1[B_ymin:B_ymax,B_xmin:B_xmax][mom1[B_ymin:B_ymax,B_xmin:B_xmax]!=0.0]
				mom2_array = mom2[B_ymin:B_ymax,B_xmin:B_xmax][mom2[B_ymin:B_ymax,B_xmin:B_xmax]!=0.0]
				Nthin_array = N[B_ymin:B_ymax,B_xmin:B_xmax][N[B_ymin:B_ymax,B_xmin:B_xmax]!=0.0]
				if molecule=="Acetaldehyde":
					EU_A.append(EU_K)
					for i in range(len(mom0_array)):
						mom0_B_A.append(mom0_array[i])
						mom1_B_A.append(mom1_array[i])
						mom2_B_A.append(mom2_array[i])
						# for energy versus moment 1/2 plots
						EU_all_A.append(EU_K)
						mom1_all_A.append(mom1_array[i])
						mom2_all_A.append(mom2_array[i])
					if len(Nthin_array)>0:
						for i in range(len(Nthin_array)):
							Nthin_B_A.append(Nthin_array[i])
				elif molecule=="Methanol":
					EU_M.append(EU_K)
					for i in range(len(mom0_array)):
						mom0_B_M.append(mom0_array[i])
						mom1_B_M.append(mom1_array[i])
						mom2_B_M.append(mom2_array[i])
						# for energy versus moment 1/2 plots
						EU_all_M.append(EU_K)
						mom1_all_M.append(mom1_array[i])
						mom2_all_M.append(mom2_array[i])
					if len(Nthin_array)>0:
						for i in range(len(Nthin_array)):
							Nthin_B_M.append(Nthin_array[i])
				elif molecule=="MethylFormate":
					EU_MF.append(EU_K)
					for i in range(len(mom0_array)):
						mom0_B_MF.append(mom0_array[i])
						mom1_B_MF.append(mom1_array[i])
						mom2_B_MF.append(mom2_array[i])
						# for energy versus moment 1/2 plots
						EU_all_MF.append(EU_K)
						mom1_all_MF.append(mom1_array[i])
						mom2_all_MF.append(mom2_array[i])
					if len(Nthin_array)>0:
						for i in range(len(Nthin_array)):
							Nthin_B_MF.append(Nthin_array[i])


				#SNR = (abs(mom0[B_ymin:B_ymax,B_xmin:B_xmax])/(np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]))).max()
				#detectionfile.writelines(str(molecule)+" "+str(species)+" "+str(nu_rest)+" "+str(QN)+" "+str(format_e(Aij))+" "+str(EU_K)+" "+str(g)+" "+str(round(SNR,1))+" "+str(linelist)+"\n")

	#detectionfile.close()

print "Acetaldehyde:", np.min(EU_A), np.max(EU_A), np.median(EU_A)
print "Methanol:", np.min(EU_M), np.max(EU_M), np.median(EU_M)
print "MethylFormate:", np.min(EU_MF), np.max(EU_MF), np.median(EU_MF)



#'''
'''
				fig = plt.figure(figsize=(11,8.5))

				# contours
				if molecule=="Acetaldehyde":
					nsig=5; nsjump= 1
				elif molecule=="Methanol":
					nsig=5; nsjump= 5
				elif molecule=="MethylFormate":
					nsig=5; nsjump= 2
				c_mom0 = max(abs(mom0[B_ymin:B_ymax,B_xmin:B_xmax].min()),mom0[B_ymin:B_ymax,B_xmin:B_xmax].max())
				#levels = np.array((np.arange(nsig*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]),c_mom0,nsjump*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]))*(-1))[::-1].tolist() + np.arange(nsig*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax]),c_mom0,nsjump*np.mean(noise[B_ymin:B_ymax,B_xmin:B_xmax])).tolist())
				levels = np.array((np.arange(nsig*delta_mom0,c_mom0,nsjump*delta_mom0)*(-1))[::-1].tolist() + np.arange(nsig*delta_mom0,c_mom0,nsjump*delta_mom0).tolist())
				if len(levels)==0:
					levels=np.array([0.0,c_mom0])

				# size in AU: 0.5/3600.*(np.pi/180.)*160.*constants.parsec/constants.au
				size_arcsec = 0.0625
				size_pixels = size_arcsec / (abs(header["CDELT2"])*3600.)
				size_AU = size_arcsec/3600.*(np.pi/180.)*d_parsec*constants.parsec/constants.au


				########################################################################################################################
				#                                                                                                                      #
				#                                                  MOMENT 0                                                            #
				#                                                    MAP                                                               #
				#                                                                                                                      #
				########################################################################################################################

				ax1 = pwg.subplot(321,header=header)
				ax1.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
				# plot YSO position
				ax1["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
				# plot bright pixel location
				plt.scatter(bright_x,bright_y,marker="s",edgecolor="black",facecolor="none",s=50,linewidth=1)
				divider = make_axes_locatable(ax1)
				cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
				fig.add_axes(cax)

				# moment 0 map either has exlusive positive values, exclusively negative values, or both
				# each has a different plotting routine for optimal visualization: 
				# oranges for positive values; blues for negative values, reds/white/blues for both
				if mom0[B_ymin:B_ymax,B_xmin:B_xmax].min()<0.0 and mom0[B_ymin:B_ymax,B_xmin:B_xmax].max()>0.0:
					# center image about the zero point (white) and plot positive values in reds with negatives in blues
					im = ax1.imshow(mom0, cmap=plt.cm.bwr, origin="lower", interpolation="nearest", vmin=-c_mom0, vmax=c_mom0)
				elif mom0[B_ymin:B_ymax,B_xmin:B_xmax].min()==0.0 and mom0[B_ymin:B_ymax,B_xmin:B_xmax].max()>0.0:
					# if only positive values, only plot positive values in oranges
					im = ax1.imshow(mom0, cmap=plt.cm.Oranges, origin="lower", interpolation="nearest", vmin=0.0, vmax=c_mom0)
				elif mom0[B_ymin:B_ymax,B_xmin:B_xmax].min()<0.0 and mom0[B_ymin:B_ymax,B_xmin:B_xmax].max()==0.0:
					# if only negative values, only plot negative values in blues
					im = ax1.imshow(mom0, cmap=plt.cm.Blues_r, origin="lower", interpolation="nearest", vmin=-c_mom0, vmax=0.0)

				# set contours
				cont = ax1.contour(mom0, levels, colors='#880000', alpha=0.5)
				# show colourbar
				cbar = plt.colorbar(im, cax=cax)
				# add contour levels to colourbar
				cbar.add_lines(cont)
				# set colourbar units
				cax.set_ylabel("$\mathrm{Jy}{\,}\mathrm{beam}^{-1}{\,}\mathrm{km}{\,}\mathrm{s}^{-1}$")
				# add reference scale bar
				ax1.add_size_bar(size_pixels, r"${0}$".format(int(size_AU))+"${\,}\mathrm{AU}$", loc=8)
				# Beam size
				ax1.add_beam_size(beam_maj_deg/abs(header["CDELT1"]), beam_min_deg/abs(header["CDELT2"]), beam_pos, loc=3)
				ax1.add_compass(loc=1)
				# turn off axis labels
				#ax1.axis["bottom","left"].toggle(all=True, label=False)
				# specicy moment
				#plt.text(0.08,0.85,r"${\mu}_0$",transform=ax1.transAxes)
				plt.text(0.08,0.85,r"$\mathrm{a)}$",transform=ax1.transAxes)
				# set plot limits
				ax1.set_xlim(B_xmin,B_xmax)
				ax1.set_ylim(B_ymin,B_ymax)
				
				########################################################################################################################
				#                                                                                                                      #
				#                                                  MOMENT 1                                                            #
				#                                                    MAP                                                               #
				#                                                                                                                      #
				########################################################################################################################

				ax2 = pwg.subplot(322,header=header)
				ax2.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
				# plot YSO position
				ax2["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
				# plot bright pixel location
				plt.scatter(bright_x,bright_y,marker="s",edgecolor="black",facecolor="none",s=50,linewidth=1)
				divider = make_axes_locatable(ax2)
				cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
				fig.add_axes(cax)

				#im = ax2.imshow(mom1, cmap=plt.cm.OrRd, origin="lower", interpolation="nearest",vmin=(mom1[B_ymin:B_ymax,B_xmin:B_xmax]).min(), vmax=(mom1[B_ymin:B_ymax,B_xmin:B_xmax]).max())
				if molecule=="Acetaldehyde":
					vmin, vmax = 0.0, 5.0
				elif molecule=="MethylFormate" or "Methanol":
					vmin, vmax = 2.0, 7.0
				im = ax2.imshow(mom1, cmap=plt.cm.OrRd, origin="lower", interpolation="nearest",vmin=vmin, vmax=vmax)
				# same contours as moment 0 map
				cont = ax2.contour(mom0, levels, colors='#880000', alpha=0.5)
				# show colourbar
				cbar = plt.colorbar(im, cax=cax)
				# set colourbar units
				cax.set_ylabel("$\mathrm{km}{\,}\mathrm{s}^{-1}$")
				# turn off axis labels
				#ax2.axis["bottom","left"].toggle(all=True, label=False)
				# specify moment
				#plt.text(0.08,0.85,r"${\mu}_1$",transform=ax2.transAxes)
				plt.text(0.08,0.85,r"$\mathrm{b)}$",transform=ax2.transAxes)
				# set plot limits
				ax2.set_xlim(B_xmin,B_xmax)
				ax2.set_ylim(B_ymin,B_ymax)

				########################################################################################################################
				#                                                                                                                      #
				#                                                  MOMENT 2                                                            #
				#                                                    MAP                                                               #
				#                                                                                                                      #
				########################################################################################################################

				ax3 = pwg.subplot(323,header=header)
				ax3.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
				# plot YSO position
				ax3["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
				# plot bright pixel location
				plt.scatter(bright_x,bright_y,marker="s",edgecolor="black",facecolor="none",s=50,linewidth=1)
				divider = make_axes_locatable(ax3)
				cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
				fig.add_axes(cax)

				#im = ax3.imshow(mom2, cmap=plt.cm.OrRd, origin="lower", interpolation="nearest",vmin=(mom2[B_ymin:B_ymax,B_xmin:B_xmax]).min(), vmax=(mom2[B_ymin:B_ymax,B_xmin:B_xmax]).max())
				if molecule=="Acetaldehyde" or molecule=="MethylFormate":
					vmin, vmax = 0.0, 2.0
				elif molecule=="Methanol":
					vmin, vmax = 0.0, 5.0			
				im = ax3.imshow(mom2, cmap=plt.cm.OrRd, origin="lower", interpolation="nearest",vmin=vmin, vmax=vmax)
				# same contours as moment 0 map
				cont = ax3.contour(mom0, levels, colors='#880000', alpha=0.5)
				# show colourbar
				cbar = plt.colorbar(im, cax=cax)
				cax.set_ylabel("$\mathrm{km}{\,}\mathrm{s}^{-1}$")
				# turn off axis labels
				#ax3.axis["bottom","left"].toggle(all=True, label=False)
				#specify moment
				#plt.text(0.08,0.85,r"${\mu}_2$",transform=ax3.transAxes)
				plt.text(0.08,0.85,r"$\mathrm{c)}$",transform=ax3.transAxes)
				# set plot limits
				ax3.set_xlim(B_xmin,B_xmax)
				ax3.set_ylim(B_ymin,B_ymax)

				########################################################################################################################
				#                                                                                                                      #
				#                                               COLUMN DENSITY                                                         #
				#                                                    MAP                                                               #
				#                                                                                                                      #
				########################################################################################################################

				ax4 = pwg.subplot(324,header=header)
				ax4.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
				# plot YSO position
				ax4["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
				# plot bright pixel location
				plt.scatter(bright_x,bright_y,marker="s",edgecolor="black",facecolor="none",s=50,linewidth=1)
				divider = make_axes_locatable(ax4)
				cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
				fig.add_axes(cax)

				im = ax4.imshow(np.log10(N), cmap=plt.cm.Oranges, origin="lower", interpolation="nearest",vmin=1.,vmax=18.)#,vmin=(N[B_ymin:B_ymax,B_xmin:B_xmax]).min(), vmax=(N[B_ymin:B_ymax,B_xmin:B_xmax]).max())

				# set contours
				cont = ax4.contour(mom0, levels, colors='#880000', alpha=0.5)
				# show colourbar
				cbar = plt.colorbar(im, cax=cax)
				# set colourbar units
				cax.set_ylabel(r"$\mathrm{log}_{10}\,N_\mathrm{ul}^\mathrm{thin}{\;}(\mathrm{cm^{-2}})$", fontsize=17)
				# turn off axis labels
				#ax4.axis["bottom","left"].toggle(all=True, label=False)
				# set plot limits
				plt.text(0.08,0.85,r"$\mathrm{d)}$",transform=ax4.transAxes)
				ax4.set_xlim(B_xmin,B_xmax)
				ax4.set_ylim(B_ymin,B_ymax)

				########################################################################################################################
				#                                                                                                                      #
				#                                                     NOISE                                                            #
				#                                                      MAP                                                             #
				#                                                                                                                      #
				########################################################################################################################

				ax5 = pwg.subplot(325,header=header)
				ax5.set_ticklabel_type("delta", center_pixel=(140.839,146.088))
				# plot YSO position
				ax5["fk5"].plot([B.ra.value], [B.dec.value], "*", color="white", ms=9, mew=2, zorder=3, linewidth=1)
				# plot bright pixel location
				plt.scatter(bright_x,bright_y,marker="s",edgecolor="black",facecolor="none",s=50,linewidth=1)
				divider = make_axes_locatable(ax5)
				cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
				fig.add_axes(cax)

				im = ax5.imshow(noise, cmap=plt.cm.Oranges, origin="lower", interpolation="nearest",vmin=(noise[B_ymin:B_ymax,B_xmin:B_xmax]).min(), vmax=(noise[B_ymin:B_ymax,B_xmin:B_xmax]).max())
				# same contours as moment 0 map
				cont = ax5.contour(mom0, levels, colors='#880000', alpha=0.5)
				cbar = plt.colorbar(im, cax=cax)
				cax.set_ylabel("$\mathrm{Jy}{\,}\mathrm{beam}^{-1}$")
				# turn off axis labels
				#ax4.axis["bottom","left"].toggle(all=True, label=False)
				# specify rms
				#plt.text(0.08,0.85,r"$\mathrm{rms}$",transform=ax5.transAxes)
				plt.text(0.08,0.85,r"$\mathrm{e)}$",transform=ax5.transAxes)
				# set plot limits
				ax5.set_xlim(B_xmin,B_xmax)
				ax5.set_ylim(B_ymin,B_ymax)

				########################################################################################################################
				#                                                                                                                      #
				#                                                BRIGHT PIXEL                                                          #
				#                                                  SPECTRUM                                                            #
				#                                                                                                                      #
				########################################################################################################################
				
				ax6 = plt.subplot(326)

				plt.step(bright_spectrum_axis,bright_spectrum,color="black")
				plt.axhline(0,color="black",linestyle="dashed")
				plt.axvline(v_LSR,color="blue")
				# 1 sigma linewidth range
				plt.axvline(v_LSR+1.*sigma_v,linestyle="dashed",color="red")
				plt.axvline(v_LSR-1.*sigma_v,linestyle="dashed",color="red")
				plt.xlabel(r"$v_\mathrm{LSR}{\;}(\mathrm{km{\,}s^{-1}})$", fontsize=15)
				plt.ylabel(r"$F_{\nu}{\;}(\mathrm{Jy}{\,}\mathrm{beam}^{-1})$", fontsize=15)
				plt.xlim(bright_spectrum_axis[0],bright_spectrum_axis[-1])
				plt.text(0.08,0.85,r"$\mathrm{f)}$",transform=ax6.transAxes)

				plt.subplots_adjust(bottom=0.07, top=0.99, left=0.01, right=0.99, hspace=0.23, wspace=0.1)

				plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/"+str(molecule)+"/"+str(molecule)+"_"+str(QN)+".pdf",bbox_inches="tight")
				#plt.show()
				plt.close()

'''





########################################################################################################################
#                                                                                                                      #
#                                                     PLOT                                                             #
#                                                  HISTOGRAMS                                                          #
#                                                                                                                      #
########################################################################################################################

'''
### moment 0 histograms ###
fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)

pylab.hist(mom0_B_A,bins=np.linspace(np.min([np.min(mom0_B_A),np.min(mom0_B_M),np.min(mom0_B_MF)]),np.max([np.max(mom0_B_A),np.max(mom0_B_M),np.max(mom0_B_MF)]),30),histtype="step", color="red", edgecolor="red", alpha=0.3, linewidth=3, label="$\mathrm{CH_3CHO}$")
pylab.hist(mom0_B_M,bins=np.linspace(np.min([np.min(mom0_B_A),np.min(mom0_B_M),np.min(mom0_B_MF)]),np.max([np.max(mom0_B_A),np.max(mom0_B_M),np.max(mom0_B_MF)]),30),histtype="step", color="blue", edgecolor="blue", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OH}$")
pylab.hist(mom0_B_MF,bins=np.linspace(np.min([np.min(mom0_B_A),np.min(mom0_B_M),np.min(mom0_B_MF)]),np.max([np.max(mom0_B_A),np.max(mom0_B_M),np.max(mom0_B_MF)]),30),histtype="step", color="purple", edgecolor="purple", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OCHO}$")

plt.legend(prop={"size":20})
#ax.tick_params(axis='both', which='major', labelsize=13)
ax.tick_params(axis='both', which='major', labelsize=20)
#plt.xlabel("$W{\;}(\mathrm{Jy{\,}beam^{-1}\,km{\,}s^{-1}})$", fontsize=17)
#plt.ylabel("$N$", fontsize=17)
plt.xlabel("$W{\;}(\mathrm{Jy{\,}beam^{-1}\,km{\,}s^{-1}})$", fontsize=20)
plt.ylabel("$N$", fontsize=20)
plt.text(0.08,0.85,r"$\mathrm{a)}$",transform=ax.transAxes, fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Histograms/mom0_hist.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/mom0_hist.pdf", bbox_inches="tight")
plt.show()

### moment 1 histograms ###
fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)

pylab.hist(mom1_B_A,bins=np.linspace(np.min([np.min(mom1_B_A),np.min(mom1_B_M),np.min(mom1_B_MF)]),np.max([np.max(mom1_B_A),np.max(mom1_B_M),np.max(mom1_B_MF)]),30),histtype="step", color="red", edgecolor="red", alpha=0.3, linewidth=3, label="$\mathrm{CH_3CHO}$")
pylab.hist(mom1_B_M,bins=np.linspace(np.min([np.min(mom1_B_A),np.min(mom1_B_M),np.min(mom1_B_MF)]),np.max([np.max(mom1_B_A),np.max(mom1_B_M),np.max(mom1_B_MF)]),30),histtype="step", color="blue", edgecolor="blue", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OH}$")
pylab.hist(mom1_B_MF,bins=np.linspace(np.min([np.min(mom1_B_A),np.min(mom1_B_M),np.min(mom1_B_MF)]),np.max([np.max(mom1_B_A),np.max(mom1_B_M),np.max(mom1_B_MF)]),30),histtype="step", color="purple", edgecolor="purple", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OCHO}$")

plt.axvline(2.7,linestyle="dashed",color="black")
plt.axvline(4.7,linestyle="dashed",color="black")

#plt.legend(prop={"size":20})
#ax.tick_params(axis='both', which='major', labelsize=13)
ax.tick_params(axis='both', which='major', labelsize=20)
#plt.xlabel("$<v>{\;}(\mathrm{km{\,}s^{-1}})$", fontsize=17)
#plt.ylabel("$N$", fontsize=17)
plt.xlabel("$<v>{\;}(\mathrm{km{\,}s^{-1}})$", fontsize=20)
plt.ylabel("$N$", fontsize=20)
plt.text(0.08,0.85,r"$\mathrm{b)}$",transform=ax.transAxes, fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Histograms/mom1_hist.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/mom1_hist.pdf", bbox_inches="tight")
plt.show()

### moment 2 histograms ###
fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)

pylab.hist(mom2_B_A,bins=np.linspace(np.min([np.min(mom2_B_A),np.min(mom2_B_M),np.min(mom2_B_MF)]),np.max([np.max(mom2_B_A),np.max(mom2_B_M),np.max(mom2_B_MF)]),30),histtype="step", color="red", edgecolor="red", alpha=0.3, linewidth=3, label="$\mathrm{CH_3CHO}$")
pylab.hist(mom2_B_M,bins=np.linspace(np.min([np.min(mom2_B_A),np.min(mom2_B_M),np.min(mom2_B_MF)]),np.max([np.max(mom2_B_A),np.max(mom2_B_M),np.max(mom2_B_MF)]),30),histtype="step", color="blue", edgecolor="blue", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OH}$")
pylab.hist(mom2_B_MF,bins=np.linspace(np.min([np.min(mom2_B_A),np.min(mom2_B_M),np.min(mom2_B_MF)]),np.max([np.max(mom2_B_A),np.max(mom2_B_M),np.max(mom2_B_MF)]),30),histtype="step", color="purple", edgecolor="purple", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OCHO}$")

plt.axvline(1.9/(2.*np.sqrt(2.*np.log(2.))),linestyle="dashed", color="black")

#plt.legend(prop={"size":20})
#ax.tick_params(axis='both', which='major', labelsize=13)
ax.tick_params(axis='both', which='major', labelsize=20)
#plt.xlabel("$<v>^{1/2}{\;}(\mathrm{km{\,}s^{-1}})$", fontsize=17)
#plt.ylabel("$N$", fontsize=17)
plt.xlabel("$<v>^{1/2}{\;}(\mathrm{km{\,}s^{-1}})$", fontsize=20)
plt.ylabel("$N$", fontsize=20)
plt.text(0.08,0.85,r"$\mathrm{c)}$",transform=ax.transAxes, fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Histograms/mom2_hist.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/mom2_hist.pdf", bbox_inches="tight")
plt.show()

### Nthin histograms ###

# turn column densities into numpy arrays for masking
Nthin_B_A = np.array(Nthin_B_A)
Nthin_B_A = np.log10(Nthin_B_A[Nthin_B_A>0.01])
Nthin_B_M = np.array(Nthin_B_M)
Nthin_B_M = np.log10(Nthin_B_M[Nthin_B_M>0.01])
Nthin_B_MF = np.array(Nthin_B_MF)
Nthin_B_MF = np.log10(Nthin_B_MF[Nthin_B_MF>0.01])


fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)

pylab.hist(Nthin_B_A,bins=np.linspace(np.min([np.min(Nthin_B_A),np.min(Nthin_B_M),np.min(Nthin_B_MF)]),np.max([np.max(Nthin_B_A),np.max(Nthin_B_M),np.max(Nthin_B_MF)]),30),histtype="step", color="red", edgecolor="red", alpha=0.3, linewidth=3, label="$\mathrm{CH_3CHO}$")
pylab.hist(Nthin_B_M,bins=np.linspace(np.min([np.min(Nthin_B_A),np.min(Nthin_B_M),np.min(Nthin_B_MF)]),np.max([np.max(Nthin_B_A),np.max(Nthin_B_M),np.max(Nthin_B_MF)]),30),histtype="step", color="blue", edgecolor="blue", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OH}$")
pylab.hist(Nthin_B_MF,bins=np.linspace(np.min([np.min(Nthin_B_A),np.min(Nthin_B_M),np.min(Nthin_B_MF)]),np.max([np.max(Nthin_B_A),np.max(Nthin_B_M),np.max(Nthin_B_MF)]),30),histtype="step", color="purple", edgecolor="purple", alpha=0.3, linewidth=3, label="$\mathrm{CH_3OCHO}$")

plt.legend(prop={"size":20})
#ax.tick_params(axis='both', which='major', labelsize=13)
ax.tick_params(axis='both', which='major', labelsize=20)
#plt.xlabel("$\mathrm{log}\,N_u^{\mathrm{thin}}{\;}(\mathrm{cm^{-2}})$", fontsize=17)
#plt.ylabel("$N$", fontsize=17)
plt.xlabel("$\mathrm{log}\,N_u^{\mathrm{thin}}{\;}(\mathrm{cm^{-2}})$", fontsize=20)
plt.ylabel("$N$", fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Histograms/Nthin_hist.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/Nthin_hist.pdf", bbox_inches="tight")
plt.show()
'''


########################################################################################################################
#                                                                                                                      #
#                                               CORRELATION                                                            #
#                                                  PLOTS                                                               #
#                                                                                                                      #
########################################################################################################################

# moment 1 versus energy
#'''
fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)
# plots
plt.plot(mom1_all_A,EU_all_A,"ro",color="red",alpha=0.3,label="$\mathrm{CH_3CHO}$")
plt.plot(mom1_all_M,EU_all_M,"ro",color="blue",alpha=0.3,label="$\mathrm{CH_3OH}$")
plt.plot(mom1_all_MF,EU_all_MF,"ro",color="purple",alpha=0.3,label="$\mathrm{CH_3OCHO}$")
# average moment 1 values
plt.axvline(np.mean(mom1_all_A),linestyle="dashed",color="red")
plt.axvline(np.mean(mom1_all_M),linestyle="dashed",color="blue")
plt.axvline(np.mean(mom1_all_MF),linestyle="dashed",color="purple")
plt.xlabel("$<v>\,(\mathrm{km\,s^{-1}})$",fontsize=20)
plt.ylabel("$E_U\,(\mathrm{K})$",fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.legend(numpoints=1)
plt.text(0.08,0.85,r"$\mathrm{a)}$",transform=ax.transAxes, fontsize=20)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Energy_mom1.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/Energy_mom1.pdf", bbox_inches="tight")
plt.show()


# moment 2 versus energy

fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)
plt.plot(mom2_all_A,EU_all_A,"ro",color="red",alpha=0.3,label="$\mathrm{CH_3CHO}$")
plt.plot(mom2_all_M,EU_all_M,"ro",color="blue",alpha=0.3,label="$\mathrm{CH_3OH}$")
plt.plot(mom2_all_MF,EU_all_MF,"ro",color="purple",alpha=0.3,label="$\mathrm{CH_3OCHO}$")
# average moment 2 values
plt.axvline(np.mean(mom2_all_A),linestyle="dashed",color="red")
plt.axvline(np.mean(mom2_all_M),linestyle="dashed",color="blue")
plt.axvline(np.mean(mom2_all_MF),linestyle="dashed",color="purple")
plt.xlabel("$<v>^{1/2}\,(\mathrm{km\,s^{-1}})$",fontsize=20)
plt.ylabel("$E_U\,(\mathrm{K})$",fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.text(0.08,0.85,r"$\mathrm{b)}$",transform=ax.transAxes, fontsize=20)

#plt.legend(numpoints=1)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/Energy_mom2.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/Energy_mom2.pdf", bbox_inches="tight")
plt.show()


# moment 2 versus moment 1
fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)
# plots
plt.plot(mom1_all_A,mom2_all_A,"ro",color="red",alpha=0.3,label="$\mathrm{CH_3CHO}$")
plt.plot(mom1_all_M,mom2_all_M,"ro",color="blue",alpha=0.3,label="$\mathrm{CH_3OH}$")
plt.plot(mom1_all_MF,mom2_all_MF,"ro",color="purple",alpha=0.3,label="$\mathrm{CH_3OCHO}$")
# average moment 1 values
plt.axvline(np.mean(mom1_all_A),linestyle="dashed",color="red")
plt.axvline(np.mean(mom1_all_M),linestyle="dashed",color="blue")
plt.axvline(np.mean(mom1_all_MF),linestyle="dashed",color="purple")
# average moment 2 values
plt.axhline(np.mean(mom2_all_A),linestyle="dashed",color="red")
plt.axhline(np.mean(mom2_all_M),linestyle="dashed",color="blue")
plt.axhline(np.mean(mom2_all_MF),linestyle="dashed",color="purple")

plt.xlabel("$<v>\,(\mathrm{km\,s^{-1}})$",fontsize=20)
plt.ylabel("$<v>^{1/2}\,(\mathrm{km\,s^{-1}})$",fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.text(0.08,0.85,r"$\mathrm{c)}$",transform=ax.transAxes, fontsize=20)
#plt.legend(numpoints=1)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/Maps/Plots/mom2_mom1.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/mom2_mom1.pdf", bbox_inches="tight")
plt.show()
#'''



