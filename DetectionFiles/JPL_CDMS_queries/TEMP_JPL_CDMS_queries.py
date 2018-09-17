from linedb import *
import numpy as np
from scipy import constants


# spw files
files = np.array(["0a","0b","1a","1b","2a","2b","3a","3b"])
# lower-limits of spw frequencies
freq_lower = np.array([703.249883899,704.176545691,691.228863879,690.302202086,689.429050373,688.50238858,687.42925759,686.502595797])
# upper-limit of spw frequencies
freq_upper = np.array([704.176545691,705.103207484,692.155525672,691.228863879,690.355712166,689.429050373,688.355919383,687.42925759])
# molecular species
species_names = ["CH3OH"]
molecule_names = ["Methanol"]
# upper-level energy limit
Eu_max = 1000
# lower-limit on Einstein coefficient
Aij_min = -5.0

# constants and units
h = constants.h
c = constants.c
k = constants.k
cm_K = (h*c/k)*1E2 # cm^(-1) --> K
K_cm = 1./cm_K

i=0
while i<len(species_names):
	# iterate through molecules
	molecule = molecule_names[i]
	species = species_names[i]
	# iterate through spw files
	j=0
	while j<len(files):
		spw = files[j]
		freqLower = freq_lower[j]
		freqUpper = freq_upper[j]

		print molecule, spw, freqLower, freqUpper
		"""
		Arguments:
        fmin   -- the minimum frequency in MHz
        fmax   -- the maximum frequency in MHz
        specie -- the specie name (default All)
        Eu_max -- maximum upper level energy expressed in cm-1 (default None)
        """
		CDMS_lines		= cdms.search( freqLower*1E3, freqUpper*1E3, specie=species, Eu_max=Eu_max*K_cm )
		JPL_lines		= jpl.search(  freqLower*1E3, freqUpper*1E3, specie=species, Eu_max=Eu_max*K_cm )

		print JPL_lines

		j+=1


	i+=1



