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
species_names = ["CH3CHO","CH3OH","CH3OCHO"]
#species_names = ["CH3CHO"]
molecule_names = ["Acetaldehyde","Methanol","MethylFormate"]
#molecule_names = ["Acetaldehyde"]
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

		"""
		Arguments:
        fmin   -- the minimum frequency in MHz
        fmax   -- the maximum frequency in MHz
        specie -- the specie name (default All)
        Eu_max -- maximum upper level energy expressed in cm-1 (default None)
        """
		JPL_queryLines		= jpl.search(  freqLower*1E3, freqUpper*1E3, specie=species, Eu_max=Eu_max*K_cm )
		JPL_lines = []
		for JPL_line in JPL_queryLines:
			restFreq = JPL_line.frequency/1E3 # MHz --> GHz
			Aij = JPL_line.einstein_coefficient
			g = JPL_line.upper_level.statistical_weight
			Eu = JPL_line.upper_level.energy*cm_K # cm^(-1) --> K
			QN = JPL_line.upper_level.quantum_numbers
			# append to list
			JPL_lines.append([molecule,species,restFreq,QN,Aij,Eu,g])
			print molecule, spw, restFreq, QN


		if molecule=="Methanol":
			CDMS_QueryLines		= cdms.search( freqLower*1E3, freqUpper*1E3, specie="CH3OH, vt=0,1", Eu_max=Eu_max*K_cm )
			CDMS_lines = []

			for CDMS_line in CDMS_QueryLines:
				restFreq = CDMS_line.frequency/1E3 # MHz --> GHz
				Aij = CDMS_line.einstein_coefficient
				g = CDMS_line.upper_level.statistical_weight
				Eu = JPL_line.upper_level.energy*cm_K # cm^(-1) --> K
				QN = CDMS_line.upper_level.quantum_numbers
				# append to list
				CDMS_lines.append([molecule,species,restFreq,QN,Aij,Eu,g])
				print molecule, spw, restFreq, QN

		# write to text files
		#'''
		#if spw=="0a":
		#	f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/"+str(molecule)+".txt", "w+")
		#	f.writelines("Chemical Name | Species | Rest Frequency (GHz) | Resolved Quantum Numbers | Aij | E_U (K) | Statistical Weight | Linelist"+"\n\n")

		#f.writelines("# "+str(spw)+" "+str(freqLower)+" -- "+str(freqUpper)+" (GHz)\n")
		# write JPL results to text file if non-empty
		#if len(JPL_lines)>0:
		#	for JPL_line in JPL_lines:
		#		f.writelines(str(JPL_line[0])+" "+str(JPL_line[1])+" "+str(JPL_line[2])+" ("+str(JPL_line[3])+") "+str(JPL_line[4])+" "+str(JPL_line[5])+" "+str(JPL_line[6])+" JPL\n")
		# write CDMS results to text file if non-empty and if molecule if Methanol
		#if molecule=="Methanol" and len(CDMS_lines)>0:
		#	for CDMS_line in CDMS_lines:
		#		f.writelines(str(CDMS_line[0])+" "+str(CDMS_line[1])+" "+str(CDMS_line[2])+" ("+str(CDMS_line[3])+") "+str(CDMS_line[4])+" "+str(CDMS_line[5])+" "+str(CDMS_line[6])+" CDMS\n")
		#'''

		#print molecule, spw, freqLower, freqUpper

		j+=1

	#f.close()


	i+=1



