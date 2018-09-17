########################################################################################
#
#                                   About this code: 
#
#  This code queries the Splatalogue database (http://www.cv.nrao.edu/php/splat/) for
#  pure rotational transitions of Acetaldehyde (CH3CHO), Methanol (CH3OH), and Methyl
#  Formate (CH3OCHO), to obtain rest frequencies (and uncertainties), resolved Quantum
#  Numbers (QNs), Einstein Coefficients, upper and lower energy states in units of both
#  both inverse cm and Kelvins, upper level degeneracies (statistical weights), and
#  source catalogues (linelists). The results are then saved to text files.
#
#  Written by: Jessica Campbell 
#
########################################################################################

########### ALMA passbands ###########
#      File | passband (GHz)
#      0a 703.249883899 704.176545691
#      0b 704.176545691 705.103207484
#      1a 691.228863879 692.155525672
#      1b 690.302202086 691.228863879
#      2a 689.429050373 690.355712166
#      2b 688.50238858 689.429050373
#      3a 687.42925759 688.355919383
#      3b 686.502595797 687.42925759

# import modules
from astroquery.splatalogue import Splatalogue
import astropy.units as units
import numpy as np 

# spw files
files = np.array(["0a","0b","1a","1b","2a","2b","3a","3b"])
# lower-limits of spw frequencies
freq_lower = np.array([703.249883899,704.176545691,691.228863879,690.302202086,689.429050373,688.50238858,687.42925759,686.502595797])
# upper-limit of spw frequencies
freq_upper = np.array([704.176545691,705.103207484,692.155525672,691.228863879,690.355712166,689.429050373,688.355919383,687.42925759])
# molecular species
species = ["CH3CHO","CH3OH","CH3OCHO"]
molecule_names = ["Acetaldehyde","Methanol","MethylFormate"]
# upper-level energy limit
energy_max = 1000
# lower-limit on Einstein coefficient
Aij_min = -5.0

j=0
# iterate through molecular species for Splatalogue query
while j<len(species):
	# open text file to save query results
	f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/"+str(molecule_names[j])+"Query.txt", "w+")
	f.writelines("Species | Chemical Name | Freq (GHz) +/- error | Resolved QNs | Aij | E_L (1/cm) | E_L (K) | E_U (1/cm) | E_U (K) | Upper Level Degeneracy | Linelist \n")

	# will count QN occurences for each molecular species to make sure there are no duplicate transitions
	QNs = []

	# for each molecular species, iterate through the spw files
	i=0
	# will count total number of emission lines for each molecular species
	emission_lines = 0
	while i<len(freq_lower):
		filename = files[i]
		f_lower = freq_lower[i] # lower-limit of spw frequency axis
		f_upper = freq_upper[i] # upper-limit of spw frequency axis
		# write spw frequency range to text file
		f.writelines("\n"+str(filename)+" ("+str(f_lower)+" -- "+str(f_upper)+" GHz ) \n")
		# query Splatalogue -- restricted by upper energy level
		lines = Splatalogue.query_lines(f_lower*units.GHz,f_upper*units.GHz,chemical_name=species[j],energy_max=energy_max,energy_type="eu_k",show_upper_degeneracy=True,version="v2.0")[("Species","Chemical Name","Freq-GHz","Freq Err","Resolved QNs","Log<sub>10</sub> (A<sub>ij</sub>)","E_L (cm^-1)","E_L (K)","E_U (cm^-1)","E_U (K)","Upper State Degeneracy","Linelist")]
		# sort query results by Quantum Number
		lines.sort("Resolved QNs")
		# iterate through query results
		for line in lines:
			# narrow results down by pure rotational transitions and lower-limit Einstein coefficient
			if (line[0]==str(species[j])+"v=0" or line[0]==str(species[j])+"vt=0") and line[5]>=(Aij_min):
				if line[-1] == "SLAIM" or (molecule_names[j]=="MethylFormate" and line[2]==690.3326 and filename=="2a"):
					pass # skip -- these are duplicates
				else:
					# count emission line
					emission_lines+=1
					# count QN
					QNs.append(line[4])
					# write final result to text file
					if molecule_names[j]=="MethylFormate":
						# write Methyl Formate as MethylFormate in text file -- easier for programming to read the files
						f.writelines(str(line[0])+" MethylFormate "+str(line[2])+" "+str(line[3])+" "+str(line[4])+" "+str(line[5])+" "+str(line[6])+" "+str(line[7])+" "+str(line[8])+" "+str(line[9])+" "+str(line[10])+" "+str(line[11])+"\n")
					else:
						f.writelines(str(line[0])+" "+str(line[1])+" "+str(line[2])+" "+str(line[3])+" "+str(line[4])+" "+str(line[5])+" "+str(line[6])+" "+str(line[7])+" "+str(line[8])+" "+str(line[9])+" "+str(line[10])+" "+str(line[11])+"\n")


		i+=1
	f.close()


	# ensure there are no duplicate QNs for a given molecular species
	duplicates = 0
	for QN in QNs:
		if QNs.count(QN)>1:
			QNs.remove(QN)
			duplicates+=1

	# output molecular species, number of duplicates, and total number of emission lines
	print str(molecule_names[j]), "has", duplicates, "duplicate emission lines."
	print str(molecule_names[j]), "has", emission_lines, "emission lines."

	j+=1



