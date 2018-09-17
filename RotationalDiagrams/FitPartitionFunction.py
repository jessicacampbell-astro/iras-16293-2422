#####################################################################################
#
#                              About this program:
#
#  This program fits the partition function for Acetaldehyde, Methanol, and Methyl
#  Formate, to obtain the partition function constants for each molecule. The data
#  used for the fit originates the JPL catalogue via Splatalogue.
#
#####################################################################################

##########################################
#             import modules
#
import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
#
##########################################


def fQ(T,a,b):
	return a*T**b
def round_to_1(x):
	return round(x, -int(math.floor(math.log10(x))))

molecules = ["MethylFormate","Acetaldehyde","Methanol"]

### partition function data from the JPL catalogue via Splatalogue
# Acetaldehyde
A_T = np.array([9.375,18.75,37.50,75.0,150.0,225.0,300.0])
A_Q = np.array([270.1167,760.2458,2154.5221,6495.8039,22892.2045,50049.7231,86841.0823])
# Methanol
M_T = np.array([9.375,18.75,37.50,75.0,150.0,225.0,300.0])
M_Q = np.array([19.5433,68.7464,230.2391,731.0698,2437.7654,5267.8635,9473.1198])
# Methyl Formate
MF_T = np.array([9.375,18.75,37.50,75.0,150.0,225.0,300.0])
MF_Q = np.array([720.82,2030.84,5772.42,17548.82,59072.96,121102.02,199602.70])

### results of partition function fitting routine
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



i = 0
# for subplots:
#fig = plt.figure(figsize=(11,8.5))
# for one plot
fig = plt.figure(figsize=(15,5))
ax = fig.add_subplot(111)
while i<len(molecules):
	molecule=molecules[i]
	if molecule=="Acetaldehyde":
		T_array = A_T
		Q_array = A_Q
	elif molecule=="Methanol":
		T_array = M_T
		Q_array = M_Q
	elif molecule=="MethylFormate":
		T_array = MF_T
		Q_array = MF_Q

	# fit partition function
	popt, pcov = curve_fit(fQ, T_array, Q_array)
	# obtain partition function constants
	a = popt[0]
	b = popt[1]
	# aquire fit
	T_array_fit = np.linspace(T_array[0],T_array[-1],100)
	fit = fQ(T_array_fit,a,b)
	# uncertainties
	popt_errors = np.sqrt(pcov.diagonal())
	a_error = popt_errors[0]
	b_error = popt_errors[1]

	# get values and uncertainties rounded properly
	a_error_rounded = round_to_1(a_error)
	b_error_rounded = round_to_1(b_error)
	a_rounded = round(a,len(str(a_error_rounded).split(".")[1]))
	b_rounded = round(b,len(str(b_error_rounded).split(".")[1]))



	# output molecular species name along with partition function constants
	print i, molecule, a, b

	##### SUBPLOTS #####
	'''
	ax = fig.add_subplot(3,1,i)
	# plot original data
	plt.plot(T_array,Q_array,"ro",color="red")
	# plot partition function fit data
	plt.plot(T_array_fit,fit,color="purple")
	# display results on plot
	plt.text(0.05,0.85,"$a={0} \pm \, {1}$".format(a_rounded,a_error_rounded),transform=ax.transAxes,fontsize=15)
	plt.text(0.05,0.7,"$b={0} \pm \, {1}$".format(b_rounded,b_error_rounded),transform=ax.transAxes,fontsize=15)
	if molecule=="Acetaldehyde":
		plt.text(0.95,0.05,"$(a)$",transform=ax.transAxes,fontsize=15)
	elif molecule=="Methanol":
		plt.text(0.95,0.05,"$(b)$",transform=ax.transAxes,fontsize=15)
	elif molecule=="MethylFormate":
		plt.text(0.95,0.05,"$(c)$",transform=ax.transAxes,fontsize=15)
	#plt.title(str(molecule))

	if i==2:
		plt.ylabel("$Q(T) = aT^b$", fontsize=18)
	elif i==0:
		plt.xlabel("$T{\,}(\mathrm{K})$", fontsize=18)

	i+=1
	'''
	##### ONE PLOT #####

	# data
	if molecule=="Acetaldehyde":
		plt.plot(T_array,Q_array,"ro",color="red",label="$CH_3CHO$")
	elif molecule=="Methanol":
		plt.plot(T_array,Q_array,"ro",color="blue",label="$CH_3OH$")
	elif molecule=="MethylFormate":
		plt.plot(T_array,Q_array,"ro",color="purple",label="$CH_3OCHO$")
	# fit
	plt.plot(T_array_fit,fit,color="black")

	plt.xlabel("$T{\,}(\mathrm{K})$", fontsize=20)
	plt.ylabel("$Q(T) = aT^b$", fontsize=20)
	ax.tick_params(axis='both', which='major', labelsize=20)
	plt.legend(loc="upper left",numpoints=1)

	i+=1

#plt.subplots_adjust(hspace=0.5)
#plt.savefig("/Users/jessica/Documents/LEAPS2015/RotationalDiagrams/PartitionFunction.pdf", bbox_inches="tight")
#plt.savefig("/Users/jessica/Documents/LEAPS2015/FinalReport/PartitionFunction.pdf", bbox_inches="tight")
plt.show()






