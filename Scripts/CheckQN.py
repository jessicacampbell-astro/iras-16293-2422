from os import listdir
from astropy.io.fits import getheader

maps_dir = "/Users/jessica/Documents/LEAPS2015/Maps/Data/"
molecules = ["Acetaldehyde","Methanol","MethylFormate"]

for molecule in molecules:
	for file in listdir(str(maps_dir)+str(molecule)+"/"):
		if file!= ".DS_Store":
			header = getheader(str(maps_dir)+str(molecule)+"/"+str(file))
			QN_filename = file[:-5]
			QN_head = header["QN"]
			if QN_filename != QN_head:
				print QN_filename, QN_head