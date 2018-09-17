import numpy as np 

##################################################################################################
#                                                                                                #
#                                        Acetaldehyde                                            #
#                                                                                                #
##################################################################################################
'''
f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/DetectedAcetaldehyde.txt", "r")
lines = f.readlines()
f.close()

f = open("/Users/jessica/Documents/LEAPS2015/FinalReport/AcetaldehydeTable.tex", "w+")
f.writelines(r"\begin{deluxetable*}{lcccccc}"+"\n")
f.writelines(r"\tablecolumns{7}"+"\n")
f.writelines(r"\tablecaption{Detected Acetaldehyde Transitions}"+"\n")
f.writelines(r"\tablehead{"+"\n")
f.writelines(r"\colhead{\restfreq} & \colhead{QN} & \colhead{log($A_{ij}$)} & \colhead{$E_U$} & \colhead{g} & \colhead{SNR$_\mathrm{max}$} & \colhead{linelist}"+"\\\\\n")
f.writelines(r"\colhead{(GHz)} & \colhead{} & \colhead{(s$^{-1}$)} & \colhead{(K)} & \colhead{} & \colhead{} & \colhead{}"+"\\\\\n")
f.writelines(r"\colhead{(1)} & \colhead{(2)} & \colhead{(3)} & \colhead{(4)} & \colhead{(5)} & \colhead{(6)} & \colhead{(7)}"+"\n")
f.writelines(r"}"+"\n")
f.writelines(r"\startdata"+"\n")
f.writelines(r"spw 0a"+"\\\\\n")
for line in lines[3:26]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 0b"+"\\\\\n")
for line in lines[29:34]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 1a"+"\\\\\n")
for line in lines[37:43]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 1b"+"\\\\\n")
for line in lines[46:69]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 2a"+"\\\\\n")
for line in lines[72:77]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 2b"+"\\\\\n")
for line in lines[80:95]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 3a"+"\\\\\n")
for line in lines[98:107]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 3b"+"\\\\\n")
for line in lines[110:121]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")

f.writelines(r"\enddata"+"\n")
f.writelines(r"\label{table:acetaldehyde}"+"\n")
f.writelines(r"\end{deluxetable*}"+"\n")
f.close()
'''

##################################################################################################
#                                                                                                #
#                                          Methanol                                              #
#                                                                                                #
##################################################################################################
'''
f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/DetectedMethanol.txt", "r")
lines = f.readlines()
f.close()

f = open("/Users/jessica/Documents/LEAPS2015/FinalReport/MethanolTable.tex", "w+")
f.writelines(r"\begin{deluxetable*}{lcccccc}"+"\n")
f.writelines(r"\tablecolumns{7}"+"\n")
f.writelines(r"\tablecaption{Detected Methanol Transitions}"+"\n")
f.writelines(r"\tablehead{"+"\n")
f.writelines(r"\colhead{\restfreq} & \colhead{QN} & \colhead{log($A_{ij}$)} & \colhead{$E_U$} & \colhead{g} & \colhead{SNR$_\mathrm{max}$} & \colhead{linelist}"+"\\\\\n")
f.writelines(r"\colhead{(GHz)} & \colhead{} & \colhead{(s$^{-1}$)} & \colhead{(K)} & \colhead{} & \colhead{} & \colhead{}"+"\\\\\n")
f.writelines(r"\colhead{(1)} & \colhead{(2)} & \colhead{(3)} & \colhead{(4)} & \colhead{(5)} & \colhead{(6)} & \colhead{(7)}"+"\n")
f.writelines(r"}"+"\n")
f.writelines(r"\startdata"+"\n")
#f.writelines(r"spw 0a"+"\\\\\n")
#for line in lines[3:17]:
#	line = line[:-1].split(" ")
#	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 0b"+"\\\\\n")
for line in lines[6:11]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 1a"+"\\\\\n")
for line in lines[14:16]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 1b"+"\\\\\n")
for line in lines[19:20]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 2a"+"\\\\\n")
for line in lines[23:26]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 2b"+"\\\\\n")
for line in lines[29:31]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 3a"+"\\\\\n")
for line in lines[34:35]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 3b"+"\\\\\n")
for line in lines[38:42]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")

f.writelines(r"\enddata"+"\n")
f.writelines(r"\label{table:methanol}"+"\n")
f.writelines(r"\end{deluxetable*}"+"\n")
f.close()
'''

##################################################################################################
#                                                                                                #
#                                      Methyl Formate                                            #
#                                                                                                #
##################################################################################################
'''
f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/DetectedMethylFormate.txt", "r")
lines = f.readlines()
f.close()

f = open("/Users/jessica/Documents/LEAPS2015/FinalReport/MethylFormateTable.tex", "w+")
f.writelines(r"\begin{deluxetable*}{lcccccc}"+"\n")
f.writelines(r"\tablecolumns{7}"+"\n")
f.writelines(r"\tablecaption{Detected Methyl Formate Transitions}"+"\n")
f.writelines(r"\tablehead{"+"\n")
f.writelines(r"\colhead{\restfreq} & \colhead{QN} & \colhead{log($A_{ij}$)} & \colhead{$E_U$} & \colhead{g} & \colhead{SNR$_\mathrm{max}$} & \colhead{linelist}"+"\\\\\n")
f.writelines(r"\colhead{(GHz)} & \colhead{} & \colhead{(s$^{-1}$)} & \colhead{(K)} & \colhead{} & \colhead{} & \colhead{}"+"\\\\\n")
f.writelines(r"\colhead{(1)} & \colhead{(2)} & \colhead{(3)} & \colhead{(4)} & \colhead{(5)} & \colhead{(6)} & \colhead{(7)}"+"\n")
f.writelines(r"}"+"\n")
f.writelines(r"\startdata"+"\n")
f.writelines(r"spw 0a"+"\\\\\n")
for line in lines[3:5]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 0b"+"\\\\\n")
for line in lines[8:17]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 1a"+"\\\\\n")
for line in lines[20:21]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 1b"+"\\\\\n")
for line in lines[24:30]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 2a"+"\\\\\n")
for line in lines[33:39]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 2b"+"\\\\\n")
for line in lines[42:58]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 3a"+"\\\\\n")
for line in lines[61:68]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")
f.writelines(r"spw 3b"+"\\\\\n")
for line in lines[71:74]:
	line = line[:-1].split(" ")
	f.writelines(str(line[2])+" & "+str(line[3])+" & "+"{:.2E}".format(float(line[4]))+ " & "+str(round(float(line[5]),2))+" & "+str(line[6])+ " & "+str(line[7])+" & "+str(line[8])+"\\\\\n")

f.writelines(r"\enddata"+"\n")
f.writelines(r"\label{table:methylformate}"+"\n")
f.writelines(r"\end{deluxetable*}"+"\n")
f.close()
'''




