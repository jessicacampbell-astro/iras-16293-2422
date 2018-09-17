import numpy as np 

detections_A = 0
EU_A = []
Aij_A = []
g_A = []

detections_M = 0
EU_M = []
Aij_M = []
g_M = []

detections_MF = 0
EU_MF = []
Aij_MF = []
g_MF = []


# acetaldehyde
f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/DetectedAcetaldehyde.txt")
lines = f.readlines()
f.close()

# 0a
for line in lines[3:26]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 0b
for line in lines[29:34]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 1a
for line in lines[37:43]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 1b
for line in lines[46:69]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 2a
for line in lines[72:77]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 2b
for line in lines[80:95]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 3a
for line in lines[98:107]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1
# 3b
for line in lines[110:121]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_A.append(EU)
	Aij_A.append(Aij)
	g_A.append(g)
	detections_A +=1

print "Acetaldehyde"
print "Detections:", detections_A
print "E_U (K):", np.mean(EU_A), np.min(EU_A), np.max(EU_A)
print "Aij:", np.mean(Aij_A), np.min(Aij_A), np.max(Aij_A)
print "g:", np.mean(g_A), np.min(g_A), np.max(g_A)


# methanol
f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/DetectedMethanol.txt")
lines = f.readlines()
f.close()

# 0a
#for line in lines[3:26]:
#	line = line[:-1].split(" ")
#	Aij = float(line[4])
#	EU = float(line[5])
# 	g = float(line[6])
#	EU_M.append(EU)
#	Aij_M.append(Aij)
# 	g_M.append(g)
#	detections_M +=1
# 0b
for line in lines[6:11]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1
# 1a
for line in lines[14:16]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1
# 1b
for line in lines[19:20]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1
# 2a
for line in lines[23:26]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1
# 2b
for line in lines[29:31]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1
# 3a
for line in lines[34:35]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1
# 3b
for line in lines[38:42]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_M.append(EU)
	Aij_M.append(Aij)
	g_M.append(g)
	detections_M +=1

print "\nMethanol"
print "Detections:", detections_M
print "E_U (K):", np.mean(EU_M), np.min(EU_M), np.max(EU_M)
print "Aij:", np.mean(Aij_M), np.min(Aij_M), np.max(Aij_M)
print "g:", np.mean(g_M), np.min(g_M), np.max(g_M)


# Methyl Formate
f = open("/Users/jessica/Documents/LEAPS2015/DetectionFiles/JPL_CDMS_queries/DetectedMethylFormate.txt")
lines = f.readlines()
f.close()

# 0a
for line in lines[3:5]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 0b
for line in lines[8:17]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 1a
for line in lines[20:21]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 1b
for line in lines[24:30]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 2a
for line in lines[33:39]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 2b
for line in lines[42:58]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 3a
for line in lines[61:68]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1
# 3b
for line in lines[71:74]:
	line = line[:-1].split(" ")
	Aij = float(line[4])
	EU = float(line[5])
	g = float(line[6])
	EU_MF.append(EU)
	Aij_MF.append(Aij)
	g_MF.append(g)
	detections_MF +=1

print "\nMethyl Formate"
print "Detections:", detections_MF
print "E_U (K):", np.mean(EU_MF), np.min(EU_MF), np.max(EU_MF)
print "Aij:", np.mean(Aij_MF), np.min(Aij_MF), np.max(Aij_MF)
print "g:", np.mean(g_MF), np.min(g_MF), np.max(g_MF)




