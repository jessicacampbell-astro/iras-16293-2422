# IRAS 16293-2422

This repository contains the files for the Leiden-ESA Astrophysics Summer Program (LEAPS) project on IRAS 16293-2422.

## Steps

1. Queried JPL and CDMS database for Acetaldehyde, Methanol, and Methyl Formate pure rotational transitions in the spw passbands (`DetectionFiles/SplatalogueQuery.py`)
	- spw 0a/b was temporarily ignored due to an apparent velocity offset of -2 km/s in the reduced data -- found to not have an offset by Mihkel via LTE model
    - used a lower limit of log10(Aij)=-5 for the Einstein coefficients; checked this by confirming detections with the lowest Einstein coefficient in each molecule (since this value might be dependent on the molecule)
    - used an upper limit of 1E3 K for the upper energy level
    - removed duplicates; CMDS > SPLAIN source catalogs (linelist); was still getting duplicates, but this turned out to be the result of v1.0 and v2.0 both being searched; Splatalogue source code was edited (by Magnus) and this fixed duplicate issues (yay)
    - there remained one duplicate for Methyl Formate due to the slight passband overlap between spw 1b and 2a -- assuming local-standard-of-rest source velocity of 2.7 km/s, line should be found in both data files -- there was no obvious detection in either, so the line was removed from its second occurrence (i.e., in spw 2a)
    - going back to spw 0a/b files, used a source velocity of 1 km/s instead of 3 km/s; used strong Acetaldehyde emission lines to find the shift
2. Created moment 0 (integrated line intensity), 1 (velocity), and 2 (linewidth) maps of the "B" source (`Maps/Maps.py`)
3. Created histograms of moment 0, 1, and 2 maps for each molecule in the "B" source
4. Found the brightest pixel in the moment 0 map of each emission line for all molecular species to visualize each emission line as a sanity check.
5. Created noise maps of each emission line for all molecular species -- tried various baselines to reduce source detection.
6. Partition function was fit for each molecule, using data from the CDMS and JPL catalogs via Splatalogue (`RotationalDiagrams/FitPartitionFunction.py`)
7. Calculated column densities of each rotational transition at each pixel where there was only emission and no absorption (`RotationalDiagrams/RotationalDiagrams.py`)
	- population diagram analysis only used where emission is assumed to be optically thin
	- data was converted from Jy/beam --> Jy
	- data was then converted from Jy --> K
8. Created histograms of Nthin column densities of each molecule in the "B" source region (`RotationalDiagrams/RotationalDiagrams.py`)
9. Created rotational temperature maps using population diagrams for each molecule using partition function and column densities of each transition (`RotationalDiagrams/RotationalDiagrams.py`)
	- compared rotational temperature maps using all transitions of methanol versus methanol transitions with the same K (K=3) quantum number
