# partition_function.py -- Functions to get partition functions
# 
# Copyright (c) 2009-2010 Sebastien Maret and Pierre Hily-Blant
# 
# This file is part of Weeds.

import urllib2
#./home/sw-astro/python/lib/python2.4/urllib2.py
import votable
import cache

class NotFoundError(Exception):
    pass

def partition_function_votable(self, url):
    """
    Returns the partition function at different temperatures
    
    This function fetches a VO table containing the partition function
    at different temeperatures from the given URL, and then parse it.
    
    Arguments:
    url -- the URL of the VO table
    
    """
    
    vot = votable.VOTable()
    vot.parse (urllib2.urlopen(url))
    
    temperature_index = vot.getColumnIdx('temperature')
    value_index = vot.getColumnIdx('value')
    
    temperature = []
    partfunc = []
    
    for r in vot.getDataRows():
        temperature.append(float(vot.getData(r)[temperature_index]))
        partfunc.append(float(vot.getData(r)[value_index]))

    return temperature, partfunc

def partition_function_jpl(self, specie):
    """
    Returns the partition function at different temperatures

    This function fetches a text file (in the format used in the JPL
    database) containing the partition function of the given specie
    for different temperatures and then parse it. The file content is
    kept in memory to avoid re-fetching if the function is called
    several times.
    
    Arguments:
    specie -- the specie name
    
    """

    global weeds_jpl_partition_function_file

    catdir_url = self.url + "/ftp/pub/catalog/catdir.cat"
    temperature = [300., 225., 150., 75., 37.5, 18.75, 9.375]

    try:
        weeds_jpl_partition_function_file
    except NameError:
        f = urllib2.urlopen(catdir_url)
        weeds_jpl_partition_function_file = f.readlines()
        f.close()

    for l in weeds_jpl_partition_function_file:
        s = l[7:20].rstrip()
        if s == specie:
            q1 = 10**float(l[26:33])
            q2 = 10**float(l[33:40])
            q3 = 10**float(l[40:47])
            q4 = 10**float(l[47:54])
            q5 = 10**float(l[54:61])
            q6 = 10**float(l[61:68])
            q7 = 10**float(l[68:75])
            partition_function = [q1, q2, q3, q4, q5, q6, q7]            
            return temperature, partition_function
        
    raise NotFoundError, "No partition function found for %s." % specie

def partition_function_cdms(self, specie):
    """
    Returns the partition function at different temperatures

    This function fetches a text file (in the format used in the CDMS
    database) containing the partition function of the given specie
    for different temperatures and then parse it. The file content is
    kept in memory to avoid re-fetching if the function is called
    several times.
    
    Arguments:
    specie -- the specie name
    
    """

    global weeds_cdms_partition_function_file

    # FixMe: give relative URL
    partfunc_url = "http://www.astro.uni-koeln.de/site/vorhersagen/catalog/partition_function.html"
    temp = [1000., 500., 300., 225., 150., 75., 37.5, 18.75, 9.375]

    try:
        weeds_cdms_partition_function_file
    except NameError:
        f =  urllib2.urlopen(partfunc_url)
        weeds_cdms_partition_function_file = f.readlines()
        f.close()

    temperature = []
    partition_function = []

    for l in weeds_cdms_partition_function_file:
        if l[0] == "<":
            continue
        try:
            spec = l[7:28].rstrip()
            if spec == "":
                continue
            if spec == specie:
                field = l[40:].split()
                for i in range(len(temp)):
                    if field[i] == "---":
                        continue
                    temperature.append(temp[i])
                    partition_function.append(10**float(field[i]))
        except:
            continue

    if partition_function == []:
        raise NotFoundError, "No partition function found for %s." % specie

    return temperature, partition_function

def partition_function_cache(self, specie):
    """
    Returns the partition function at different temperatures

    This function fetches the partition function values at different
    temperatures from the cache for the given specie.
    
    Arguments:
    specie -- the specie name
    
    """

    c = cache.cache(self.cache_file)

    return c.partition_function(specie)

