# search.py -- Functions for searching lines
# 
# Copyright (c) 2009-2010 Sebastien Maret and Pierre Hily-Blant
# 
# This file is part of Weeds.

import urllib
import urllib2
#./home/sw-astro/python/lib/python2.4/urllib2.py
import urlparse
import math
import votable
import cache
import line

timeout = 60. # Connection timeout for online searches, in seconds

# A few constants, in MKS units

speed_of_light = 299792458.0
boltzmann_constant = 1.38065040e-23
planck_constant = 6.62606896e-34

cm_K = planck_constant * speed_of_light / boltzmann_constant * 1e2 # cm-1 -> K

def search_slap(self, fmin, fmax, specie = "All", Eu_max = None):
    """
    Search lines in a database using SLAP

    This function make a query in a VO-compliant database using the
    Simple Line Access Protocol (SLAP).

    Arguments:
    fmin   -- the minimum frequency in MHz
    fmax   -- the maximum frequency in MHz
    specie -- the specie name (default All)
    Eu_max -- maximum upper level energy expressed
              in cm-1 (default None)
    
    """

    lines = []

    # Convert the frequencies into wavelenghts

    wlmin = speed_of_light / (fmax * 1e6) # m
    wlmax = speed_of_light / (fmin * 1e6)

    # Construct the URL for the database query
    # FixMe: Not all database support the CHEMICAL_ELEMENT
    # keyword. Filter this by hand?
    
    slap_url = self.url + "WAVELENGTH=%.9e/%.9e" % (wlmin, wlmax)
    if specie != "All":
        slap_url = slap_url + "&CHEMICAL_ELEMENT=%s" % specie

    slap_url = slap_url + "&REQUEST=queryData"
        
    # Get the results in the VO table format and parse it

    try:
        response = urllib2.urlopen(slap_url, timeout = timeout)
    except Exception, error:
        raise Exception, "Could not connect to database: %s" % error

    try:
        vot = votable.VOTable()
        vot.parse(response)

        # Find-out what fields are in the VO table. Each field is
        # identified with its utype or its name. Note that the utype
        # is defined by the SLAP standard, but the name may vary from
        # one database to the other.

        indexes = {}
        fields = vot.getFieldsAttrs()
        for i in range(len(fields)):
            if fields[i].has_key('utype'):
                field_id = fields[i]['utype']
            elif fields[i].has_key('name'):
                field_id = fields[i]['name']
            else:
                continue

            # Patch for voparis wrong utypes
            if field_id == "ldm:Line.wavelength":
                field_id = "ssldm:Line.wavelength.value"
            if field_id == "ldm:Line.einsteinA" and \
                    fields[i]['name'] == "log10_A_Einstein_coefficient":
                continue # ignore Einstein A in log units
            if field_id.split(":", 1)[0] == "ldm":
                field_id = "ssldm:" + field_id.split(":", 1)[1]

            indexes[field_id] = i

        # Read each row of the VO table an fill the line object.

        for r in vot.getDataRows():

            l = line.line()

            # Mandatory fields
            if indexes.has_key("ssldm:Line.species.name"):
                l.specie = vot.getData(r)[indexes["ssldm:Line.species.name"]]
            elif indexes.has_key("ssldm:Line.title"):
                l.specie = vot.getData(r)[indexes["ssldm:Line.title"]]
            else:
                raise Exception, "missing species name"
            if indexes.has_key("ssldm:Line.wavelength.value"):
                wl = float(vot.getData(r)[indexes["ssldm:Line.wavelength.value"]])
                l.frequency = speed_of_light / wl * 1e-6 # MHz
            else:
                raise Exception, "missing wavelenght"

            # Optional fields
            if indexes.has_key("ssldm:Line.einsteinA"):
                l.einstein_coefficient = float(vot.getData(r)[indexes["ssldm:Line.einsteinA"]])
            if indexes.has_key("ssldm:Line.initialLevel.energy"):
                l.upper_level.energy = float(vot.getData(r)[indexes["ssldm:Line.initialLevel.energy"]])
            if indexes.has_key("ssldm:Line.initialLevel.name"):
                l.upper_level.quantum_numbers = vot.getData(r)[indexes["ssldm:Line.initialLevel.name"]]
            if indexes.has_key("ssldm:Line.finalLevel.name"):
                l.lower_level.quantum_numbers = vot.getData(r)[indexes["ssldm:Line.finalLevel.name"]]

            # The following fields are not defined by the SLAP
            # standard, but may be returned by some databases (and are
            # needed for LTE modeling)
            if indexes.has_key("partition_function_link"):
                l.partition_function_link = vot.getData(r)[indexes["partition_function_link"]]
            if indexes.has_key("initial_statistical_weight"):
                l.upper_level.statistical_weight = float(vot.getData(r)[indexes["initial_statistical_weight"]])
            if indexes.has_key("final_statistical_weight_index"):
                l.lower_level.statistical_weight = float(vot.getData(r)[indexes["final_statistical_weight_index"]])

            # Skip lines whose upper level energy is lower than Eu
            if Eu_max:
                if l.upper_level.energy > Eu_max:
                    continue

            lines.append(l)

    except Exception, error:
        raise Exception, "Can't parse the response from the database: %s" % error

    return lines

def search_jpl_post(self, fmin, fmax, specie = "All", Eu_max = None):
    """
    Search lines in a the JPL database using HTTP/POST method

    Arguments:
    fmin   -- the minimum frequency in MHz
    fmax   -- the maximum frequency in MHz
    specie -- the specie name (default All)
    Eu_max -- maximum upper level energy expressed
              in cm-1 (default None)
    
    """

    def __decode_jpl_stat_weight(stat_weight):
        """
        Decode the statistical weight value from the JPL

        For some lines, the statistical weight is coded with letter,
        which allow to write number that would not otherwise fit on
        three digit. The letter correspond to the hundreds, i.e. A is
        100, B is 200, etc. Decode the value, and returns an integer.

        Arguments:
        stat_weight -- string coding the statistical weight

        """

        ascii_code = ord(stat_weight[0])
        if  65 <= ascii_code <= 90:
            hundred_digit = "%s" % (ascii_code - 55)
            stat_weight = hundred_digit + stat_weight[1:]

        return int(stat_weight)

    lines = []

    # In the JPL database, species are identified by a tag, so we must
    # first first find what is the tag of each specie. We can get
    # these from the catdir.cat file. This file also contains the
    # partition function values at 300 K which we need to compute the
    # Einstein coefficient of each line.

    try:

        specie_tag = self.data.specie_tag
        tag_specie = self.data.tag_specie
        specie_Q = self.data.specie_Q

    except AttributeError:
        
        specie_tag, tag_specie = {}, {}
        specie_Q = {}
        catdir_url = self.url + "/ftp/pub/catalog/catdir.cat"

        response = urllib2.urlopen(catdir_url)
        for l in response.readlines():
            try:
                t = int(l[0:6])
                s = l[7:20].rstrip()
                q = 10**float(l[26:33]) # Partition function at 300 K
            except (IndexError, TypeError), details:
                raise Exception, "Error while reading catalog file: %s" % details
            specie_tag[s] = t
            tag_specie[t] = s
            specie_Q[t] = q
        response.close()

        # Store the results to speed-up future searches

        self.data.specie_tag = specie_tag 
        self.data.tag_specie = tag_specie
        self.data.specie_Q = specie_Q
    
    # Make a HTTP/POST query on the server. Note that for the JPL
    # database, the keyword order does matter, so we must use a tuple
    # for the form values.

    fmin_s = "%.9f" % (fmin * 1e-3) # MHz -> GHz
    fmax_s = "%.9f" % (fmax * 1e-3) # MHz -> GHz

    if specie != "All":
        if specie in specie_tag.keys():
            tag = specie_tag[specie]
        else:
            return [] # specie is not in the database
    else:
        tag = "All"
        
    form_values = (('MinNu', fmin_s), ('MaxNu', fmax_s), ('UnitNu', "GHz"),
                   ('Mol', tag), ('StrLim', -500))

    try:

        search_url = self.url + "/cgi-bin/catform"
        data = urllib.urlencode(form_values)
        req = urllib2.Request(search_url, data)
        response = urllib2.urlopen(req, timeout = timeout)

    except Exception, error:
        raise Exception, "Could not connect to database: %s" % error

    # Parse the results

    for l in response.readlines():
        try:

            # Skip lines containing the specie name and tag
            if len(l) == 35:
                tg = int(l.split()[0])
                continue
            
            # If the search is incomplete, split the search frequency
            # range in two.
            # FixMe: maybe we should make sure that the lines at the
            # edge of the two frequency range don't appear twice in
            # the list because of rounding.
            # FixMe: this could be implemented at an higher level --
            # this might be useful for other databases too.
            if l[0:9] == "THIS form":
                lines1 = search_jpl_post(self, fmin, (fmin + fmax) / 2., 
                                         specie = specie, Eu_max = Eu_max)
                lines2 = search_jpl_post(self, (fmin + fmax) / 2., fmax,
                                         specie = specie, Eu_max = Eu_max)
                return lines1 + lines2
            
            freq = float(l[0:13]) # MHz
            wavelength = speed_of_light / (freq * 1e6) * 1e2 # cm
            log_intensity = float(l[21:29])
            lower_level_energy = float(l[31:41]) # cm-1
            upper_level_statistical_weight = __decode_jpl_stat_weight(l[41:44])
            upper_level_quantum_numbers = l[55:67].rstrip().lstrip()
            lower_level_quantum_numbers = l[67:].rstrip().lstrip()

            # Get the specie name and partition function
            spec = tag_specie[tg]
            Q = specie_Q[tg]

            # Compute the Einstein coefficient from the line intensity
            # at 300 K. For the formulae used see Eq. 5 on:
            # http://www.astro.uni-koeln.de/cdms/catalog
            upper_level_energy = lower_level_energy + 1. / wavelength # cm-1
            upper_level_energy_K = upper_level_energy * cm_K # K
            lower_level_energy_K = lower_level_energy * cm_K # K
            einstein_coefficient = 2.7964e-16 * 10**log_intensity * freq**2 \
                * Q / upper_level_statistical_weight \
                / (math.exp(-lower_level_energy_K / 300.) - math.exp(-upper_level_energy_K / 300.))

            sl =  line.line()
            sl.specie = spec
            sl.frequency = freq
            sl.einstein_coefficient = einstein_coefficient
            sl.upper_level.energy = upper_level_energy
            sl.upper_level.statistical_weight = upper_level_statistical_weight
            sl.upper_level.quantum_numbers = upper_level_quantum_numbers
            sl.lower_level.energy = lower_level_energy
            sl.lower_level.statistical_weight = upper_level_statistical_weight
            sl.lower_level.quantum_numbers = lower_level_quantum_numbers
        
            # Skip lines whose upper level energy is lower than Eu
            if Eu_max:
                if sl.upper_level.energy > Eu_max:
                    continue

            lines.append(sl)

        except Exception, error:
            raise Exception, "Can't parse the response from the database: %s" % error

    return lines

def search_cdms_post(self, fmin, fmax, specie = "All", Eu_max = None):
    """
    Search lines in a the CDMS database using HTTP/POST method

    Arguments:
    fmin   -- the minimum frequency in MHz
    fmax   -- the maximum frequency in MHz
    specie -- the specie name (default All)
    Eu_max -- maximum upper level energy expressed
              in cm-1 (default None)
    
    """

    lines = []


    # Make a HTTP/POST query on the server. Results of the query are
    # stored in a cache file, so the request is a two step process:
    # first we get the URL of the cache file, and then we read this
    # file.
    
    fmin = "%.9f" % (fmin * 1e-3) # MHz -> GHz
    fmax = "%.9f" % (fmax * 1e-3) # MHz -> GHz

    # FixMe: It looks that there is a bug in the database. If one
    # select a large (a few GHz) frequency range and request the
    # einstein coefficients (temp=0), then one get an internal server
    # error (HTTP 500). Let's request for the S*\mu^2 values instead
    # (temp=1) until they fix that.

    form_values = {'MinNu': fmin, 'MaxNu': fmax , 'UnitNu': "GHz",
                   'Moleculesgrp': "all species", 'StrLim': "-10",
                   "temp": "0", "output": "text", "sort": "frequency",
                   "mol_sort_query": "tag", "logscale": "yes",
                   "autoscale": "yes"}

    try:

        data = urllib.urlencode(form_values)
        req = urllib2.Request(self.url, data)
        response = urllib2.urlopen(req, timeout = timeout)
        
        base_url = urlparse.urlsplit(self.url).scheme + "://" \
            + urlparse.urlsplit(self.url).netloc
        cache_url = base_url \
            + response.read().split("\n")[4].split('"')[1]
        
        req = urllib2.Request(cache_url)
        response = urllib2.urlopen(req, timeout = timeout)

    except Exception, error:
        raise Exception, "Could not connect to database: %s" % error

    # Parse the results.
    for l in response.readlines()[10:-1]:
        try:

            freq = float(l[0:13]) # MHz
            wavelength = speed_of_light / (freq * 1e6) * 1e2 # cm
            einstein_coefficient = 10**float(l[24:35])
            lower_level_energy = float(l[37:47]) # cm-1
            upper_level_statistical_weight = int(l[47:50])
            upper_level_quantum_numbers = l[61:73].split(None, 0)[0].rsplit(None, 0)[0]
            lower_level_quantum_numbers = l[73:88].split(None, 0)[0].rsplit(None, 0)[0]
            spec = l[88:].split(None, 0)[0].rsplit(None, 0)[0]

            # Drop the asterisk at the beginning of some species names
            if spec[0] == "*":
                spec = spec[1:]

            sl =  line.line()
            sl.specie = spec
            sl.frequency = freq
            sl.einstein_coefficient = einstein_coefficient
            sl.upper_level.energy = lower_level_energy + 1. / wavelength # cm-1
            sl.upper_level.statistical_weight = upper_level_statistical_weight
            sl.upper_level.quantum_numbers = upper_level_quantum_numbers
            sl.lower_level.energy = lower_level_energy
            sl.lower_level.statistical_weight = upper_level_statistical_weight
            sl.lower_level.quantum_numbers = lower_level_quantum_numbers
        
            # Skip lines whose upper level energy is lower than Eu
            if Eu_max:
                if sl.upper_level.energy > Eu_max:
                    continue

            # Skip lines that does not correspond to the requested specie
            # FixMe: Widcards are ignored for the moment
            if specie != "All":
                if sl.specie != specie:
                    continue

            lines.append(sl)

        except Exception, error:

            # FixMe: Some species have missing entries. Ignore them
            # for the moment.

            continue
#        raise Exception, "Can't parse the response from the database: %s" % error

    return lines
         
def search_cache(self, fmin, fmax, specie = "All", Eu_max = None):
    """
    Search lines in the cache

    Arguments:
    fmin   -- the minimum frequency in MHz
    fmax   -- the maximum frequency in MHz
    specie -- the specie name (default All)
    Eu_max -- maximum upper level energy expressed
              in cm-1 (default None)
    
    """
    c = cache.cache(self.cache_file)
    return c.search(fmin,fmax,specie.replace('*','%'),Eu_max)

