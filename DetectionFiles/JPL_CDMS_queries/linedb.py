# linedb.py -- Line database classes and instances
# 
# Copyright (c) 2009-2010 Sebastien Maret and Pierre Hily-Blant
# 
# This file is part of Weeds.

#import numpy
import os
from search import *
from partition_function import *
import cache
import progressbar as pb

class linedb:
    """
    Molecular and atomic line database

    This class defines a molecular and line database that can be
    searched online or offline. Database can be searched online
    through the SLAP protocol, or directly with HTTP POST or GET
    requests. Databases can also be fetched entirely and cached on a
    disk as an a SQLite database for offline searches.

    """

    def __init__(self, url, cache_file, protocol, Q_format, online = True):
        """
        Create a database instance

        This function create a database instance which is accessed
        through a given protocol. Supported protocols are "slap",
        "cdms_post" and "jpl_post".

        Arguments:
        url         -- The URL of the database
        protocol    -- The database protocol
        cache_file  -- The name of cache file
        Q_format    -- The partition function file format
        online      -- Search the online database (default True)
        
        """

        if not protocol in ["slap", "cdms_post", "jpl_post"]:
            raise ValueError, "Unknown protocol"
        if not Q_format in ["votable", "cdms", "jpl"]:
            raise ValueError, "Unknown partition function file format"

#        class linedb_data_class():
        class linedb_data_class:
            """
            Storage class for linedb run time data
            
            """
#            def __init__(self, contents=None): 
#                self.contents = contents or []

        self.url = url
        self.protocol = protocol
        self.cache_file = os.path.expanduser(cache_file)
        self.Q_format = Q_format
        self.online = online
        self.data = linedb_data_class()

    def setonline(self, online):
        """
        """
        if not online:
            if os.path.exists(self.cache_file):
                self.online = online
            else:
                raise ValueError, "Offline cache not available"
        else:
            # Should try a request to see whether it is available
            self.online = online

    def search(self, fmin, fmax, specie = 'All', Eu_max = None):
        """
        Search for lines in the database
        
        Arguments:
        fmin   -- the minimum frequency in MHz
        fmax   -- the maximum frequency in MHz
        specie -- the specie name (default All)
        Eu_max -- maximum upper level energy expressed
                  in cm-1 (default None)
            
        """

        # Select the search method

        if self.online:
            if self.protocol == "slap":
                s = search_slap
            elif self.protocol == "jpl_post":
                s = search_jpl_post
            elif self.protocol == "cdms_post":
                s = search_cdms_post
            else:
                raise ValueError, "Unknown protocol"
        else:
            s = search_cache

        # Make the search
        
        lines = s(self, fmin, fmax, specie, Eu_max)
            
        return lines
        
    def cache(self, fmin = 0, fmax = 2e6, fstep = 1e4):
        """
        Cache an online database
        
        Arguments:
        fmin   -- the minimum frequency in MHz
        fmax   -- the maximum frequency in MHz
        fstep  -- the frequency step in MHz between each query
        
        """

        c = cache.cache(self.cache_file)
        c.create( "0.10", fmin, fmax, True)

        # Add the lines to the database

        species = []
        f = fmin
        if (fstep > fmax - fmin):
            fstep = fmax - fmin

        widgets = ['Downloading lines:      ', pb.Percentage(), ' ', pb.Bar(), ' ',
                   pb.ETA(), ' ']
        pbar = pb.ProgressBar(widgets=widgets, maxval= fmax - fmin).start()

        while f < fmax:
            self.online = True # Make sure to search online!
            lines = self.search(fmin = f, fmax = f + fstep + 1e-4)
            for l in lines:
                if not(l.specie) in species:
                    species.append(l.specie)
            c.add_lines(lines)
            f += fstep
            pbar.update(f - fmin)
            
        pbar.finish()

        # Add the partition functions
        # FixMe: implement a function to get the partition function
        # values as a function of the temperature.

        if self.Q_format == "cdms":
            pf = partition_function_cdms
        elif self.Q_format == "jpl":
            pf = partition_function_jpl
        elif self.Q_format == "votable":
            raise NotImplementedError, "Caching of partition function " \
                "in VO table format is not implemented."
        else:
            raise ValueError, "Unknown partition function format"

        widgets = ['Downloading part.func.: ', pb.Percentage(), ' ',
                   pb.Bar(), ' ', pb.ETA(), ' ']
        pbar = pb.ProgressBar(widgets=widgets, maxval= len(species)).start()

        for i in range(len(species)):
            try:
                temperature, partfunc = pf(self, species[i])
                c.add_partfunc(species[i], temperature, partfunc)
            except NotFoundError:
                continue
            pbar.update(i)
            
        pbar.finish()

    def update(self):
        """
        Update the cache file of a database
        
        """

        pass

    def partition_function(self, T, *arg, **kw):
        """
        Returns the partition function
    
        This function returns the partition function of a specie at a
        given temperature. The partition function values at different
        temperatures is fetched from the database (online searches) or
        queried in the cache (offline searches). The value of the
        partition function at the requested temperature is obtained by
        a linear fit of the partition function values at the the two
        closest temperature from the desired value.
        
        Arguments:
        url    -- the URL of partition function file, for votable format
        specie -- the name of the specie, for other formats and offline searches
        T      -- the temperature
        
        """

        if self.online and self.Q_format == "votable":
            if not "url" in kw.keys():
                raise ValueError, "The URL of the partition function must be provided."
        else:
            if not "specie" in kw.keys():
                raise ValueError, "The specie name must be provided."

        if not(self.online):
            temperature, partfunc = partition_function_cache(self, kw["specie"])
        elif self.Q_format == "votable":
            temperature, partfunc = partition_function_votable(self, kw["url"])
        elif self.Q_format == "cdms":
            temperature, partfunc = partition_function_cdms(self, kw["specie"])
        elif self.Q_format == "jpl":
            temperature, partfunc = partition_function_jpl(self, kw["specie"])
        else:
            raise ValueError, "Unknown partition function format"

        # Make a linear fit of the two closest values
        #temperature = numpy.array(temperature)
        #partfunc = numpy.array(partfunc)
        index = abs((temperature - T)).argsort()
        t0 = temperature[index[0]]
        t1 = temperature[index[1]]
        p0 = partfunc[index[0]]
        p1 = partfunc[index[1]]
        result = (p1 - p0) / (t1 - t0) * T + p0 - (p1 - p0) / (t1 - t0) * t0

        return result

voparis = linedb(url = "http://linelists.obspm.fr/transitions.php?base=cdms&",
                 cache_file = "/Users/mihkelkama/Codes/HWeeds/voparis.db", protocol = "slap",
                 Q_format = "votable", online = True)

vamdc = linedb(url = "http://casx019-zone1.ast.cam.ac.uk/vamdc-slap/slap?",
               cache_file = "/Users/mihkelkama/Codes/HWeeds/vamdc.db", protocol = "slap",
               Q_format = "votable", online = True)

splatalogue = linedb(url = "http://find.nrao.edu/splata-slap/slap?",
                      cache_file = "/Users/mihkelkama/Codes/HWeeds/splatalogue.db", protocol = "slap",
                      Q_format = "votable", online = True)

cdms = linedb(url = "http://www.astro.uni-koeln.de/cgi-bin/cdmssearch",
              cache_file = "/Users/mihkelkama/Codes/HWeeds/weedsCache_cdms.db", protocol = "cdms_post",
              Q_format = "cdms", online = True)

jpl = linedb(url = "http://spec.jpl.nasa.gov",
             cache_file = "/Users/mihkelkama/Codes/HWeeds/weedsCache_jpl.db", protocol = "jpl_post",
             Q_format = "jpl", online = True)

default = cdms
