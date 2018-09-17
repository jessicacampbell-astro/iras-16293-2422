# cache.py -- Classes and methods to cache a line database
# 
# Copyright (c) 2009 Emmanuel Reynier and Sebastien Maret
# 
# This file is part of Weeds.

import sqlite3
import array
import os
import line

class cache:
    """
    Cache management class for weeds database

    This class allows to create a local copy of a line database. This
    is done by fetching the entire (or part of) a database and storing
    it on the disk as a SQLite database.
    
    """

    def __init__(self, db_file):
        """
        Create a cache instance
        
        Arguments:
        db_file   -- name of the SQLite dabase file

        """

        self.db_file = db_file

    def create(self, version, fmin, fmax, overwrite = False):
        """
        Create a SQL database containing the line catalog and
        partition functions

        Arguments:
        version   -- version of database
        fmin      -- minimal frequency
        fmax      -- maximal frequency
        overwrite -- overwrite existing db file (default False)

        """

        if os.path.isfile(self.db_file):
            if overwrite:
                os.remove(self.db_file)
            else:
                raise ValueError, "Database file %s already exists." % self.db_file

        db_connect = sqlite3.connect(self.db_file)
        db_cursor = db_connect.cursor()
        db_cursor.execute("create table line ("
                          "'specie' char(32),"
                          "'frequency' real,"
                          "'einstein_coefficient' real,"
                          "'upper_level_energy' real,"
                          "'upper_level_statistical_weight' real,"
                          "'upper_level_quantum_numbers' char(32),"
                          "'lower_level_energy' real,"
                          "'lower_level_statistical_weight' real,"
                          "'lower_level_quantum_numbers' char(32)"
                          ");")
        db_cursor.execute("create table partfunc ("
                          "'specie' char(32),"
                          "'temperature' blob,"
                          "'partfunc' blob"
                          ");")
        db_cursor.execute("create index 'specie' on line('specie');")
        db_cursor.execute("create index 'frequency' on line('frequency');")
        db_cursor.execute("create table info ("
                          "'version' char(32),"
                          "'min_frequency' real,"
                          "'max_frequency' real"
                          ");")
        db_cursor.execute("insert into info "
                          "values('%s',%f,%f);" % (version, fmin, fmax)
                          )
        db_connect.commit()
        db_cursor.close()

    @staticmethod
    def __execute_insert_line(db_cursor, line):
 
        db_cursor.execute("insert into line "
                          "values (\"%s\",%f,%f,%f,%f,\"%s\",%f,%f,\"%s\");"
                          % (line.specie,
                             line.frequency,
                             line.einstein_coefficient,
                             line.upper_level.energy,
                             line.upper_level.statistical_weight,
                             line.upper_level.quantum_numbers,
                             line.lower_level.energy,
                             line.lower_level.statistical_weight,
                             line.lower_level.quantum_numbers)
                          )

    @staticmethod
    def __execute_insert_partfunc(db_cursor, specie, temperature, partfunc):

        t = sqlite3.Binary(array.array('f', temperature).tostring())
        p = sqlite3.Binary(array.array('f', partfunc).tostring())
        
        query = u'''insert into partfunc values(?,?,?)'''
        db_cursor.execute(query, (specie, t, p))

    def add_line(self, line):
        """
        Add a line to the database

        Arguments:
        line -- line object

        """

        db_connect = sqlite3.connect(self.db_file)
        db_cursor = db_connect.cursor()
        self.__execute_insert_line(db_cursor, line)
        db_connect.commit()
        db_cursor.close()

    def add_lines(self, lines):
        """
        Add a list of lines to the database

        Arguments:
        lines -- line list

        """

        db_connect = sqlite3.connect(self.db_file)
        db_cursor = db_connect.cursor()
        for line in lines:
            self.__execute_insert_line(db_cursor, line)
        db_connect.commit()
        db_cursor.close()

    def add_partfunc(self, specie, temperature, partfunc):
        """
        Add a partition function to the database

        Arguments:
        specie      -- the specie name
        temperature -- list of temperatures
        partfunc    -- list of partition function values at
                       these temperatures

        """

        db_connect = sqlite3.connect(self.db_file)
        db_cursor = db_connect.cursor()
        self.__execute_insert_partfunc(db_cursor, specie, temperature, 
                                       partfunc)
        db_connect.commit()
        db_cursor.close()

    def add_partfuncs(self, species, temperatures, partfuncs):
        """
        Add a partition function to the database

        Arguments:
        species      -- list of specie names
        temperatures -- list of list of temperatures
        partfuncs    -- list of list of partition function values at
                        these temperatures

        """

        db_connect = sqlite3.connect(self.db_file)
        db_cursor = db_connect.cursor()
        for specie, temperature, partfunc in zip(species, temperatures, partfuncs):
            self.__execute_insert_partfunc(db_cursor, specie, temperature, 
                                           partfunc)
        db_connect.commit()
        db_cursor.close()

    def info(self):
        """
        Display informations on the database
        
        """

        db_connect = sqlite3.connect(self.db_file)
        db_connect.row_factory = sqlite3.Row
        db_cursor = db_connect.cursor()
        db_cursor.execute( "select * from info")
        row = db_cursor.fetchone()
        print
        print '****** Database informations ******'
        print '* version:           %12s *' % row['version']
        print '* minimal frequency: %12s *' % row['min_frequency']
        print '* maximal frequency: %12s *' % row['max_frequency']
        print '***********************************'
        print
        db_connect.commit()
        db_cursor.close()

    def remove(self, freq, freq_precision = 1e-4, specie = "All"):

        self.remove_range(freq - freq_precision, freq + freq_precision, specie)

    def remove_range(self, freq_min, freq_max, specie = "All"):

        db_connect = sqlite3.connect(self.db_file)
        db_connect.row_factory = sqlite3.Row
        db_cursor = db_connect.cursor()

        lines = []
        query = "delete from line where frequency >= %f and frequency <= %f" % \
            (freq_min, freq_max)
        if specie != "All":
            query += " and specie like '%s'" % specie

        db_cursor.execute( query)

        db_connect.commit()
        db_cursor.close()

    def search(self, fmin, fmax, specie = "All", Eu_max = None):
        """
        Search lines in the cache

        Arguments:
        fmin   -- the minimum frequency in MHz
        fmax   -- the maximum frequency in MHz
        specie -- the specie name (default All)
        Eu_max -- maximum upper level energy expressed
                  in cm-1 (default None)        

        """

        db_connect = sqlite3.connect(self.db_file)
        db_connect.row_factory = sqlite3.Row
        db_cursor = db_connect.cursor()

        lines = []
        query = "select * from line where frequency >= %f and frequency <= %f" % (fmin, fmax)
        if specie != "All":
            query += " and specie like '%s'" % specie
        if Eu_max:
            query += " and upper_level_energy <= %f" % Eu_max

        db_cursor.execute(query)
        for row in db_cursor:
            l = line.line()
            l.specie = row ['specie']
            l.frequency = row ['frequency']
            l.einstein_coefficient = row ['einstein_coefficient']
            l.upper_level.energy = row ['upper_level_energy']
            l.upper_level.statistical_weight = row ['upper_level_statistical_weight']
            l.upper_level.quantum_numbers = row ['upper_level_quantum_numbers']
            l.lower_level.energy = row ['lower_level_energy']
            l.lower_level.statistical_weight = row ['lower_level_statistical_weight']
            l.lower_level.quantum_numbers = row ['lower_level_quantum_numbers']
            lines.append(l)

        db_connect.commit()
        db_cursor.close()

        return lines

    def partition_function(self, specie):
        """
        Returns the partition function at different temperatures
        
        This function get the partition function of the given specie
        for different temperatures from the cache.
        
        Arguments:
        specie -- the specie name
    
        """

        db_connect = sqlite3.connect(self.db_file)
        db_connect.row_factory = sqlite3.Row
        db_cursor = db_connect.cursor()        

        query = "select * from partfunc where specie = '%s';" % specie
        db_cursor.execute(query)

        temperature = array.array('f', [])
        partfunc = array.array('f', [])

        for row in db_cursor:
            temperature.fromstring(row['temperature'])
            partfunc.fromstring(row['partfunc'])
            break

        if len(partfunc) == 0:
            raise Exception, "No partition function found for this specie."

        return temperature, partfunc
