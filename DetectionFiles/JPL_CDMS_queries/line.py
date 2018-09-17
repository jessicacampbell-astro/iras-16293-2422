# line.py -- Classes and methods for spectral lines
# 
# Copyright (c) 2009 Sebastien Maret and Pierre Hily-Blant
# 
# This file is part of Weeds.

class line:
    """
    Spectral line

    """

    def __init__(self):
        self.specie = ''
        self.frequency = 0.
        self.einstein_coefficient = 0.
        self.upper_level = level()
        self.lower_level = level()

    def __repr__(self):
        return "%-16s | %12.6f | %6.1f | %s -- %s" % \
        (self.specie, self.frequency*1e-3, self.upper_level.energy * 1.44,
        self.upper_level.quantum_numbers, self.lower_level.quantum_numbers)

class level:
    """
    Energy level

    """

    def __init__(self):
        self.energy = 0.
        self.statistical_weight = 0.
        self.quantum_numbers = ''

if __name__ == "__main__":
    l = line()
    print l


