""" Module with some simple stats classes
"""

import math

class WeightedStats(object):
    """ Perform simple weighted statistics
    """
    def __init__(self, data):
        self.data = data
        self.mean = None
        self.std_dev = None
        self.std_dev2 = None
        self.std_err_mean = None
        self.std_err_std_dev = None
        self.sumw = 0
        self.sumwerr2 = 0
        self.sumwx = 0
        self.sumwx2 = 0
        self.effective_entries = None
        self.calculate_stats()

    def calculate_stats(self):
        """ Calculates the relevant statistical variables
        """
        self.sumw = self.data["w"].sum()
        if self.sumw <= 0:
            return
        self.sumwx = (self.data["w"] * self.data["x"]).sum()
        self.sumwx2 = (self.data["w"] * self.data["x"] ** 2).sum()
        self.mean = self.sumwx / self.sumw
        
        self.std_dev2 = self.sumwx2 / self.sumw - self.mean * self.mean
        #in exceptional cases it may be so small that it gets negative because of rounding
        if self.std_dev2 > 0:
            self.std_dev = math.sqrt(self.std_dev2)
        else:
            self.std_dev = 0

        self.sumwerr2 = (self.data["werr"] ** 2).sum()
        if self.sumwerr2 <= 0:
            return
        self.effective_entries = self.sumw * self.sumw / self.sumwerr2
        sqrt_effective_entries = math.sqrt(self.effective_entries)
        self.std_err_mean = self.std_dev / sqrt_effective_entries
        self.std_err_std_dev = self.std_err_mean / math.sqrt(2)

        print(self.mean, self.std_err_mean)
        print(self.sumw, self.sumwerr2, self.sumwx, self.sumwx2)
