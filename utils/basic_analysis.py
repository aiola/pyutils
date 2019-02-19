""" This module analyses the ntuples obtained from the Sbt software.
Usage: ipython -i basic_analysis.py config.yaml
"""

import os
import math
import exceptions
import ROOT

class MeasuredQuantity(object):
    def __init__(self, value, error, units=""):
        if error < 0:
            raise exceptions.ValueError()
        if math.isnan(value) or math.isnan(error):
            raise exceptions.ValueError()
        if units == "%":
            value *= 100
            error *= 100
        self.value = value
        self.error = error
        self.units = units

    def is_significant(self):
        if self.value == 0:
            return self.error == 0
        return self.error / math.fabs(self.value) < 1.0

    def calculate_precision(self):
        if self.error == 0:
            return int("inf")
        
        vErrLog10 = math.log10(self.error)
        if vErrLog10 < 0:
            return int(math.floor(vErrLog10))
        else:
            return int(math.ceil(vErrLog10))

    def to_string(self, pm = "#pm"):
        precision = self.calculate_precision()
        if math.isinf(precision):
            precision = 6
        abs_precision = abs(precision)
        if abs_precision > 3:
            if precision > 0:
                exp = precision - 1
            else:
                exp = precision + 1
            result = "({value:.1f} {pm} {error:.1f}) #times 10^{exp}".format(value = self.value / (10**exp), error = self.error / (10**exp), exp=exp, pm=pm)
        else:
            format_string = "{{value:.{prec}f}} {{pm}} {{error:.{prec}f}}".format(prec = abs_precision)
            result = format_string.format(value = self.value, error=self.error, pm=pm)
        if self.units:
            result += " {}".format(self.units)
        return result

    def __str__(self):
        return self.to_string("+/-")

    def __add__(self, other):
        if isinstance(other, MeasuredQuantity):
            if self.units != other.units:
                return NotImplemented
            return MeasuredQuantity(self.value + other.value, math.sqrt(self.error**2 + other.error**2), self.units)
        elif isinstance(other, (int, float)):
            if self.units:
                return NotImplemented
            return MeasuredQuantity(self.value + other, self.error)
        return NotImplemented

    def __neg__(self):
        return MeasuredQuantity(-self.value, self.error, self.units)

    def __sub__(self, other):
        return self + (-other)
    
    def __mul__(self, other):
        if isinstance(other, MeasuredQuantity):
            value = self.value * other.value
            error = math.sqrt(self.error**2 / self.value**2 + other.error**2 / other.value**2) * value
            units = self.units
            if units and other.units:
                if units == other.units:
                    if units == "%":
                        value /= 100
                        error /= 100
                    else:
                        units = "({})^2".format(units)
                else:
                    units += " " + other.units
            return MeasuredQuantity(value, error, units)
        elif isinstance(other, (int, float)):
            return MeasuredQuantity(self.value * other, self.error * other)
        else:
            return NotImplemented
    
    def __invert__(self):
        if self.units != "%":
            return MeasuredQuantity(1 / self.value, self.error / self.value**2, self.units)
        return NotImplemented

    def __div__(self, other):
        return self * (~other)

    def __lt__(self, other):
        if isinstance(other, MeasuredQuantity):
            return self.value + self.error < other.value - other.error
        elif isinstance(other, (int, float)):
            return self.value + self.error < other
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, MeasuredQuantity):
            return self.value + self.error > other.value - other.error
        elif isinstance(other, (int, float)):
            return self.value + self.error > other
        return NotImplemented

    def ___le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

    def __ne__(self, other):
        return self > other or self < other

    def __eq__(self, other):
        return not self != other

class BasicHistogramContainer(object):
    """ Histogram container
    """
    def __init__(self, name):
        self.name = name
        self.histograms = dict()
        self.canvases = []
        self.histograms_on_canvases = []

    def save_plots(self, path):
        """ Save the plots on disk
        """
        for canvas in self.canvases:
            canvas.SaveAs("{}/{}.pdf".format(path, canvas.GetName()))

    def basic_plot(self):
        """ Plots the results
        """
        for hist in self.histograms.itervalues():
            if not isinstance(hist, ROOT.TH1):
                continue
            cname = "{}".format(hist.GetName())
            canvas = ROOT.TCanvas(cname, cname)
            self.canvases.append(canvas)
            canvas.cd()
            hist.Draw()
            canvas.Update()

class BasicAnalysis(object):
    """ Base class for any analysis
    """
    def __init__(self, config, debug=False):
        self.path = config["path"]
        self.tree_name = config["tree_name"]
        self.current_tree = None
        self.current_file = None
        self.debug = debug
        self.histogram_containers = dict()
        self.canvases = []
        self.histograms_on_canvases = []
        self.output_path = config["output_path"]

    def open_tree(self, file_name):
        """ Open the Sbt ntuple
        """
        print("Opening file '{}'...".format(file_name))
        self.current_file = ROOT.TFile(file_name)
        if not self.current_file or self.current_file.IsZombie():
            raise exceptions.RuntimeError("Could not open file '{}'".format(file_name))
        self.current_tree = self.current_file.Get(self.tree_name)
        if not self.current_tree:
            raise exceptions.RuntimeError(
                "Could not find tree '{}' in file '{}'".format(self.tree_name, file_name))

    def build_histograms(self):
        """ Build the histograms where the ntuple will be projected
        """
        for histogram_container in self.histogram_containers.itervalues():
            if isinstance(histogram_container, BasicHistogramContainer):
                histogram_container.build_histograms()
            elif isinstance(histogram_container, list):
                for container_inside in histogram_container:
                    container_inside.build_histograms()

    def basic_plot(self):
        """ Plotting
        """
        for histogram_container in self.histogram_containers.itervalues():
            if isinstance(histogram_container, BasicHistogramContainer):
                histogram_container.basic_plot()
            elif isinstance(histogram_container, list):
                for container_inside in histogram_container:
                    container_inside.basic_plot()

    def save_plots(self):
        """ Save the plots on disk
        """
        output_path = "{}/../{}".format(self.path, self.output_path)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        for canvas in self.canvases:
            canvas.SaveAs("{}/{}.pdf".format(output_path, canvas.GetName()))
        for histogram_container in self.histogram_containers.itervalues():
            if isinstance(histogram_container, BasicHistogramContainer):
                histogram_container.save_plots(output_path)
            elif isinstance(histogram_container, list):
                for container_inside in histogram_container:
                    container_inside.save_plots(output_path)

class SingleFileAnalysis(BasicAnalysis):
    """ Base class for analysis based on a single ROOT input file
    """
    def __init__(self, config, debug=False):
        BasicAnalysis.__init__(self, config, debug)
        self.file_name = config["file_name"]

    def open_tree_from_file(self):
        """ Open the Sbt ntuple
        """
        full_file_name = "{}/{}.root".format(self.path, self.file_name)
        self.open_tree(full_file_name)

class MultiFileAnalysis(BasicAnalysis):
    """ Base class for analysis based on multiple ROOT input files
    """
    def __init__(self, config, debug=False):
        BasicAnalysis.__init__(self, config, debug)
        self.file_names = config["file_names"]
        self.current_file_name = None
        self.current_file_title = None

    def open_next_file(self):
        """ Opens the next file
        """
        for title, file_name in self.file_names.iteritems():
            full_file_name = "{}/{}.root".format(self.path, file_name)
            self.current_file_name = full_file_name
            self.current_file_title = title
            self.open_tree(full_file_name)
            yield title, file_name

def add_histogram_to_dict(histograms, histogram):
    """ Add an histogram to a dictionary using its own name as key
    """
    histograms[histogram.GetName()] = histogram
