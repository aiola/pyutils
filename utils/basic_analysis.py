""" This module analyses the ntuples obtained from the Sbt software.
Usage: ipython -i basic_analysis.py config.yaml
"""

import os
import math
import exceptions
import uproot
import ROOT

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
        self.data_frame = None
        self.data_frames = dict()
        self.data_frame_definitions = config["data_frames"]
        self.flatten = config["flatten"]

    def open_tree(self, file_name):
        """ Open the Sbt ntuple
        """
        print("Opening file '{}'...".format(file_name))
        self.current_file = uproot.open(file_name)
        if not self.current_file:
            raise exceptions.RuntimeError("Could not open file '{}'".format(file_name))
        self.current_tree = self.current_file[self.tree_name]
        if not self.current_tree:
            raise exceptions.RuntimeError("Could not find tree '{}' in file '{}'".format(self.tree_name, file_name))
        for df_def in self.data_frame_definitions:
            df = self.current_tree.pandas.df(df_def["branches"], flatten=self.flatten)
            self.data_frames[df_def["name"]] = df
        if "default" in self.data_frames:
            self.data_frame = self.data_frames["default"]

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
