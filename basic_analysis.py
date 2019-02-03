""" This module analyses the ntuples obtained from the Sbt software.
Usage: ipython -i basic_simulation_analysis.py config.yaml
"""

import os
import exceptions
import ROOT

class BasicAnalysis(object):
    """ Main class
    """
    def __init__(self, config):
        self.path = config["path"]
        self.file_names = config["file_names"]
        self.tree_name = config["tree_name"]
        self.min_n_det = config["min_n_det"]
        self.current_file_name = None
        self.current_file_title = None
        self.current_file = None
        self.current_tree = None
        self.histograms = dict()
        self.canvases = []

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

    def basic_plot(self):
        """ Plots the results
        """
        for title, histograms in self.histograms:
            for hist in self.histograms.itervalues():
                cname = "{}_{}".format(title, hist.GetName())
                canvas = ROOT.TCanvas(cname, cname)
                self.canvases.append(canvas)
                canvas.cd()
                hist.Draw()
                canvas.Update()
    
    def save_plots(self):
        """ Save the plots on disk
        """
        for canvas in self.canvases:
            canvas.SaveAs("{}/{}.pdf".format(self.path, canvas.GetName()))
