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
        self.file_name = config["file_name"]
        self.tree_name = config["tree_name"]
        self.min_n_det = config["min_n_det"]
        self.file = None
        self.tree = None
        self.histograms = dict()
        self.canvases = []

    def open_tree(self):
        """ Open the Sbt ntuple
        """
        self.file = ROOT.TFile(self.file_name)
        if not self.file or self.file.IsZombie():
            raise exceptions.RuntimeError("Could not open file '{}'".format(self.file_name))
        self.tree = self.file.Get(self.tree_name)
        if not self.tree:
            raise exceptions.RuntimeError(
                "Could not find tree '{}' in file '{}'".format(self.tree_name, self.file_name))

    def plot(self):
        """ Plots the results
        """
        for hist in self.histograms.itervalues():
            canvas = ROOT.TCanvas(hist.GetName(), hist.GetName())
            self.canvases.append(canvas)
            canvas.cd()
            hist.Draw()
            canvas.Update()
    
    def save_plots(self):
        """ Save the plots on disk
        """
        output_path = self.file_name[:-5]
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        for canvas in self.canvases:
            canvas.SaveAs("{}/{}.pdf".format(output_path, canvas.GetName()))
