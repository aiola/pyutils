#!/usr/bin/env python
# compare histograms

import math
import random
from enum import Enum
import ROOT
import root_utils
from physics import MeasuredQuantity

class CompareHistograms:
    """ Compare histograms
    """
    def __init__(self, name):
        self.name = name
        self.baseline_histogram = None
        self.histograms = None
        self.ratios = []
        self.fit_results = []
        self.opt_spectrum_baseline = ""
        self.opt_spectrum = ""
        self.opt_ratio = ""
        self.y_axis_ratio = "Ratio"
        self.do_spectra_plot = "logy"
        self.do_ratio_plot = "lineary"
        self.canvas_spectra = None
        self.canvas_ratio = None
        self.legend_spectra = None
        self.legend_ratio = None
        self.baseline_ratio = None
        self.marker_size = 1.2
        self.colors = [ROOT.kBlack, ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kAzure + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2, ROOT.kYellow + 2]
        self.markers = [ROOT.kOpenCircle, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kStar, ROOT.kFullCross, ROOT.kMultiply, ROOT.kPlus]
        self.lines = [1, 2, 9, 5, 7, 10, 4, 3, 6, 8]
        self.line_widths = [2] * 10
        self.fill_styles = [3001, 3002, 3003, 3004, 3005, 3006, 3007]
        self.main_histogram = None
        self.main_ratio_histogram = None
        self.max_spectrum = None
        self.min_spectrum = None
        self.max_ratio = None
        self.min_ratio = None
        self.results = None
        self.n_cols_leg_ratio = 1
        self.n_cols_leg_spectrum = 1
        self.x1_leg_ratio = 0.55
        self.x2_leg_ratio = 0.90
        self.y1_leg_ratio = 0.87
        self.x1_leg_spectrum = 0.55
        self.x2_leg_spectrum = 0.90
        self.y1_leg_spectrum = 0.87
        self.log_upper_space = 10  # this factor will be used to adjust the y axis in log scale
        self.log_lower_space = 2  # this factor will be used to adjust the y axis in log scale
        self.lin_upper_space = 0.9  # this factor will be used to adjust the y axis in linear scale
        self.lin_lower_space = 0.2  # this factor will be used to adjust the y axis in linear scale
        self.leg_text_size = 18
        self.leg_line_height = 0.06
        self.fixed_lower_ratio_bound = None

        self.baseline_for_ratio = None
        self.separate_baseline_uncertainty = False
        self.no_error_in_baseline = False
        self.ratio_relative_uncertainty = None
        self.ratio_relative_uncertainty_title = "Rel. Unc."
        self.grid_y_ratio = True

        self.grid_y_spectrum = False

        self.fit_function = "expo(0)+expo(2)"
        self.do_spectrum_legend = True
        self.do_ratio_legend = True
        self.units = ""
        self.minimum_limit = 0

    def set_ratio_relative_uncertainty_from_histogram(self, hist):
        self.ratio_relative_uncertainty = hist.Clone("{0}_unc".format(hist.GetName()))
        self.ratio_relative_uncertainty.SetTitle(self.ratio_relative_uncertainty_title)
        for ibin in range(0, self.ratio_relative_uncertainty.GetNbinsX() + 2):
            if self.ratio_relative_uncertainty.GetBinContent(ibin) != 0:
                self.ratio_relative_uncertainty.SetBinError(ibin, self.ratio_relative_uncertainty.GetBinError(ibin) / self.ratio_relative_uncertainty.GetBinContent(ibin))
                self.ratio_relative_uncertainty.SetBinContent(ibin, 1)
            else:
                self.ratio_relative_uncertainty.SetBinError(ibin, 0)
                self.ratio_relative_uncertainty.SetBinContent(ibin, 0)

    def add_stat(self, h):
        mean = MeasuredQuantity(h.GetMean(), h.GetMeanError(), self.units)
        sigma = MeasuredQuantity(h.GetStdDev(), h.GetStdDevError(), self.units)
        entry = self.legend_spectra.AddEntry(ROOT.nullptr, "#mu = {}".format(mean.to_string()), "")
        entry.SetTextSize(self.leg_text_size * 0.8)
        entry = self.legend_spectra.AddEntry(ROOT.nullptr, "#sigma = {}".format(sigma.to_string()), "")
        entry.SetTextSize(self.leg_text_size * 0.8)

    def prepare_spectra_canvas(self):
        if not self.canvas_spectra:
            print("Creating new canvas {0}".format(self.name))
            self.canvas_spectra = ROOT.TCanvas(self.name, self.name)

        if self.do_spectrum_legend:
            if self.do_spectrum_legend == "stat":
                self.leg_line_height *= 2.5
            if self.legend_spectra:
                y1 = self.legend_spectra.GetY1() - self.leg_line_height * (len(self.histograms) + 1) / self.n_cols_leg_spectrum
                if y1 < 0.2: y1 = 0.2
                self.legend_spectra.SetY1(y1)
            else:
                y1 = self.y1_leg_spectrum - self.leg_line_height * (len(self.histograms) + 1) / self.n_cols_leg_spectrum
                if y1 < 0.2: y1 = 0.2
                self.legend_spectra = ROOT.TLegend(self.x1_leg_spectrum, y1, self.x2_leg_spectrum, self.y1_leg_spectrum)
                self.legend_spectra.SetName("{0}_legend".format(self.canvas_spectra.GetName()))
                self.legend_spectra.SetNColumns(self.n_cols_leg_spectrum)
                self.legend_spectra.SetFillStyle(0)
                self.legend_spectra.SetBorderSize(0)
                self.legend_spectra.SetTextFont(43)
                self.legend_spectra.SetMargin(0.1)
                self.legend_spectra.SetTextSize(self.leg_text_size)

        if self.grid_y_spectrum:
            self.canvas_spectra.SetGridy()

        if "hist" in self.opt_spectrum_baseline:
            self.baseline_histogram.SetLineColor(self.colors[0])
            self.baseline_histogram.SetLineWidth(self.line_widths[0])
            self.baseline_histogram.SetLineStyle(self.lines[0])
            if self.do_spectrum_legend:
                self.legend_spectra.AddEntry(self.baseline_histogram, self.baseline_histogram.GetTitle(), "l")
                if self.do_spectrum_legend == "stat":
                    self.add_stat(self.baseline_histogram)
        elif "e2" in self.opt_spectrum_baseline:
            self.baseline_histogram.SetLineColor(self.colors[0])
            self.baseline_histogram.SetFillColor(self.colors[0])
            self.baseline_histogram.SetLineWidth(1)
            self.baseline_histogram.SetLineStyle(self.lines[0])
            self.baseline_histogram.SetFillStyle(self.fill_styles[0])
            if self.do_spectrum_legend:
                self.legend_spectra.AddEntry(self.baseline_histogram, self.baseline_histogram.GetTitle(), "f")
                if self.do_spectrum_legend == "stat":
                    self.add_stat(self.baseline_histogram)
        else:
            self.baseline_histogram.SetMarkerColor(self.colors[0])
            self.baseline_histogram.SetLineColor(self.colors[0])
            self.baseline_histogram.SetMarkerStyle(self.markers[0])
            self.baseline_histogram.SetMarkerSize(self.marker_size)
            if self.do_spectrum_legend:
                self.legend_spectra.AddEntry(self.baseline_histogram, self.baseline_histogram.GetTitle(), "pe")
                if self.do_spectrum_legend == "stat":
                    self.add_stat(self.baseline_histogram)

        if isinstance(self.baseline_histogram, ROOT.TGraph):
            if not "a" in self.opt_spectrum_baseline or not "A" in self.opt_spectrum_baseline:
                self.opt_spectrum_baseline += "A"
        
        if len(self.baseline_histogram.GetListOfFunctions()) > 0:
            for obj in self.baseline_histogram.GetListOfFunctions():
                if isinstance(obj, ROOT.TF1):
                    obj.SetLineColor(self.colors[0])
                    obj.SetLineStyle(self.lines[0])
                    obj.SetLineWidth(self.line_widths[0])

        print("Plotting histogram '{0}' with option '{1}'".format(self.baseline_histogram.GetName(), self.opt_spectrum_baseline))
        self.canvas_spectra.cd()
        self.baseline_histogram.Draw(self.opt_spectrum_baseline)

        if "frac" in self.baseline_histogram.GetYaxis().GetTitle():
            self.canvas_spectra.SetLeftMargin(0.12)
            self.baseline_histogram.GetYaxis().SetTitleOffset(1.4)
        if "same" not in self.opt_spectrum:
            self.opt_spectrum += "same"

        if not self.main_histogram:
            if isinstance(self.baseline_histogram, ROOT.TH1):
                self.main_histogram = self.baseline_histogram
            elif isinstance(self.baseline_histogram, ROOT.TGraph):
                self.main_histogram = self.baseline_histogram.GetHistogram()
            else:
                print("Type of object '{}' not recognized!".format(self.baseline_histogram))
        self.baseline_for_ratio = self.baseline_histogram.Clone("{0}_copy".format(self.baseline_histogram.GetName()))

        minimum = root_utils.find_minimum(self.baseline_histogram, self.minimum_limit, "hist" not in self.opt_spectrum_baseline)
        if not minimum is None:
            if self.min_spectrum is None:
                self.min_spectrum = minimum
            else:
                self.min_spectrum = min(self.min_spectrum, minimum)
        maximum = root_utils.find_maximum(self.baseline_histogram, self.minimum_limit, "hist" not in self.opt_spectrum_baseline)
        if not maximum is None:
            if self.max_spectrum is None:
                self.max_spectrum = maximum
            else:
                self.max_spectrum = max(self.max_spectrum, maximum)

    def prepare_ratio_canvas(self):
        cname = "{0}_Ratio".format(self.name)
        if not self.canvas_ratio:
            self.canvas_ratio = ROOT.TCanvas(cname, cname)
        self.canvas_ratio.cd()

        n = len(self.histograms)
        if self.ratio_relative_uncertainty or self.separate_baseline_uncertainty:
            n += 1

        if self.do_ratio_legend:
            if self.legend_ratio:
                y1 = self.legend_ratio.GetY1() - self.leg_line_height * n / self.n_cols_leg_ratio
                if y1 < 0.2: y1 = 0.2
                self.legend_ratio.SetY1(y1)
            else:
                y1 = self.y1_leg_ratio - self.leg_line_height * n / self.n_cols_leg_ratio
                if y1 < 0.2: y1 = 0.2
                self.legend_ratio = ROOT.TLegend(self.x1_leg_ratio, y1, self.x2_leg_ratio, self.y1_leg_ratio)
                self.legend_ratio.SetName("{0}_legend".format(self.canvas_ratio.GetName()))
                self.legend_ratio.SetNColumns(self.n_cols_leg_ratio)
                self.legend_ratio.SetFillStyle(0)
                self.legend_ratio.SetBorderSize(0)
                self.legend_ratio.SetTextFont(43)
                self.legend_ratio.SetTextSize(self.leg_text_size)

        if self.grid_y_ratio:
            self.canvas_ratio.SetGridy()

        if self.separate_baseline_uncertainty or self.no_error_in_baseline:
            for ibin in range(0, self.baseline_for_ratio.GetNbinsX() + 2):
                self.baseline_for_ratio.SetBinError(ibin, 0)

        if self.separate_baseline_uncertainty:
            self.set_ratio_relative_uncertainty_from_histogram(self.baseline_histogram)
            opt = "e2"
            if "same" in self.opt_ratio:
                opt += "same"
            h = self.ratio_relative_uncertainty.DrawCopy(opt)
            h.SetFillColor(self.colors[0])
            h.SetFillStyle(self.fill_styles[0])
            h.SetLineColor(self.colors[0])
            h.GetYaxis().SetTitle(self.y_axis_ratio)
            if self.do_ratio_legend: self.legend_ratio.AddEntry(h, h.GetTitle(), "f")
            self.results.append(h)
            if "same" not in self.opt_ratio:
                self.opt_ratio += "same"
            if not self.main_ratio_histogram:
                self.main_ratio_histogram = h
            minimum = root_utils.find_minimum(h, self.minimum_limit, True)
            if not minimum is None:
                if self.min_ratio is None:
                    self.min_ratio = minimum
                else:
                    self.min_ratio = min(self.min_ratio, minimum)
            maximum = root_utils.find_maximum(h, self.minimum_limit, True)
            if not maximum is None:
                if self.max_ratio is None:
                    self.max_ratio = maximum
                else:
                    self.max_ratio = max(self.max_ratio, maximum)

    def plot_histogram(self, color, marker, line, lwidth, h):
        minimum = root_utils.find_minimum(h, self.minimum_limit, not "hist" in self.opt_spectrum)
        if not minimum is None:
            if self.min_spectrum is None:
                self.min_spectrum = minimum
            else:
                self.min_spectrum = min(self.min_spectrum, minimum)
        maximum = root_utils.find_maximum(h, self.minimum_limit, not "hist" in self.opt_spectrum)
        if not maximum is None:
            if self.max_spectrum is None:
                self.max_spectrum = maximum
            else:
                self.max_spectrum = max(self.max_spectrum, maximum)

        print("Plotting histogram '{0}' with option '{1}'".format(h.GetName(), self.opt_spectrum))
        self.canvas_spectra.cd()
        h.Draw(self.opt_spectrum)

        if "hist" in self.opt_spectrum:
            h.SetLineColor(color)
            h.SetLineWidth(lwidth)
            h.SetLineStyle(line)
            if self.do_spectrum_legend:
                self.legend_spectra.AddEntry(h, h.GetTitle(), "l")
                if self.do_spectrum_legend == "stat":
                    self.add_stat(h)
        else:
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetMarkerStyle(marker)
            h.SetMarkerSize(self.marker_size)
            if self.do_spectrum_legend:
                self.legend_spectra.AddEntry(h, h.GetTitle(), "pe")
                if self.do_spectrum_legend == "stat":
                    self.add_stat(h)

        if len(h.GetListOfFunctions()) > 0:
            for obj in h.GetListOfFunctions():
                if isinstance(obj, ROOT.TF1):
                    obj.SetLineColor(color)
                    obj.SetLineStyle(line)
                    obj.SetLineWidth(lwidth)

    def fit_and_make_consistent(self, h, templateH):
        fit_func = ROOT.TF1("{0}_fit".format(h.GetName()), self.fit_function, h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        fit_func.SetParameter(0, 1)
        fit_func.SetParameter(1, -0.5)
        fit_func.SetParameter(2, 0)
        fit_func.SetParameter(3, 0)
        fit_func.SetParameter(0, 1. / fit_func.Eval(h.GetXaxis().GetBinCenter(1)))
        fit_func.SetParameter(2, 1)
        fit_func.SetParameter(3, -0.2)
        fit_func.SetParameter(2, 1. / fit_func.Eval(h.GetXaxis().GetBinCenter(int(h.GetNbinsX() / 2))))
        fitR = h.Fit(fit_func, "NS")
        fitOk = int(fitR)
        if not fitOk == 0: return None
        h_fit = root_utils.soft_clone(templateH, "{0}_fith".format(h.GetName()))
        for ibin in range(1, h_fit.GetNbinsX() + 1):
            valErr = fit_func.IntegralError(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
            val = fit_func.Integral(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
            print("integral = {0:.5f}, central = {1:.5f}".format(val, fit_func.Eval((h_fit.GetXaxis().GetBinCenter(ibin)))))
            h_fit.SetBinContent(ibin, val)
            h_fit.SetBinError(ibin, valErr)
        self.fit_results.append(fit_func)
        self.fit_results.append(h_fit)
        fit_func.SetLineColor(h.GetLineColor())
        self.canvas_spectra.cd()
        fit_func.Draw("same")
        return h_fit

    def rebin_and_make_consistent(self, h, templateH):
        return root_utils.rebin_1D(h, templateH.GetXaxis())

    def plot_ratio(self, color, marker, line, lwidth, h):
        compBinning = root_utils.AxisCompare.check_consistency(h.GetXaxis(), self.baseline_for_ratio.GetXaxis())
        print("Result of the binning comparison between {} and {} is: {}".format(h.GetName(), self.baseline_for_ratio.GetName(), compBinning))
        if compBinning == root_utils.AxisCompare.Identical:
            hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
        elif compBinning == root_utils.AxisCompare.ContainsSameBinning or compBinning == root_utils.AxisCompare.IsContainedSameBinning or compBinning == root_utils.AxisCompare.OverlapsSameBinning:
            print("Trying to rebin histogram {0}".format(h.GetName()))
            hRatio = self.rebin_and_make_consistent(h, self.baseline_for_ratio)
            if not hRatio:
                print("Rebin unsuccessfull!")
                return
        elif compBinning == root_utils.AxisCompare.Contains or compBinning == root_utils.AxisCompare.IsContained or compBinning == root_utils.AxisCompare.Overlaps:
            print("Need to rebin.")
            bins = "["
            for x in h.GetXaxis().GetXbins(): bins += "{}, ".format(x)
            bins = bins[:-2]
            bins += "]"
            print("Original binning: {}".format(bins))
            bins = "["
            for x in self.baseline_for_ratio.GetXaxis().GetXbins(): bins += "{}, ".format(x)
            bins = bins[:-2]
            bins += "]"
            print("Final binning: {}".format(bins))
            print("Trying to fit histogram {0} with function {1}".format(h.GetName(), self.fit_function))
            hRatio = self.fit_and_make_consistent(h, self.baseline_for_ratio)
            if not hRatio:
                print("Fit unsuccessfull!")
                return
        elif compBinning == root_utils.AxisCompare.NoOverlap:
            print("The two histograms {}, {} have no overlap. Unable to generate a ratio.".format(h.GetName(), self.baseline_for_ratio.GetName()))
            return
        else:
            print("compare_histograms, plot_ration: Should not end up here!")
            exit(1)

        hRatio.GetYaxis().SetTitle(self.y_axis_ratio)
        if not self.baseline_ratio:
            self.baseline_ratio = hRatio
        if "hist" in self.opt_ratio:
            hRatio.SetLineColor(color)
            hRatio.SetLineWidth(lwidth)
            hRatio.SetLineStyle(line)
            if self.do_ratio_legend: self.legend_ratio.AddEntry(hRatio, h.GetTitle(), "l")
        else:
            hRatio.SetMarkerColor(color)
            hRatio.SetLineColor(color)
            hRatio.SetMarkerStyle(marker)
            hRatio.SetLineStyle(1)
            hRatio.SetMarkerSize(self.marker_size)
            if self.do_ratio_legend: self.legend_ratio.AddEntry(hRatio, h.GetTitle(), "pe")
        self.ratios.append(hRatio)
        hRatio.SetTitle("{0} Ratio".format(h.GetTitle()))
        hRatio.Divide(self.baseline_for_ratio)
        self.canvas_ratio.cd()
        hRatio.Draw(self.opt_ratio)
        if not self.main_ratio_histogram:
            self.main_ratio_histogram = hRatio
        minimum = root_utils.find_minimum(hRatio, self.minimum_limit, not "hist" in self.opt_ratio)
        if not minimum is None:
            if self.min_ratio is None:
                self.min_ratio = minimum
            else:
                self.min_ratio = min(self.min_ratio, minimum)
        maximum = root_utils.find_maximum(hRatio, self.minimum_limit, not "hist" in self.opt_ratio)
        if not maximum is None:
            if self.max_ratio is None:
                self.max_ratio = maximum
            else:
                self.max_ratio = max(self.max_ratio, maximum)

        if not "same" in self.opt_ratio:
            self.opt_ratio += " same"

    def compare_spectra(self, baseline, histos):
        while len(histos) + 1 > len(self.colors): self.colors += random.sample(self.colors[1:], len(self.colors) - 1)
        while len(histos) + 1 > len(self.markers): self.markers += random.sample(self.markers[1:], len(self.markers) - 1)
        while len(histos) + 1 > len(self.lines): self.lines += random.sample(self.lines[1:], len(self.lines) - 1)
        while len(histos) + 1 > len(self.line_widths): self.line_widths += random.sample(self.line_widths[1:], len(self.line_widths) - 1)
        while len(histos) + 1 > len(self.fill_styles): self.fill_styles += random.sample(self.fill_styles[1:], len(self.fill_styles) - 1)
        self.results = []
        print("compare_spectra: {0}".format(self.name))
        self.baseline_histogram = baseline
        self.histograms = histos
        if not isinstance(baseline, ROOT.TH1) and self.do_ratio_plot:
            print("Ratio is only available for histograms. Option is disabled.")
            self.do_ratio_plot = False
        print("Baseline: {0}".format(self.baseline_histogram.GetName()))
        for s in self.histograms:
            print(s.GetName())

        if self.do_spectra_plot:
            self.prepare_spectra_canvas()

        if self.do_ratio_plot:
            self.prepare_ratio_canvas()

        for color, marker, line, linew, h in zip(self.colors[1:], self.markers[1:], self.lines[1:], self.line_widths[1:], self.histograms):
            if self.do_spectra_plot:
                self.plot_histogram(color, marker, line, linew, h)
            if self.do_ratio_plot:
                self.plot_ratio(color, marker, line, linew, h)
        self.adjust_y_limits()
        self.generate_results()
        if self.main_histogram not in self.histograms:
            self.results.append(self.main_histogram)
        return self.results

    def compare_uncertainties(self, baseline, histos):
        baseline_unc = root_utils.get_relative_uncertainty(baseline)
        histos_unc = [root_utils.get_relative_uncertainty(h) for h in histos]
        self.CompareSpectra(baseline_unc, histos_unc)
        self.results.append(baseline_unc)
        self.results.extend(histos_unc)
        return self.results

    def generate_results(self):
        self.results.extend(self.ratios)
        self.results.extend(self.fit_results)
        if self.canvas_spectra:
            self.results.append(self.canvas_spectra)
        if self.canvas_ratio:
            self.results.append(self.canvas_ratio)
        if self.legend_spectra:
            self.results.append(self.legend_spectra)
        if self.legend_ratio:
            self.results.append(self.legend_ratio)
        if self.ratio_relative_uncertainty:
            self.results.append(self.ratio_relative_uncertainty)

    def adjust_y_limits(self):
        if not self.max_ratio is None and not self.min_ratio is None and self.do_ratio_plot:
            print("Adjusting y limits for Ratio")
            if self.min_ratio <= 0 and self.do_ratio_plot == "logy":
                print("{}: requested logy ratio, but minimum is <= 0. Switching to linear scale.".format(self.name))
                self.do_ratio_plot = "lineary"
            if self.do_ratio_plot == "logy":
                max = self.max_ratio * self.log_upper_space
                min = self.min_ratio / self.log_lower_space
                self.canvas_ratio.SetLogy()
            else:
                max = self.max_ratio + (self.max_ratio - self.min_ratio) * self.lin_upper_space
                min = self.min_ratio - (self.max_ratio - self.min_ratio) * self.lin_lower_space
                if min < 0 and self.min_ratio > 0: min = 0
            self.main_ratio_histogram.GetYaxis().UnZoom()
            self.main_ratio_histogram.SetMinimum(min)
            self.main_ratio_histogram.SetMaximum(max)
            if self.do_ratio_legend:
                self.canvas_ratio.cd()
                self.legend_ratio.Draw()
        if not self.fixed_lower_ratio_bound is None:
            self.main_ratio_histogram.SetMinimum(self.fixed_lower_ratio_bound)

        if not self.max_spectrum is None and not self.min_spectrum is None and self.do_spectra_plot:
            print("Adjusting y limits for Spectrum")
            if self.min_spectrum <= 0 and self.do_spectra_plot == "logy":
                print("{}: requested logy spectra, but minimum is <= 0. Switching to linear scale.".format(self.name))
                self.do_spectra_plot = "lineary"
            if self.do_spectra_plot == "logy":
                max = self.max_spectrum * self.log_upper_space
                min = self.min_spectrum / self.log_lower_space
                self.canvas_spectra.SetLogy()
            else:
                max = self.max_spectrum + (self.max_spectrum - self.min_spectrum) * self.lin_upper_space
                min = self.min_spectrum - (self.max_spectrum - self.min_spectrum) * self.lin_lower_space
                if min < 0 and self.min_spectrum > 0: min = 0
            self.main_histogram.SetMinimum(min)
            self.main_histogram.SetMaximum(max)
            if self.do_spectrum_legend:
                self.canvas_spectra.cd()
                self.legend_spectra.Draw()
