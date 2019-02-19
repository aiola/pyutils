#!/usr/bin/env python
# pyroot utilities

import math
import array
from collections import OrderedDict
from enum import Enum
import ROOT

class AxisCompare(Enum):
    """ Compare axis of two ROOT histograms
    """ 
    Identical = 0
    SameLimits = 1
    IsContainedSameBinning = 2
    ContainsSameBinning = 3
    OverlapsSameBinning = 4
    IsContained = 5
    Contains = 6
    Overlaps = 7
    NoOverlap = 8

    @classmethod
    def check_consistency(cls, axis1, axis2):
        """ Compares the axisi of two histograms
        """
        is_contained = False
        contains = False
        overlaps = False
        if axis1.GetBinLowEdge(1) <= axis2.GetBinLowEdge(1):
            if axis1.GetBinUpEdge(axis1.GetNbins()) >= axis2.GetBinUpEdge(axis2.GetNbins()):
                contains = True
            else:
                overlaps = True
        if axis2.GetBinLowEdge(1) <= axis1.GetBinLowEdge(1):
            if axis2.GetBinUpEdge(axis2.GetNbins()) >= axis1.GetBinUpEdge(axis1.GetNbins()):
                is_contained = True
            else:
                overlaps = True

        if not contains and not is_contained and not overlaps:
            return cls.NoOverlap

        same_binning = True
        for ibin1 in range(1, axis1.GetNbins() + 1):
            if axis1.GetBinLowEdge(ibin1) >= axis2.GetBinLowEdge(1): 
                break
            ibin1 += 1

        for ibin2 in range(1, axis2.GetNbins() + 1):
            if axis2.GetBinLowEdge(ibin2) >= axis1.GetBinLowEdge(1): 
                break
            ibin2 += 1

        while(ibin1 <= axis1.GetNbins() and ibin2 <= axis2.GetNbins()):
            if axis1.GetBinLowEdge(ibin1) != axis2.GetBinLowEdge(ibin2):
                same_binning = False
                break
            if axis1.GetBinUpEdge(ibin1) != axis2.GetBinUpEdge(ibin2):
                same_binning = False
                break
            ibin1 += 1
            ibin2 += 1

        if contains and is_contained:
            if same_binning:
                return cls.Identical
            else:
                return cls.SameLimits
        elif contains:
            if same_binning:
                return cls.ContainsSameBinning
            else:
                return cls.Contains
        elif is_contained:
            if same_binning:
                return cls.IsContainedSameBinning
            else:
                return cls.IsContained
        else:
            if same_binning:
                return cls.OverlapsSameBinning
            else:
                return cls.Overlaps


def soft_clone(origin, name, title=None, y_axis_title=None):
    """ Makes a soft clone of a histogram. It copies only the structure but not the content
    """
    if not title:
        title = name
    if not y_axis_title:
        y_axis_title = origin.GetYaxis().GetTitle()
    new_histogram = ROOT.TH1D(name, title, origin.GetNbinsX(), origin.GetXaxis().GetXbins().GetArray())
    new_histogram.GetXaxis().SetTitle(origin.GetXaxis().GetTitle())
    new_histogram.GetYaxis().SetTitle(y_axis_title)
    return new_histogram

def get_relative_uncertainty(origin):
    """ Creates a new histogram with the relative uncertainty from another histogram
    """
    relative_uncertainty = soft_clone(origin, "{0}_unc".format(origin.GetName()), "{0} Rel. Unc.".format(origin.GetTitle()), "rel. unc.")
    for ibin in xrange(0, origin.GetNbinsX() + 2):
        if origin.GetBinContent(ibin) == 0:
            continue
        relative_uncertainty.SetBinContent(ibin, origin.GetBinError(ibin) / origin.GetBinContent(ibin))
    return relative_uncertainty

def get_root_object(obj, name):
    """ Obtain a ROOT object form a ROOT file or collection
    """
    slash = 0
    while slash >= 0:
        slash = name.find("/", slash)
        if name[slash + 1] == '/':
            slash += 2
        else:
            break

    if slash < 0:
        name_lookup = name.replace("//", "/")
        name = None
    else:
        name_lookup = name[:slash].replace("//", "/")
        name = name[slash + 1:]
    res = None
    name_obj = None
    if isinstance(obj, ROOT.TCollection):
        res = obj.FindObject(name_lookup)
        name_obj = obj.GetName()
    elif isinstance(obj, dict) or isinstance(obj, OrderedDict):
        res = obj[name_lookup]
        name_obj = obj
    elif isinstance(obj, ROOT.TDirectory):
        res = obj.Get(name_lookup)
        name_obj = obj.GetName()
    if not res:
        if name_obj:
            print("Could not find object {0} in collection '{1}'".format(name_lookup, name_obj))
            if isinstance(obj, ROOT.TObject): obj.ls()
        else:
            print("Could not find object {} in following collection".format(name_lookup))
            print(obj)
        return None
    if name:
        return get_root_object(res, name)
    else:
        return res

def get_object_and_merge(file_list, name):
    """ Obtains objects from a list of files and merges them
    """
    res = None
    for current_file in file_list:
        obj = get_root_object(current_file, name)
        if res:
            res.Add(obj)
        else:
            res = obj.Clone()
    return res

def generate_multi_canvas(cname, n_pads):
    """ Generate a ROOT canvas and divides it in several pads
    """
    rows = int(math.floor(math.sqrt(n_pads)))
    cols = int(math.ceil(float(n_pads) / rows))
    canvas = ROOT.TCanvas(cname, cname, cols * 400, rows * 400)
    canvas.Divide(cols, rows)
    return canvas

def find_minimum(histogram, limit=0., errors=True):
    """ Finds the minimum in a histogram, taking into account also errors
    Useful for plotting
    """
    minimum = None
    if histogram.GetDimension() == 1:
        for ibin in xrange(1, histogram.GetNbinsX() + 1):
            if errors:
                cont = histogram.GetBinContent(ibin) - histogram.GetBinError(ibin)
            else:
                cont = histogram.GetBinContent(ibin)
            if cont <= limit:
                continue
            if minimum is None or cont < minimum:
                minimum = cont
    elif histogram.GetDimension() == 2:
        for xbin in xrange(1, histogram.GetNbinsX() + 1):
            for ybin in xrange(1, histogram.GetNbinsY() + 1):
                if errors:
                    cont = histogram.GetBinContent(xbin, ybin) - histogram.GetBinError(xbin, ybin)
                else:
                    cont = histogram.GetBinContent(xbin, ybin)
                if cont <= limit:
                    continue
                if minimum is None or cont < minimum:
                    minimum = cont
    return minimum

def find_maximum(histogram, limit=0., errors=True):
    """ Finds the maximum in a histogram, taking into account also errors
    Useful for plotting
    """
    maximum = None
    if histogram.GetDimension() == 1:
        for ibin in xrange(1, histogram.GetNbinsX() + 1):
            if errors:
                cont = histogram.GetBinContent(ibin) + histogram.GetBinError(ibin)
            else:
                cont = histogram.GetBinContent(ibin)
            if cont <= limit:
                continue
            if maximum is None or cont > maximum:
                maximum = cont
    elif histogram.GetDimension() == 2:
        for xbin in xrange(1, histogram.GetNbinsX() + 1):
            for ybin in xrange(1, histogram.GetNbinsY() + 1):
                if errors:
                    cont = histogram.GetBinContent(xbin, ybin) + histogram.GetBinError(xbin, ybin)
                else:
                    cont = histogram.GetBinContent(xbin, ybin)
                if cont <= limit: continue
                if maximum is None or cont > maximum:
                    maximum = cont
    return maximum

def divide_no_errors(ratio, den):
    """ Divides two histograms, ignoring errors
    """
    if not ratio.GetNbinsX() == den.GetNbinsX():
        print("DMesonJetUtils.DivideNoErrors: histograms have different number of bins!")
        return False

    for ibin in xrange(0, ratio.GetNbinsX() + 2):
        if den.GetBinContent(ibin) == 0:
            ratio.SetBinContent(ibin, 0)
        else:
            ratio.SetBinContent(ibin, ratio.GetBinContent(ibin) / den.GetBinContent(ibin))

    return True

def V2TH1(vect):
    """ Transforms a vector into a TH1
    """
    result = ROOT.TH1D("vect", "vect", len(vect) - 2, 1, len(vect) - 2)
    for ibin in xrange(0, result.GetNbinsX() + 2):
        result.SetBinContent(ibin, vect[ibin])
    return result

def build_histogram(axis, name, y_axis_title):
    """ Builds a histogram given an axis
    """
    if len(axis) == 1:
        hist = ROOT.TH1D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(y_axis_title)
        hist.Sumw2()
    elif len(axis) == 2:
        hist = ROOT.TH2D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(y_axis_title)
        hist.Sumw2()
    else:
        hist = ROOT.TH3D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins), len(axis[2].fBins) - 1, array.array('d', axis[2].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(axis[2].GetTitle())
        hist.Sumw2()
    return hist

def rebin_1D(hist, xaxis, warnings=False):
    """ Rebin a histogram
    """
    return rebin_1D_from_bins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), warnings)

def rebin_1D_from_bins(hist, name, nbinsX, binsX, warnings=False):
    """ Rebin a histogram
    """
    axis = ROOT.TAxis(nbinsX, binsX)
    compare_axis = AxisCompare.check_consistency(hist.GetXaxis(), axis)
    if compare_axis == AxisCompare.Identical:
        result = hist.Clone(name)
        return result
    result = ROOT.TH1D(name, name, nbinsX, binsX)
    result.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    result.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    for xbin in xrange(0, hist.GetXaxis().GetNbins() + 2):
        xbin_center = hist.GetXaxis().GetBinCenter(xbin)
        rxbin = result.GetXaxis().FindBin(xbin_center)
        bin_value = hist.GetBinContent(xbin) + result.GetBinContent(rxbin)
        bin_error = math.sqrt(hist.GetBinError(xbin) ** 2 + result.GetBinError(rxbin) ** 2)
        result.SetBinContent(rxbin, bin_value)
        result.SetBinError(rxbin, bin_error)
        if bin_value > 0:
            rel_err = bin_error / bin_value
            if rel_err > 0.9 and warnings:
                print("Bin ({0}) has rel stat err = {1}. This is VERY dangerous!".format(xbin, rel_err))
    return result

def rebin_2D(hist, xaxis, yaxis, warnings=False):
    """ Rebin a histogram
    """
    return rebin_2D_from_bins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), yaxis.GetNbins(), yaxis.GetXbins().GetArray(), warnings)

def rebin_2D_from_bins(hist, name, nbinsX, binsX, nbinsY, binsY, warnings=False):
    """ Rebin a histogram
    """
    xaxis = ROOT.TAxis(nbinsX, binsX)
    x_compare_axis = AxisCompare.check_consistency(hist.GetXaxis(), xaxis)
    yaxis = ROOT.TAxis(nbinsY, binsY)
    y_compare_axis = AxisCompare.check_consistency(hist.GetYaxis(), yaxis)
    if x_compare_axis == AxisCompare.Identical and y_compare_axis == AxisCompare.Identical:
        result = hist.Clone(name)
        return result
    result = ROOT.TH2D(name, name, nbinsX, binsX, nbinsY, binsY)
    result.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    result.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    result.GetZaxis().SetTitle(hist.GetZaxis().GetTitle())
    result.Sumw2()
    hist_getxaxis_getbincenter = hist.GetXaxis().GetBinCenter
    hist_getyaxis_getbincenter = hist.GetYaxis().GetBinCenter
    result_getxaxis_getbinupedge = result.GetXaxis().GetBinUpEdge
    result_getyaxis_getbinupedge = result.GetYaxis().GetBinUpEdge
    hist_at = hist.At
    result_at = result.At
    result_setat = result.SetAt
    if hist.GetSumw2().GetSize() > 0:
        hist_getsumw2_at = hist.GetSumw2().At
    else:
        hist_getsumw2_at = hist.At
    result_getsumw2_at = result.GetSumw2().At
    result_getsumw2_setat = result.GetSumw2().SetAt
    hist_getxaxis_getnbins_2 = hist.GetXaxis().GetNbins() + 2
    hist_getyaxis_getnbins_2 = hist.GetYaxis().GetNbins() + 2
    hist_bin = 0
    rybin = 0
    result_addbincontent = result.AddBinContent
    for ybin in xrange(0, hist_getyaxis_getnbins_2):
        y_bin_center = hist_getyaxis_getbincenter(ybin)
        while (y_bin_center > result_getyaxis_getbinupedge(rybin) and rybin < nbinsY + 1):
            rybin += 1
        rxbin = 0
        for xbin in xrange(0, hist_getxaxis_getnbins_2):
            xbin_center = hist_getxaxis_getbincenter(xbin)
            while (xbin_center > result_getxaxis_getbinupedge(rxbin) and rxbin < nbinsX + 1):
                rxbin += 1
            r_bin = rxbin + (nbinsX + 2) * rybin
            result_setat(result_at(r_bin) + hist_at(hist_bin), r_bin)
            result_getsumw2_setat(result_getsumw2_at(r_bin) + hist_getsumw2_at(hist_bin), r_bin)
            hist_bin += 1

    result.SetEntries(hist.GetEntries())

    if warnings:
        for xbin in xrange(0, result.GetXaxis().GetNbins() + 2):
            for ybin in xrange(0, result.GetYaxis().GetNbins() + 2):
                bin_value = result.GetBinContent(xbin, ybin)
                bin_error = result.GetBinError(xbin, ybin)
                if bin_value > 0:
                    rel_err = bin_error / bin_value
                    if rel_err > 0.9:
                        print("Bin ({0},{1}) has rel stat err = {2}. This is VERY dangerous!".format(xbin, ybin, rel_err))
    return result

def frange(start, stop, step, closed=False):
    """ Float range
    """
    i = start
    if closed:
        stop += step
    while i < stop:
        yield i
        i += step

def fit_and_rebin(hstat, new_axis):
    """ Fit histogram and rebin
    """
    results = []
    xmin = new_axis.GetBinLowEdge(1) - 5
    if xmin < 5:
        xmin = 5
    xmax = new_axis.GetBinUpEdge(new_axis.GetNbins()) + 5
    fit_func = ROOT.TF1("myfit", "expo(0)+expo(2)", xmin, xmax)
    results.append(fit_func)
    cname = "Fit_{}".format(hstat.GetName())
    canvas = ROOT.TCanvas(cname, cname)
    canvas.SetLogy()
    results.append(canvas)
    hcopy = hstat.DrawCopy()
    hcopy.GetXaxis().SetRangeUser(new_axis.GetBinLowEdge(1), new_axis.GetBinUpEdge(new_axis.GetNbins()))
    results.append(hcopy)
    fit_func.SetParameter(0, 1)
    fit_func.SetParameter(1, -1)
    fit_func.SetParameter(2, -99999999)
    fit_func.SetParameter(3, 0)
    fit_func.SetParameter(0, 1 + math.log(hstat.GetBinContent(hstat.GetXaxis().FindBin(xmin)) / fit_func.Eval(xmin)))
    fit_func.SetParameter(2, 1)
    fit_func.SetParameter(3, -0.5)
    fit_func.SetParameter(2, 1 + math.log(hstat.GetBinContent(hstat.GetXaxis().FindBin((xmax - xmin) / 2)) / fit_func.Eval((xmax - xmin) / 2)))
    fit_result = hcopy.Fit(fit_func, "S", "", xmin, xmax)
    fit_status = int(fit_result)
    if fit_status != 0:
        print("The fit was unsuccessfull!")
        return None
    h_fit = ROOT.TH1D("{}_rebinned".format(hstat.GetName()), hstat.GetTitle(), new_axis.GetNbins(), new_axis.GetXbins().GetArray())
    results.append(h_fit)
    for ibin in range(1, h_fit.GetNbinsX() + 1):
        valErr = fit_func.IntegralError(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
        val = fit_func.Integral(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
        print("integral = {0:.5f}, central = {1:.5f}".format(val, fit_func.Eval((h_fit.GetXaxis().GetBinCenter(ibin)))))
        h_fit.SetBinContent(ibin, val)
        h_fit.SetBinError(ibin, valErr)
    h_fit.Draw("same")

    ratio = hcopy.Clone("ratio")
    results.append(ratio)
    ratio.Divide(fit_func)
    cname = "FitRatio_{}".format(hstat.GetName())
    canvas_ratio = ROOT.TCanvas(cname, cname)
    results.append(canvas_ratio)
    ratio.Draw()
    return h_fit, results


def is_variable_bin_size(histo):
    if not histo.GetXaxis().IsVariableBinSize():
        return False
    binWidth = histo.GetXaxis().GetBinWidth(1)
    for ibin in range(2, histo.GetNbinsX() + 1):
        if histo.GetXaxis().GetBinWidth(ibin) != binWidth:
            return False
    return True