""" This module contains utilities for fitting distributions.
"""

import math
import inspect
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import physics

# Function library
def gaus(x, C, mu, sigma):
  x = np.atleast_1d(x)
  return C * np.exp(-(x-mu)**2/(2*sigma**2))

def double_gaus(x, C1, mu1, sigma1, C2, mu2, sigma2):
  return gaus(x, C1, mu1, sigma1)+gaus(x, C2, mu2, sigma2)

def gaus_norm(x, N, mu, sigma):
  x = np.atleast_1d(x)
  return N * gaus(x, 1., mu, sigma) / sigma / np.sqrt(2 * math.pi)

def double_gaus_norm(x, N1, mu1, sigma1, N2, mu2, sigma2):
  return gaus_norm(x, N1, mu1, sigma1)+gaus_norm(x, N2, mu2, sigma2)

def double_gaus_norm_single_mean(x, mu, N, sigma1, f, sigma2):
  return gaus_norm(x, N*f, mu, sigma1) + gaus_norm(x, N*(1.-f), mu, sigma2)

def crystal_ball(x, C, mu, sigma, alpha, n):
  x = np.atleast_1d(x)
  nsigmas = (x - mu) / sigma
  A = np.power(n / np.abs(alpha), n) * np.exp(-alpha**2/2)
  B = n / np.abs(alpha) - np.abs(alpha)
  if alpha > 0:
      condition = (nsigmas < -alpha)
      part_1 = C * A * np.power((B - nsigmas[condition]), -n)
      part_2 = gaus(x[~condition], C, mu, sigma)
  else:
      condition = (nsigmas < -alpha)
      part_1 = gaus(x[condition], C, mu, sigma)
      part_2 = C * A * np.power((B + nsigmas[~condition]), -n)
  return np.concatenate((part_1, part_2))

def polynomial(x, coeff):
  x = np.atleast_1d(x)
  xpower = np.ones_like(x)
  result = 0.
  for c in coeff:
    result += c * xpower
    xpower *= x
  return result

def constant(x, C):
  return polynomial(x, [C])

def linear(x, a, b):
  return polynomial(x, [b,a])

def quadratic(x, a, b, c):
  return polynomial(x, [c,b,a])

def cubic(x, a, b, c, d):
  return polynomial(x, [d,c,b,a])

class fitter(object):
  def __init__(self, function_name, function):
    self.data = None
    self.bins = None
    self.function_name = function_name
    self.function = function
    self.density = False
    self.bounds = (-np.inf,np.inf)
    self.param_names = None
    self.histogram = None
    self.params = None
    self.cov = None
  
  def fit(self, p0):
    if self.param_names is None:
      self.param_names = [ 'p{}'.format(n) for n in range(0,len(p0)) ]
    if self.histogram is None:
      if self.data is not None:
        self.histogram = generate_histogram(self.data, self.bins, self.density)
      else:
        print("Error: no data!")
        return False
    self.bins = self.histogram[0]
    try:
      result = curve_fit(self.function, self.histogram[0][self.histogram[1]>0], self.histogram[1][self.histogram[1]>0], sigma=self.histogram[2][self.histogram[1]>0], p0=p0, bounds=self.bounds)
    except:
      self.status = 'Failed'
      self.params = p0
      self.cov = [[0]*len(p0)]*len(p0)
      return False
    self.status = 'Successful'
    self.params = result[0]
    self.cov = result[1]
    return True
  
  def generate_fit_summary(self):
    fit_sum = self.function_name + ' Fit ' + self.status
    for p,pname in zip([physics.MeasuredQuantity(self.params[i], self.cov[i][i]**0.5) for i in range(0,len(self.params))], param_names):
      fit_sum += "\n" + r'${} = {}$'.format(pname, p.to_string())
    return fit_sum

  def plot(self):
    fit_sum = self.generate_fit_summary()
    plt.plot(self.histogram[0], self.function(self.histogram[0],*(self.params)),label=fit_sum)

class inv_mass_fitter(fitter):
  def __init__(self, sig_function_name, sig_function, bkg_function_name, bkg_function):
    self.sig_function = sig_function
    self.sig_function_name = sig_function_name
    self.bkg_function = bkg_function
    self.bkg_function_name = bkg_function_name

    sig_argspec = inspect.getargspec(sig_function)
    self.n_sig_arg = len(sig_argspec.args) - 1
    if bkg_function is not None:
      bkg_argspec = inspect.getargspec(bkg_function)
      self.n_bkg_arg = len(bkg_argspec.args) - 1
      function_name = self.sig_function_name + ' + ' + self.bkg_function_name
    else:
      function_name = self.sig_function_name

    fitter.__init__(self, function_name, self.get_signal_plus_background())
  
  def plot_bkg(self):
    plt.plot(self.histogram[0], self.bkg_function(self.histogram[0],*(self.params[self.n_sig_arg:])),label='Combinatorial Background',ls='--')

  def get_signal_plus_background(self):
    def signal_plus_background(x, *arg):
      return self.sig_function(x, *(arg[:self.n_sig_arg])) + self.bkg_function(x, *(arg[self.n_sig_arg:]))
    if self.bkg_function is not None:
      return signal_plus_background
    else:
      return self.sig_function

class lambda_fitter(inv_mass_fitter):
  @staticmethod
  def double_gaus_norm_fixed(nsigma=3.5):
    def _double_gaus_norm_fixed(x, N, mu, sigma, f, mu2):
      return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu2, sigma*nsigma)
    return _double_gaus_norm_fixed

  @staticmethod
  def double_gaus_norm_fixed_v2(nsigma=2.):
    def _double_gaus_norm_fixed_v2(x, N, mu, sigma, f, sigma2):
      return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+nsigma*sigma, sigma2)
    return _double_gaus_norm_fixed_v2

  @staticmethod
  def double_gaus_norm_fixed_v3(f=0.75):
    def _double_gaus_norm_fixed_v3(x, N, mu, sigma, sigma2):
      return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+2.*sigma, sigma2)
    return _double_gaus_norm_fixed_v3

  @staticmethod
  def double_gaus_norm_fixed_v4(nsigma=3., meanshift=2.):
    def _double_gaus_norm_fixed_v4(x, N, mu, sigma, f):
      return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+meanshift*sigma, sigma*nsigma)
    return _double_gaus_norm_fixed_v4

  @staticmethod
  def sqrt():
    def _sqrt(x, x0):
      return np.where(x>=x0, np.sqrt(x-x0), 0.)
    return _sqrt

  @staticmethod
  def expo():
    def _expo(x, x0, a):
      return np.where(x>=x0, np.exp(a*(x0-x)), 0.)
    return _expo

  @staticmethod
  def lambda_bkg(x, x0, a, C):
    return C * lambda_fitter.sqrt()(x, x0) * lambda_fitter.expo()(x, x0, a)

  def __init__(self, isMC):
    if isMC:
      inv_mass_fitter.__init__(self, 'Double Gaussian', self.double_gaus_norm_fixed_v4(), None, None)
    else:
      inv_mass_fitter.__init__(self, 'Double Gaussian', self.double_gaus_norm_fixed_v4(), '$C\sqrt{M - m_0}e^{-a(M-m_0)}$', self.lambda_bkg)
    self.isMC = isMC
    self.generate_fit_summary = self.lambda_generate_fit_summary

  def lambda_generate_fit_summary(self):
    fit_sum = self.function_name + ' Fit ' + self.status
    binw = self.bins[1] - self.bins[0]
    # Signal
    N = physics.MeasuredQuantity(self.params[0], self.cov[0][0]**0.5) * (1. / binw)
    mu = physics.MeasuredQuantity(self.params[1], self.cov[1][1]**0.5, "\mathrm{MeV}/c^2")
    sigma = physics.MeasuredQuantity(self.params[2], self.cov[2][2]**0.5, "\mathrm{MeV}/c^2")
    f = physics.MeasuredQuantity(self.params[3], self.cov[3][3]**0.5)
    mu2 = mu + sigma * 2.
    sigma2 = sigma * 3.
    fit_sum += "\n" + r'$S = {}$'.format(N.to_string())
    fit_sum += "\n" + r'$\mu_1 = {}$'.format(mu.to_string())
    fit_sum += "\n" + r'$\sigma_1 = {}$'.format(sigma.to_string())
    fit_sum += "\n" + r'$\mu_2 = \mu_1 + 2\sigma_1 = {}$'.format(mu2.to_string())
    fit_sum += "\n" + r'$\sigma_2 = 3\sigma_1 = {}$'.format(sigma2.to_string())
    fit_sum += "\n" + r'$f = {}$'.format(f.to_string())
    # Background
    if not self.isMC:
      m0 = physics.MeasuredQuantity(self.params[4], self.cov[4][4]**0.5, "\mathrm{MeV}/c^2")
      a = physics.MeasuredQuantity(self.params[5], self.cov[5][5]**0.5)
      C = physics.MeasuredQuantity(self.params[6], self.cov[6][6]**0.5)
      #min_mass = mu.value - sigma.value * 3.
      #max_mass = f.value * (mu2.value + sigma.value * 3.) - (f.value-1.) * (mu.value + sigma.value * 3.)
      #def integral_of_sqrt(m, m0):
      #  return 2. / 3. * (m - m0) ** 1.5
      #Nbkg = (integral_of_sqrt(max_mass, m0.value) - integral_of_sqrt(min_mass, m0.value)) * (C.value * binw)
      #significance = N.value / math.sqrt(Nbkg + N.value)
      fit_sum += "\n" + r'$m_0 = {}$'.format(m0.to_string())
      fit_sum += "\n" + r'$a = {}$'.format(a.to_string())
      fit_sum += "\n" + r'$C = {}$'.format(C.to_string())
      #fit_sum += "\n" + r'$B = {:.3f}$ in [{:.0f},{:.0f}]'.format(Nbkg, min_mass, max_mass)
      #fit_sum += "\n" + r'$\frac{{S}}{{\sqrt{{S+B}}}} = {:.3f}$ in [{:.0f},{:.0f}]'.format(significance, min_mass, max_mass)
    return fit_sum

class lambdab_fitter(inv_mass_fitter):
  @staticmethod
  def double_gaus_norm_single_mean_fixed_sigma(nsigma=3.5):
    def _double_gaus_norm_single_mean_fixed_sigma(x, N, mu, sigma, f):
      return gaus_norm(x, N*f, mu, sigma) + gaus_norm(x, N*(1.-f), mu, sigma*nsigma)
    return _double_gaus_norm_single_mean_fixed_sigma

  def __init__(self, isMC):
    if isMC:
      inv_mass_fitter.__init__(self, 'Double Gaussian', self.double_gaus_norm_single_mean_fixed_sigma(), None, None)
    else:
      inv_mass_fitter.__init__(self, 'Double Gaussian', self.double_gaus_norm_single_mean_fixed_sigma(), 'Cubic', cubic)
    self.isMC = isMC
    self.generate_fit_summary = self.lambdab_generate_fit_summary

  def lambdab_generate_fit_summary(self):
    fit_sum = self.function_name + ' Fit ' + self.status
    binw = self.bins[1] - self.bins[0]
    # Signal
    N = physics.MeasuredQuantity(self.params[0], self.cov[0][0]**0.5) * (1. / binw)
    mu = physics.MeasuredQuantity(self.params[1], self.cov[1][1]**0.5, "\mathrm{MeV}/c^2")
    sigma = physics.MeasuredQuantity(self.params[2], self.cov[2][2]**0.5, "\mathrm{MeV}/c^2")
    f = physics.MeasuredQuantity(self.params[3], self.cov[3][3]**0.5)
    sigma2 = sigma * 3.5
    fit_sum += "\n" + r'$S = {}$'.format(N.to_string())
    fit_sum += "\n" + r'$\mu = {}$'.format(mu.to_string())
    fit_sum += "\n" + r'$\sigma_1 = {}$'.format(sigma.to_string())
    fit_sum += "\n" + r'$\sigma_2 = 3.5\sigma_1 = {}$'.format(sigma2.to_string())
    fit_sum += "\n" + r'$f = {}$'.format(f.to_string())
    # Background
    if not self.isMC:
      a = physics.MeasuredQuantity(self.params[4], self.cov[4][4]**0.5)
      b = physics.MeasuredQuantity(self.params[5], self.cov[5][5]**0.5)
      c = physics.MeasuredQuantity(self.params[6], self.cov[6][6]**0.5)
      d = physics.MeasuredQuantity(self.params[7], self.cov[7][7]**0.5)
      min_mass = mu - (sigma * f - sigma2 * (f-1.)) * 3.
      max_mass = mu + (sigma * f - sigma2 * (f-1.)) * 3.
      dimless_min_mass = min_mass * 1.
      dimless_min_mass.units = ""
      dimless_max_mass = max_mass * 1.
      dimless_max_mass.units = ""
      #Nbkg = ((a * (dimless_max_mass - dimless_min_mass) * 0.5 + b) * (dimless_max_mass - dimless_min_mass)) * (1. / binw)
      #significance = N.value / math.sqrt(Nbkg.value + N.value)
      fit_sum += "\n" + r'$a = {}$'.format(a.to_string())
      fit_sum += "\n" + r'$b = {}$'.format(b.to_string())
      fit_sum += "\n" + r'$c = {}$'.format(c.to_string())
      fit_sum += "\n" + r'$d = {}$'.format(d.to_string())
      #fit_sum += "\n" + r'$B = {}$ in [{:.0f},{:.0f}]'.format(Nbkg.to_string(), min_mass.value, max_mass.value)
      #fit_sum += "\n" + r'$\frac{{S}}{{\sqrt{{S+B}}}} = {:.3f}$ in [{:.0f},{:.0f}]'.format(significance, min_mass.value, max_mass.value)
    return fit_sum

def generate_histogram(data, bins, density=False):
    x = np.array(0.5*(bins[1:] + bins[:-1]))
    y, _ = np.histogram(data, bins=bins)
    yerr = np.sqrt(y)
    if density:
        y_pdf, _ = np.histogram(data, bins=bins, density=True)
        with np.errstate(divide='ignore'):
            yerr = yerr * np.where(y>0, (y_pdf / y), 0)
        y = y_pdf
    return x, y, yerr

def fit_single(histogram, function_name, function, p0, bounds=(-np.inf,np.inf), param_names=None, plot_function=False):
  myfitter = fitter(function_name, function)
  myfitter.histogram = histogram
  myfitter.bounds = bounds
  myfitter.param_names = param_names
  r = myfitter.fit(p0)
  if r and plot_function:
    myfitter.plot()
  return (myfitter.params,myfitter.cov), myfitter.status

def fit_multiple(histograms, function_name, function, p0, bounds=(-np.inf,np.inf), param_names=None):
    for h in histograms:
      fit_single(h, function_name, function, p0, bounds, param_names, True)
