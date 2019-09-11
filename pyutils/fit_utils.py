""" This module contains utilities for fitting distributions.
"""

import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import physics

reload(physics)

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
  try:
    result = curve_fit(function, histogram[0][histogram[1]>0], histogram[1][histogram[1]>0], sigma=histogram[2][histogram[1]>0], p0=p0, bounds=bounds)
    status = 'Successful'
  except:
    result = (p0,[[0]*len(p0)]*len(p0))
    status = 'Failed'
  if plot_function:
    if param_names is None:
      param_names = [ 'p{}'.format(n) for n in range(0,len(p0)) ]
    fit_sum = function_name + ' Fit ' + status
    for p,pname in zip([physics.MeasuredQuantity(result[0][i], result[1][i][i]**0.5) for i in range(0,len(p0))], param_names):
      fit_sum += "\n" + r'${} = {}$'.format(pname, p.to_string())
    plt.plot(histogram[0], function(histogram[0],*(result[0])),label=fit_sum)
  return result, status

def fit_gaus_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Gaussian', gaus, p0, bounds, ['C', r'\mu', r'\sigma'])

def fit_gaus_norm_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Gaussian', gaus_norm, p0, bounds, ['N', r'\mu', r'\sigma'])

def fit_double_gaus_norm_single_mean_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm_single_mean, p0, bounds, [r'\mu', 'N', r'\sigma_1', 'f', r'\sigma_2'])

def fit_double_gaus_norm_single_mean_plus_lin_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussianm + ax+b', double_gaus_norm_single_mean_plus_lin, p0, bounds, [r'\mu', 'N', r'\sigma_1', 'f', r'\sigma_2', 'a', 'b'])

def fit_double_gaus_norm_single_mean_fixed_sigma_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm_single_mean_fixed_sigma, p0, bounds, [r'\mu', 'N', r'\sigma', 'f'])

def fit_double_gaus_norm_single_mean_fixed_sigma_plus_lin_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussianm + ax+b', double_gaus_norm_single_mean_fixed_sigma_plus_lin, p0, bounds, [r'\mu', 'N', r'\sigma', 'f', 'a', 'b'])

def fit_double_gaus_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus, p0, bounds, ['C_1', r'\mu_1', r'\sigma_1', 'C_2', r'\mu_2', r'\sigma_2'])

def fit_double_gaus_norm_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm, p0, bounds, ['N_1', r'\mu_1', r'\sigma_1', 'N_2', r'\mu_2', r'\sigma_2'])

def fit_double_gaus_norm_fixed_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm_fixed, p0, bounds, ['N', r'\mu', r'\sigma', 'f', r'\mu_2'])

def fit_double_gaus_norm_fixed_v2_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm_fixed_v2, p0, bounds, ['N', r'\mu', r'\sigma', 'f', r'\sigma_2'])

def fit_double_gaus_norm_fixed_v3_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm_fixed_v3, p0, bounds, ['N', r'\mu', r'\sigma', r'\sigma_2'])

def fit_double_gaus_norm_fixed_v4_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian', double_gaus_norm_fixed_v4, p0, bounds, ['N', r'\mu', r'\sigma', 'f'])

def fit_double_gaus_norm_fixed_v4_plus_lin_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian + ax+b', double_gaus_norm_fixed_v4_plus_lin, p0, bounds, ['N', r'\mu', r'\sigma', 'f', 'a', 'b'])

def fit_double_gaus_norm_fixed_v4_plus_const_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Double Gaussian + C', double_gaus_norm_fixed_v4_plus_const, p0, bounds, ['N', r'\mu', r'\sigma', 'f', 'C'])

def fit_crystal_ball_multiple(histograms, p0, bounds=(-np.inf,np.inf)):
  fit_multiple(histograms, 'Crystal Ball', crystal_ball, p0, bounds, ['C', r'\mu', r'\sigma', r'\alpha', 'n'])

def fit_multiple(histograms, function_name, function, p0, bounds=(-np.inf,np.inf), param_names=None):
    for h in histograms:
      fit_single(h, function_name, function, p0, bounds, param_names, True)

# Function library
def gaus(x, C, mu, sigma):
  x = np.atleast_1d(x)
  return C * np.exp(-(x-mu)**2/(2*sigma**2))

def gaus_norm(x, N, mu, sigma):
  x = np.atleast_1d(x)
  return N * gaus(x, 1., mu, sigma) / sigma / np.sqrt(2 * math.pi)

def double_gaus_norm_single_mean(x, mu, N, sigma1, f, sigma2):
  return gaus_norm(x, N*f, mu, sigma1) + gaus_norm(x, N*(1.-f), mu, sigma2)

def double_gaus_norm_single_mean_fixed_sigma(x, mu, N, sigma, f):
  return gaus_norm(x, N*f, mu, sigma) + gaus_norm(x, N*(1.-f), mu, sigma*3.5)

def double_gaus_norm_single_mean_plus_lin(x, mu, N, sigma1, f, sigma2, a, b):
  return gaus_norm(x, N*f, mu, sigma1) + gaus_norm(x, N*(1.-f), mu, sigma2) + linear(x, a, b)

def double_gaus_norm_single_mean_fixed_sigma_plus_lin(x, mu, N, sigma, f, a, b):
  return gaus_norm(x, N*f, mu, sigma) + gaus_norm(x, N*(1.-f), mu, sigma*3.5) + linear(x, a, b)

def double_gaus_norm(x, N1, mu1, sigma1, N2, mu2, sigma2):
  return gaus_norm(x, N1, mu1, sigma1)+gaus_norm(x, N2, mu2, sigma2)

def double_gaus_norm_fixed(x, N, mu, sigma, f, mu2):
  return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu2, sigma*3.)

def double_gaus_norm_fixed_v2(x, N, mu, sigma, f, sigma2):
  return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+2.*sigma, sigma2)

def double_gaus_norm_fixed_v3(x, N, mu, sigma, sigma2):
  return gaus_norm(x, N*0.75, mu, sigma)+gaus_norm(x, N*0.25, mu+2.*sigma, sigma2)

def double_gaus_norm_fixed_v4(x, N, mu, sigma, f):
  return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+2.*sigma, sigma*3.)

def double_gaus_norm_fixed_v4_plus_lin(x, N, mu, sigma, f, a, b):
  return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+2.*sigma, sigma*3.)+linear(x, a, b)

def double_gaus_norm_fixed_v4_plus_const(x, N, mu, sigma, f, C):
  return gaus_norm(x, N*f, mu, sigma)+gaus_norm(x, N*(1-f), mu+2.*sigma, sigma*3.)+C

def double_gaus(x, C1, mu1, sigma1, C2, mu2, sigma2):
  return gaus(x, C1, mu1, sigma1)+gaus(x, C2, mu2, sigma2)

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
  xpower = 1.
  result = 0.
  x = np.atleast_1d(x)
  for c in coeff:
    result += c * xpower
    xpower *= x
  return result

def linear(x, a, b):
  return polynomial(x, [b,a])

def quadratic(x, a, b, c):
  return polynomial(x, [c,b,a])

def gaus_plus_lin(x, C, mu, sigma, a, b):
  return gaus(x, C, mu, sigma) + linear(x, a, b)

def gaus_plus_quad(x, C, mu, sigma, a, b, c):
  return gaus(x, C, mu, sigma) + quadratic(x, a, b, c)

def gaus_plus_poly(x, C, mu, sigma, poly_coeff):
  return gaus(x, C, mu, sigma) + polynomial(x, poly_coeff)
