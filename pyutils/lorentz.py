""" Lorentz transformations
"""

import numpy as np

def generate_decay(M, m1, m2):
    p = np.sqrt((M**2 + m1**2 - m2**2)**2 - 4*M**2*m1**2) / 2 / M
    return np.array([p, 0, 0, np.sqrt(m1**2+p**2)]), np.array([-p, 0, 0, np.sqrt(m2**2+p**2)])

def mass(P):
    return np.sqrt(P[3]**2 - P[0]**2 - P[1]**2 - P[2]**2)

def beta(P):
    p = np.sqrt(P[0]**2 + P[1]**2 + P[2]**2)
    return p / P[3]

def gamma(P):
    return 1./np.sqrt(1-beta(P)**2)

def boost(orig, n, b):
    g = 1./np.sqrt(1-b**2)
    ppar = orig[0] * n[0] + orig[1] * n[1] + orig[2] * n[2]
    E_star = orig[3] * g - ppar * g * b
    ppar_star = g*ppar - g*b*orig[3]
    p_star = []
    for n_proj, orig_proj in zip(n, orig[:-1]):
        p_new = orig_proj + (ppar_star - ppar) * n_proj
        p_star.append(p_new)
    return np.array([p_star[0], p_star[1], p_star[2], E_star])

def boost_from_mom(orig, boost_vect):
    return boost(orig, unit_vector(boost_vect), beta(boost_vect))

def four_vector_from_PXPYPZM(px, py, pz, M):
    E = np.sqrt(M**2 + px**2 + py**2 + pz**2)
    return np.array([px, py, pz, E])

def unit_vector(v):
    mag = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return (v[0] / mag, v[1] / mag, v[2] / mag)

def cos_angle(v1, v2):
    scalar_product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    norm = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) * np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    return scalar_product / norm
