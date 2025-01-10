import numpy as np
from intial_values import *


# this function returns a list with all m_i for each ket
def s_z(vec, spin, k, i):
    basis_vectors = vec.copy()
    vector = basis_vectors[k]
    m_i = vector[i] - spin
    return m_i


# this function returns a list with all pre-factors that s_plus produces for each ket
def s_plus(vec, spin, k, i):
    basis_vectors = vec.copy()
    vector = basis_vectors[k]
    m = vector[i] - spin
    sqrt = np.sqrt(spin * (spin + 1) - m * (m + 1)) * h

    return sqrt


def s_minus(vec, spin, k, i):
    basis_vectors = vec.copy()
    vector = basis_vectors[k]
    m = vector[i] - spin
    sqrt = np.sqrt(spin * (spin + 1) - m * (m - 1)) * h
    return sqrt


# this function changes the basis vector given the way s_plus * s_plus interact with the ket, i should max. be
# num_spins , since i refers to the first spin
def s_plus_s_plus_ket_change(vec, k, i, j):
    basis_vector = vec[k].copy()
    basis_vector[i] += 1
    basis_vector[j] += 1
    return basis_vector


# this function changes the basis vector given the way s_plus * s_plus interact with the ket, i should max. be
# num_spins , since i refers to the first spin
def s_plus_s_minus_ket_change(vec, k, i, j):
    basis_vector = vec[k].copy()
    basis_vector[i] += 1
    basis_vector[j] -= 1
    return basis_vector


# this function changes the basis vector given the way s_plus * s_plus interact with the ket, i should max. be
# num_spins , since i refers to the first spin
def s_minus_s_plus_ket_change(vec, k, i, j):
    basis_vector = vec[k].copy()
    basis_vector[i] -= 1
    basis_vector[j] += 1
    return basis_vector


# this function changes the basis vector given the way s_plus * s_plus interact with the ket, i should max. be
# num_spins , since i refers to the first spin
def s_minus_s_minus_ket_change(vec, k, i, j):
    basis_vector = vec[k].copy()
    basis_vector[i] -= 1
    basis_vector[j] -= 1
    return basis_vector




