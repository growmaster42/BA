import numpy as np
import time
from intial_values import *

import copy
# Code doesn't work for a single spin
# some prefactors and user inputs#####################################################################################



# 'ket entry' in this code refers to the name of the ket; a ket corresponding to a system with two
# spins 1/2 will look like this |10>
# this doesn't represent the actual entries of the vector but simply the states up and down - the actual entries will
# stay hidden/unknown - to calculate the matrix entries it is not necessary to understand the actual entries, we
# only need the corresponding states
# this function creates basis vectors in the computational basis, based on the input given############################
start_time = time.time()
def basvec():
    # generates an array with zeros, dim x and number_of_spins entries per ket
    ketname = np.zeros((dim, number_of_spins), dtype=int)
    # makes rev_ketname globally available
    global rev_ketname
    # for loop that creates the parameter a that will be divided by number_of_states using divmod
    for i in range(dim - 1, -1, -1):
        # e.g. 5 spins 1/2 will result in dim = 32 basis vectors
        a = dim - 1 - i
        # for loop that runs through every entry of the ket
        for j in range(number_of_spins - 1, -1, -1):
            # creating a tuple that corresponds to modulus, x is a tuple that which contains the result of the
            # division and the rest (a/number_of_states, rest)
            x = divmod(a, number_of_states)
            # this puts the second value of the tuple (the rest) in the previously created ketname,
            ketname[i, j] = x[1]
            # this overwrites the previously generated a with the result of the division in divmod
            a = int(x[0])
            # ketname is being flipped for no reason just to have it nicely arranged
    rev_ketname = ketname[::-1]
    # this returns the ascending basis vectors as an array and terminates the function
    return rev_ketname #gibt alle Basisvektoren aufsteigend aus und beendet die Funktion

# scalar product that simply compares if two kets are the same##########################################################
def scalarproduct(ket1, ket2):
    ket1 = np.asarray(ket1)
    ket2 = np.asarray(ket2)
    if np.array_equal(ket1, ket2):
        # returns 1 if compared kets are equal
        return 1
    else:
        # returns 0 if compared kets are unequal
        return 0
#function that creates an array with all sz_i terms for n spins, takes spin, any set of
#basis vectors and a starting index i which marks the first state in the basis vector
#that is being referred to
#we create only an array with the length dim since we only need to do calculations
#for the diagonal elements of the matrix

#function that returns values for sz with a basis ket
def sz_dipol(vec, k, i, x):
    baseket_1 = copy_array(vec)
    if np.array_equal(baseket_1[k], baseket_1[x]):
        m = baseket_1[k][i] - spin
    else:
        m = 0
    return m
#function that returns the sum of sz*sz
def sz_sz_dipol(vec, k, l, x, p):
    m = sz_dipol(vec, k, l, x) * sz_dipol(vec, k, p, x)
    return m

#creates a diagnoal matrix out of the list
#matrix = np.diag(sum_sz_i_sz_j(basvec(), number_of_spins))
#matrix[matrix == -0] = 0
#print("Sum S_Z:\n", matrix)
#make a deep copy of any given array, necessary for altering vectors:
# WE DO NOT TOUCH OUR BASIS VECTORS
def copy_array(arr):
    # Erstelle eine tiefe Kopie des Arrays
    copied_arr = np.copy(arr)
    return copied_arr
#s_plus function takes arguments: any set of basis vectors, k refers to
#the k-th basis vector, l refers to the entry of the basis vector
#s again refers to the spin

def s_plus(vec, l, s, x):
    baseket_1 = copy_array(vec)
    baseket_1_altered = copy_array(vec)
    sqrt_prefactor_plus = []  # create empty array
    if baseket_1_altered[x][l % number_of_spins] != 2 * s: #maximal value a_max = s + m_max where
        # m_max = s is maximum -> a = 2 * s
        baseket_1_altered[x][l % number_of_spins] += 1
        m = (baseket_1[x][l % number_of_spins] - s)
        sqrt = np.sqrt(s * (s + 1) - m * (m + 1)) * h

        sqrt_prefactor_plus.append(sqrt)  # Ergebnis mit append hinzufügen
    else:
        m = 0

        sqrt_prefactor_plus.append(m)
        # Ergebnis mit append hinzufügen
    return sqrt_prefactor_plus, baseket_1_altered[x]

#print("Basevec: ", basvec()[k], ",Stelle: ", l, "\n", "s_plus_function:\n", s_plus(basvec(), k, l, spin))

def s_minus(vec, l, s, x):
    baseket_1 = copy_array(vec)
    baseket_2_altered = copy_array(vec)
    sqrt_prefactor_minus = []  # create empty array
    if baseket_2_altered[x][l % number_of_spins] != 0: #minimal value a_min = s - m_max where
        # m_min = s is minimum
        baseket_2_altered[x][l % number_of_spins] -= 1
        m = (baseket_1[x][l % number_of_spins] - s)
        sqrt = np.sqrt(s * (s + 1) - m * (m - 1)) * h
        sqrt_prefactor_minus.append(sqrt)  # Ergebnis mit append hinzufügen
    else:
        m = 0 #(baseket_1[x][l % number_of_spins] - s)
        sqrt = 0# np.sqrt(s * (s + 1) - m * (m - 1)) * h
        sqrt_prefactor_minus.append(sqrt)  # Ergebnis mit append hinzufügen

    return sqrt_prefactor_minus, baseket_2_altered[x]

#print("Basevec: ", basvec()[k], ",Stelle: ", l, "\n", "s_minus_function:\n", s_minus(basvec(), k, l, spin, x))

def s_minus_s_plus(vec, k, l, x, p):
    sqrt_prefactor_plus, baseket_1_altered = s_plus(vec, p, spin, x)
    baseket_1_altered_1_empty = np.zeros((dim, number_of_spins), dtype=int)
    baseket_1_altered_1_empty[x] = np.array(baseket_1_altered)
    sqrt_prefactor_minus, baseket_2_altered = s_minus(baseket_1_altered_1_empty,  l, spin, x)
    #print(baseket_1_altered_1)
    if np.array_equal(baseket_2_altered, vec[k]):
        return sqrt_prefactor_minus[0]*sqrt_prefactor_plus[0]
    else:
        return 0



#function that returns the product of s+ s-
def s_plus_s_minus(vec, k, l, x, p):
    sqrt_prefactor_plus, baseket_1_altered = s_plus(vec,  l, spin, x)
    baseket_1_altered_1_empty = np.zeros((dim, number_of_spins), dtype=int)
    baseket_1_altered_1_empty[x] = np.array(baseket_1_altered)
    sqrt_prefactor_minus, baseket_2_altered = s_minus(baseket_1_altered_1_empty,  p, spin, x)
    if np.array_equal(baseket_2_altered, vec[k]):
        return sqrt_prefactor_minus[0]*sqrt_prefactor_plus[0]
    else:
        return 0

def s_plus_s_plus(vec, k, l, x, p):
    sqrt_prefactor_plus_i, baseket_plus_altered_i = s_plus(vec,  l, spin, x)
    baseket_1_altered_1_empty = np.zeros((dim, number_of_spins), dtype=int)
    baseket_1_altered_1_empty[x] = np.array(baseket_plus_altered_i)
    sqrt_prefactor_plus_j, baseket_plus_altered_j = s_plus(baseket_1_altered_1_empty,  p, spin, x)
    if np.array_equal(baseket_plus_altered_j, vec[k]):
        return sqrt_prefactor_plus_i[0] * sqrt_prefactor_plus_j[0]
    else:
        return 0

def s_minus_s_minus(vec, k, l, x, p):
    sqrt_prefactor_minus_i, baseket_minus_altered_i = s_minus(vec,  l, spin, x)
    baseket_1_altered_1_empty = np.zeros((dim, number_of_spins), dtype=complex)
    baseket_1_altered_1_empty[x] = np.array(baseket_minus_altered_i)
    sqrt_prefactor_minus_j, baseket_minus_altered_j = s_minus(baseket_1_altered_1_empty,  p, spin, x)
    if np.array_equal(baseket_minus_altered_j, vec[k]):
        return sqrt_prefactor_minus_i[0] * sqrt_prefactor_minus_j[0]
    else:
        return 0
#function to create a single matrix entry from summands, where k refers to row and x to column
def summand_matrix_entry(operator, vec, k, x):
    sum = 0
    for l in range(number_of_spins - 1):
        for p in range(number_of_spins):
            sum += operator(vec, k, l, x, p)
    return sum

import math

def count_couplings():
    return math.comb(number_of_spins, 2)


