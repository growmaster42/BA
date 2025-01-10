import numpy as np
import copy
import time
# Code doesn't work for a single spin
# some prefactors and user inputs#####################################################################################
from intial_values import *



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
def s_i_z_s_j_z(vec,s):
    # initialising an array with only zeros and size dim x dim,
    # allows complex numbers x+yj (where j is imaginary unit)
    sz_sz = np.zeros((dim, dim), dtype=complex)
    # calling basvec() function and saving the created basis vectors
    # in a variable named basket_1, advantage: basvec()
    # will be called once in the whole function which is more efficient
    baseket_1 = vec
    # for loop in for loop each with range dim in order to go from
    # the top left corner of the matrix systematically
    for i in range(dim):
        for j in range(dim):
            # initialising m = 0 in order to make sure m is defined
            m = 0
            # this for loop iterates through a single ket and analyses
            # every entry from left to right
            for k in range(number_of_spins):
                # sz_ij is the prefactor that comes out when we let  s^z_i*s^z_j affect the ket,
                # (found in Script J.S. page 8); in python % = modulus and x % y returns the rest
                # of the division; we use it here in order to go through every entry of the ket
                # until we reach the number_of_spins which will then make the modulus operation
                # return 0 which will refer to the first entry of the ket
                # in the second expression that is multiplied we start at k - 1 which refers to
                # the next state in the ket left, so if the ket looks like |10> the first entry
                # is being referred to with the right expression in braces, the second expression
                # is being referred to with the left expression in braces
                sz_i = baseket_1[i][k % number_of_spins] - s
                sz_j = baseket_1[j][k - 1] - s
                sz_ij = sz_i * sz_j  * scalarproduct(baseket_1[i], baseket_1[j])

            if number_of_spins == 2:
                sz_sz[i][j] +=  sz_ij

            else:
                sz_sz[i][j] += sz_ij

    return sz_sz

def s_i_plus_s_j_minus(vec, s):
    baseket_1 = basvec()
    s_plus_s_minus = np.zeros((dim,dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            # initialising m = 0 in order to make sure m is defined
            m = 0
            # this for loop iterates through a single ket and analyses
            # every entry from left to right
            for k in range(number_of_spins):
                ket_new2 = vec.copy()
                # setting n
                n = number_of_states - 1
                if ket_new2[i][k % number_of_spins] != 0 and ket_new2[i][k - 1] != n: #vorderer term spm
                    ket_new2[i][k % number_of_spins] -= 1
                    ket_new2[i][k - 1] += 1
                    m_jm2 = (baseket_1[i][k] - s)
                    m_jp2 = (baseket_1[i][k - 1] - s)
                    sqrtm2 = np.sqrt(s * (s + 1) - m_jm2 * (m_jm2 - 1)) * h         # (1)
                    sqrtp2 = np.sqrt(s * (s + 1) - m_jp2 * (m_jp2 + 1)) * h         # (1)
                    scm = (sqrtm2) * (sqrtp2) * scalarproduct(ket_new2[i], vec[j])

                else:
                    scm = 0
                m += scm
                s_plus_s_minus[i][j] = m
    return s_plus_s_minus

def s_i_minus_s_j_plus(vec, s):
    s_minus_s_plus = np.zeros((dim, dim), dtype=complex)
    baseket_1 = vec
    ket_new = vec.copy()
    n = number_of_states - 1
    for i in range(dim):
        for j in range(dim):
            # initialising m = 0 in order to make sure m is defined
            m = 0
            # this for loop iterates through a single ket and analyses
            # every entry from left to right
            for k in range(number_of_spins):
                if ket_new[i][k % number_of_spins] != n and ket_new[i][k - 1] != 0:  # hinterer Term smp
                    ket_new[i][k % number_of_spins] += 1
                    ket_new[i][k - 1] -= 1

                    m_jp1 = (baseket_1[i][
                                 k % number_of_spins] - s)  # (1): Berechnung des Wertes den der Operator aus dem Ket zieht bzw.
                    # Rücktransformation nach "a = s+m", also Berechnung von m = a-s
                    m_jm1 = (baseket_1[i][k - 1] - s)  # (1)

                    sqrtp1 = np.sqrt(s * (s + 1) - m_jp1 * (m_jp1 + 1)) * h
                    sqrtm1 = np.sqrt(s * (s + 1) - m_jm1 * (m_jm1 - 1)) * h

                    scp = (sqrtp1) * (sqrtm1) * scalarproduct(ket_new[i], vec[j])
                else:
                    scp = 0
                m += scp
                s_minus_s_plus[i][j] = m
        return s_minus_s_plus

# Heisenberg: Creates Hamiltonian that corresponds to the Heisenberg-Model####################################
def Hamiltonian(vec, s):
    spm = s_i_plus_s_j_minus(vec, s)
    smp = s_i_minus_s_j_plus(vec, s)
    szz = s_i_z_s_j_z(vec, s)
    H = -2 * (szz + 0.5 * (spm + smp))
    return(H)

    return H
#this function creates a matrix that represents the Heisenberg-Interaction. It uses the "Hamiltonian" Function.
def Heisenberg():
    H_H = J_ij * Hamiltonian(basvec(), spin)
    return H_H




