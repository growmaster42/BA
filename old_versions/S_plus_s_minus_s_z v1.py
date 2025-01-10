import numpy as np
import copy
# Code doesn't work for a single spin
# some prefactors and user inputs#####################################################################################
# hbar set to 1
h = 1

# coupling constant J_ij set to 1
J_ij = 1
# mu_0(gmu_b) ** 2 / (4 * pi * Angstrom ** 3 ) = prefactor_sum_dd: prefactor of the dipole-dipole interaction,
# value from  J.S. script
prefactor_sum_dd = 2.4917789
# Input function that asks the number of spins that are being examined
number_of_spins = 3# int(input("Type Number of Spins: "))
# Input function that asks for the type of spin: eg. Spin 1/2, spin 1 etc.
spin = 0.5 #float(input("Enter Spin: "))
# This function calculates the number of states a certain spin can be in, for spin 1/2 it would be 2 states
number_of_states = int(2 * spin + 1)
# this function calculates the dimension of the corresponding hilbert space, which is also the number of basis vectors
dim = number_of_states ** number_of_spins
# function that asks for the radius of the ring
ringrad = 1#float(input("Enter Ring Radius: "))
print("DIM =", dim)

# 'ket entry' in this code refers to the name of the ket; a ket corresponding to a system with two
# spins 1/2 will look like this |10>
# this doesn't represent the actual entries of the vector but simply the states up and down - the actual entries will
# stay hidden/unknown - to calculate the matrix entries it is not necessary to understand the actual entries, we
# only need the corresponding states
# this function creates basis vectors in the computational basis, based on the input given############################
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
def sz(s, vec, i):
    baseket = vec
    diag_sz = []
    for n in range(dim):
        m = (baseket[n][i] - s)
        diag_sz.append(m)
    return np.array(diag_sz)
#function that creates the product of sz_i*sz_j using only sz_i and shifting the index by + 1
# which will then refer to the next neighbour
def sz_i_sz_j(vec ,i, n):
    A = sz(spin, vec, i % n) * sz(spin, vec, i + 1)
    return A
#this function sums up all coefficients that were calculated in sz_i_sz_j, takes any set
# of basis vectors and the number_of_states as arguments - i is no argument any longer since
#the loop fulfills this job
def sum_sz_i_sz_j(vec, n):
    sum = 0
    for i in range(n - 1):
        sum += sz_i_sz_j(vec, i, n)
    return sum

#creates a diagnoal matrix out of the list
print(np.diag(sum_sz_i_sz_j(basvec(), number_of_spins)))
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
k = 1
l = 0
x = 2
def s_plus(vec, k, l, s, x):
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
        m = (baseket_1[x][l % number_of_spins] - s)


        sqrt = np.sqrt(s * (s + 1) - m * (m + 1)) * h

        sqrt_prefactor_plus.append(sqrt)
        # Ergebnis mit append hinzufügen
    return sqrt_prefactor_plus, baseket_1_altered[x]

#print("Basevec: ", basvec()[k], ",Stelle: ", l, "\n", "s_plus_function:\n", s_plus(basvec(), k, l, spin))

def s_minus(vec, k, l, s, x):
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
        m = (baseket_1[x][l % number_of_spins] - s)
        sqrt = np.sqrt(s * (s + 1) - m * (m - 1)) * h
        sqrt_prefactor_minus.append(sqrt)  # Ergebnis mit append hinzufügen

    return sqrt_prefactor_minus, baseket_2_altered[x]

#print("Basevec: ", basvec()[k], ",Stelle: ", l, "\n", "s_minus_function:\n", s_minus(basvec(), k, l, spin, x))

def s_minus_plus(vec, k, l,  x):
    Matrix = np.zeros((dim, dim), dtype=int)
    sqrt_prefactor_plus, baseket_1_altered = s_plus(vec, k, l + 1, spin, x)
    baseket_1_altered_1_empty = np.zeros((dim, number_of_spins), dtype=int)
    baseket_1_altered_1_empty[x] = np.array(baseket_1_altered)
    sqrt_prefactor_minus, baseket_2_altered = s_minus(baseket_1_altered_1_empty, k, l, spin, x)
    #print(baseket_1_altered_1)
    if np.array_equal(baseket_2_altered, vec[k]):
        Matrix[k][x] += sqrt_prefactor_minus[0]*sqrt_prefactor_plus[0]
    else:
        Matrix[k][x] += 0
    return Matrix


def s_plus_s_minus(vec, k, l,  x):
    Matrix = np.zeros((dim, dim), dtype=int)
    sqrt_prefactor_plus, baseket_1_altered = s_plus(vec, k, l, spin, x)
    baseket_1_altered_1_empty = np.zeros((dim, number_of_spins), dtype=int)
    baseket_1_altered_1_empty[x] = np.array(baseket_1_altered)
    sqrt_prefactor_minus, baseket_2_altered = s_minus(baseket_1_altered_1_empty, k, l + 1, spin, x)
    if np.array_equal(baseket_2_altered, vec[k]):
        Matrix[k][x] += sqrt_prefactor_minus[0]*sqrt_prefactor_plus[0]
    else:
        Matrix[k][x] += 0
    return Matrix
Matrix = s_minus_plus(basvec(), k, l, x) + s_plus_s_minus(basvec(), k, l, x)
print(Matrix)
np.savetxt('matrix.txt', Matrix, fmt='%d')
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
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
    szz = sum_sz_i_sz_j(vec, s)
    H = -2 * (szz + 0.5 * (spm + smp))
    return(H)

#this function creates a matrix that represents the Heisenberg-Interaction. It uses the "Hamiltonian" Function.
def Heisenberg():
    H_H = J_ij * Hamiltonian(basvec(), spin)
    return H_H
Heis = Heisenberg()
print("Heisenberg alt", Heis)
#this function creates the 3x3 matrices that is being used in the dipole-dipole interaction. It derives from two normal vectors that define orientation and position of the spins
def orientation():
    #initializing list for each matrix that occurs in the sum
    orientation_matrix_list = []
    #creating the angles that occur inbetween the spins by defining alpha_i and alpha_j
    for i in range(number_of_spins + 1):
        alpha_i = 2 * np.pi * (i + 1) / number_of_spins
        for j in range(i+2, number_of_spins + 1):
            alpha_j = 2 * np.pi * j / number_of_spins
            #abbreviate differences between cosines and sines of angles
            delta_cos = np.cos(alpha_j) - np.cos(alpha_i)
            delta_sin = np.sin(alpha_j) - np.sin(alpha_i)
            #defining the 3x3 matrix
            orientation_matrix = np.array([[delta_cos ** 2,            delta_cos * delta_sin, 0],
                               [delta_sin * delta_cos,     delta_sin ** 2,        0],
                               [0,                  0,                            0]])
            orientation_matrix_list.append(orientation_matrix)

    return orientation_matrix_list


# this function calculates the angles between all spins and all interactions
def spin_angles():
    spin_angles_list = []
    fraction = float(1 / number_of_spins)
    for j in range(number_of_spins):
        for i in range(1, number_of_spins - j):
            spin_angle = 2 * np.pi * fraction * i
            spin_angles_list.append(spin_angle)
    return spin_angles_list


# calling the spin_angles function, this function calculates the actual distance between all spins
def spin_distance():
    r_ij_list = []
    for i in spin_angles():
        r_ij = 2 * ringrad * np.sin(i / 2 )
        r_ij_list.append(r_ij)
    return r_ij_list
# calling the spin_distance function, this function creates the r_ij dependent prefactor that varies with
# each part of the sum
def prefactor_sum():
    prefactor_sum_list = []
    for i in spin_distance():
        R_ij = 1 / ((i / 10 ** (-10)) ** 3)
        prefactor_sum_list.append(R_ij)

    return prefactor_sum_list










