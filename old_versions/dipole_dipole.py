from intial_values import *
from S_plus_s_minus_s_z import *
import numpy as np
quarter = 0.25
quarter_i = - 0.25j
# function to create the orientation matrix that origins from the connection vectors between two spins
#returns 3x3 arrays within array with length of all possible interactions
def orientation():
    all_B_ij = []
    all_r_ij = []
    #initializing list for each matrix that occurs in the sum
    orientation_matrix_list = []
    all_r_ij = summand_prefactor(number_of_spins)
    #creating the angles that occur inbetween the spins by defining alpha_i and alpha_j
    for i in range(number_of_spins + 1):
        B_ij = []  # Liste für die aktuelle äußere Iteration
        r_ij = []
        alpha_i = 2 * np.pi * (i + 1) / number_of_spins
        for j in range(i + 2, number_of_spins + 1):
            alpha_j = 2 * np.pi * j / number_of_spins
            spin_angle = alpha_j - alpha_i
            r = 2 * ringrad * np.sin(spin_angle / 2)
            B = 1 / (r ** 3)
            frac = B * (ringrad / r) ** 2
            #abbreviate differences between cosines and sines of angles
            delta_cos = np.cos(alpha_j) - np.cos(alpha_i)
            delta_sin = np.sin(alpha_j) - np.sin(alpha_i)
            #defining the 3x3 matrix
            orientation_matrix = np.array([[frac * delta_cos ** 2,            frac * delta_cos * delta_sin, 0],
                               [frac * delta_sin * delta_cos,     frac * delta_sin ** 2,        0],
                               [0,                  0,                            0]])
            orientation_matrix_list.append(orientation_matrix)
            r_ij.append(r)
            B_ij.append(B)  # Winkel zur aktuellen Liste hinzufügen

    return orientation_matrix_list
#function that calculates the prefactors B_ij and r_ij, returns list
def summand_prefactor(n):
    all_B_ij = []
    all_r_ij = []
    # Liste für alle inneren Listen
    fraction = 1 / n  # Bruch für die Berechnung der Winkel

    for j in range(n):
        B_ij = []  # Liste für die aktuelle äußere Iteration
        r_ij = []
        for i in range(1, n - j):
            spin_angle = 2 * np.pi * fraction * i
            r = (2 * ringrad * np.sin(spin_angle / 2))
            r_squared = (2 * ringrad * np.sin(spin_angle / 2)) ** 2
            B = 1 / (r ** 3)
            r_ij.append(r_squared)
            B_ij.append(B)  # Winkel zur aktuellen Liste hinzufügen
        all_B_ij.append(B_ij)
        all_r_ij.append(r_ij)# Aktuelle Liste zu allen Listen hinzufügen

    return all_B_ij, all_r_ij
#the functions dipol_1, dipol_2, dipol_3 and dipol_4 create various terms for the whole matrix
#k refers to row, x to column
def dipol_1():
    vec = basvec()
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    r_ij = all_r_ij[0]
    orient = orientation()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            prefactor = orient[0][0, 0]
            summand_1 = quarter * summand_matrix_entry_distance(s_plus_s_plus, vec, k, x)
            summand_2 = quarter * summand_matrix_entry_distance(s_plus_s_minus, vec, k, x)
            summand_3 = quarter * summand_matrix_entry_distance(s_minus_s_plus, vec, k, x)
            summand_4 = quarter * summand_matrix_entry_distance(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += prefactor * (summand_1 + summand_2 + summand_3 + summand_4)
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def dipol_2():
    vec = basvec()
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    r_ij = all_r_ij[0]
    orient = orientation()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            prefactor = orient[0][0, 0]
            summand_1 = quarter_i * summand_matrix_entry_distance(s_plus_s_plus, vec, k, x)
            summand_2 = - quarter_i * summand_matrix_entry_distance(s_plus_s_minus, vec, k, x)
            summand_3 = quarter_i * summand_matrix_entry_distance(s_minus_s_plus, vec, k, x)
            summand_4 = - quarter_i * summand_matrix_entry_distance(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += prefactor * (summand_1 + summand_2 + summand_3 + summand_4)
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def dipol_3():
    vec = basvec()
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    r_ij = all_r_ij[0]
    orient = orientation()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            prefactor = orient[0][0, 0]
            summand_1 = quarter_i * summand_matrix_entry_distance(s_plus_s_plus, vec, k, x)
            summand_2 = quarter_i * summand_matrix_entry_distance(s_plus_s_minus, vec, k, x)
            summand_3 = - quarter_i * summand_matrix_entry_distance(s_minus_s_plus, vec, k, x)
            summand_4 = - quarter_i * summand_matrix_entry_distance(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += prefactor * (summand_1 + summand_2 + summand_3 + summand_4)
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def dipol_4():
    vec = basvec()
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    r_ij = all_r_ij[0]
    orient = orientation()
    matrix = np.zeros((dim, dim), dtype=complex)

    for k in range(dim):
        for x in range(dim):
            prefactor = orient[0][1, 1]
            summand_1 = - quarter * summand_matrix_entry_distance(s_plus_s_plus, vec, k, x)
            summand_2 = quarter * summand_matrix_entry_distance(s_plus_s_minus, vec, k, x)
            summand_3 = quarter * summand_matrix_entry_distance(s_minus_s_plus, vec, k, x)
            summand_4 = -  quarter * summand_matrix_entry_distance(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += prefactor * (summand_1 + summand_2 + summand_3 + summand_4)

            # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def heisenberg_part():
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = summand_matrix_entry_distance_heis(s_plus_s_minus, vec, k, x) #s_i_plus * s_j_minus
            summand_2 = summand_matrix_entry_distance_heis(s_minus_s_plus, vec, k, x)    #s_i_minus * s_j_plus
            summand_3 = summand_matrix_entry_distance_heis(sz_sz_dipol, vec, k, x)   #s_i_z * s_j_z

            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += 0.5 * (summand_1 + summand_2) + summand_3  # Hier wird die Berechnung in die Matrix geschrieben
    heisenberg_matrix = matrix
    return heisenberg_matrix

#this function calls prefactors from the suitable list and combines them with an operator
def summand_matrix_entry_distance(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0

    for l in range(number_of_spins - 1):
        prefactor_list = all_B_ij[l]
        for p in range(1, number_of_spins):
            sum += operator(vec, k, l, x, p)

    return sum

def summand_matrix_entry_distance_heis(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0

    for l in range(number_of_spins - 1):
        prefactor_list = all_B_ij[l]
        prefactor_list_2 = all_r_ij[l]
        for p in range(1, number_of_spins):
            sum += prefactor_list[0] * operator(vec, k, l, x, p)

    return sum



def nhdw(): #nicht heisenbergsche Dipol-Wirkung, also der hintere Teil von dem dikken Hässlon
    matrix = - 3 * (dipol_1() + dipol_2() + dipol_3() + dipol_4())
    return matrix


def dipol_dipol():
    H_dd = prefactor_sum_dd * (nhdw() + heisenberg_part())
    return H_dd

