from intial_values import *
from S_plus_s_minus_s_z import *
import numpy as np
quarter = 0.25
quarter_i = - 0.25j
def orientation():
    #initializing list for each matrix that occurs in the sum
    orientation_matrix_list = []
    all_r_ij = summand_prefactor(number_of_spins)
    #creating the angles that occur inbetween the spins by defining alpha_i and alpha_j
    for i in range(1, number_of_spins + 1):
        alpha_i = 2 * np.pi * (i) / number_of_spins
        for j in range(i + 1, number_of_spins + 1):
            alpha_j = 2 * np.pi * j / number_of_spins
            #abbreviate differences between cosines and sines of angles
            delta_cos = np.cos(alpha_j) - np.cos(alpha_i)
            delta_sin = np.sin(alpha_j) - np.sin(alpha_i)
            #defining the 3x3 matrix

            orientation_matrix = np.array([[delta_cos ** 2,            delta_cos * delta_sin, 0],
                               [delta_sin * delta_cos,     delta_sin ** 2,        0],
                               [0,                  0,                            0]], dtype=complex)
            orientation_matrix_list.append(orientation_matrix)

    return orientation_matrix_list

def summand_prefactor_1(n):
    all_b_ij = []
    all_r_ij = []
    # Liste für alle inneren Listen
    fraction = 1 / n  # Bruch für die Berechnung der Winkel
    for i in range(1, n + 1):
        b_ij = []  # list for j fix and i variable
        r_ij = []  # list for j fix and i variable
        for j in range(i + 1, n + 1):
            spin_angle = 2 * np.pi * fraction * (j - i)
            r = 2 * ringrad * np.sin(spin_angle / 2)
            b = 1 / (r ** 3)
            r_ij.append(r)
            b_ij.append(b)  # appending B_ij for fixed j and all i to list
        all_b_ij.append(b_ij)
        all_r_ij.append(r_ij)# Aktuelle Liste zu allen Listen hinzufügen

    return all_b_ij, all_r_ij
def summand_prefactor(n):
    all_b_ij = []
    all_r_ij = []
    fraction = 1 / n  # Bruch für die Berechnung der Winkel
    for i in range(1, n):  # i muss kleiner als n sein
        b_ij = []  # Liste für j fix und i variabel
        r_ij = []  # Liste für j fix und i variabel
        for j in range(i + 1, n + 1):
            spin_angle = 2 * np.pi * fraction * (j - i)
            r = 2 * ringrad * np.sin(spin_angle / 2)
            b = 1 / (r ** 3)
            r_ij.append(r)
            b_ij.append(b)
        all_b_ij.append(b_ij)
        all_r_ij.append(r_ij)

    return all_b_ij, all_r_ij

def dipol_1():
    vec = basvec()
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    r_ij = all_r_ij[0]
    orient = orientation()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):

            summand_1, orient_ij = summand_matrix_entry_distance_1(s_plus_s_plus, vec, k, x)
            summand_2, orient_ij = summand_matrix_entry_distance_1(s_plus_s_minus, vec, k, x)
            summand_3, orient_ij = summand_matrix_entry_distance_1(s_minus_s_plus, vec, k, x)
            summand_4, orient_ij = summand_matrix_entry_distance_1(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += quarter * orient_ij[0, 0] * (summand_1 + summand_2 + summand_3 + summand_4)
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def dipol_2():
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):

            summand_1, orient_ij = summand_matrix_entry_distance_2(s_plus_s_plus, vec, k, x)
            summand_2, orient_ij = summand_matrix_entry_distance_2(s_plus_s_minus, vec, k, x)
            summand_3, orient_ij = summand_matrix_entry_distance_2(s_minus_s_plus, vec, k, x)
            summand_4, orient_ij = summand_matrix_entry_distance_2(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += quarter_i * orient_ij[0, 1] * (summand_1 - summand_2 + summand_3 - summand_4)
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def dipol_3():
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1, orient_ij = summand_matrix_entry_distance_3(s_plus_s_plus, vec, k, x)
            summand_2, orient_ij = summand_matrix_entry_distance_3(s_plus_s_minus, vec, k, x)
            summand_3, orient_ij = summand_matrix_entry_distance_3(s_minus_s_plus, vec, k, x)
            summand_4, orient_ij = summand_matrix_entry_distance_3(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += quarter_i * orient_ij[1, 0] * (summand_1 + summand_2 - summand_3 - summand_4)
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def dipol_4():
    vec = basvec()
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1, orient_ij = summand_matrix_entry_distance_4(s_plus_s_plus, vec, k, x)
            summand_2, orient_ij = summand_matrix_entry_distance_4(s_plus_s_minus, vec, k, x)
            summand_3, orient_ij = summand_matrix_entry_distance_4(s_minus_s_plus, vec, k, x)
            summand_4, orient_ij = summand_matrix_entry_distance_4(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += quarter * orient_ij[1, 1] * (- summand_1 + summand_2 + summand_3 - summand_4)

            # Hier wird die Berechnung in die Matrix geschrieben
    return matrix




def summand_matrix_entry_distance_1(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0
    orient = orientation()
    for l in range(number_of_spins - 1):
        B_ij = all_B_ij[l]
        r_ij = all_r_ij[l]
        for p in range(1, number_of_spins):
            orient_ij = orient[p - 1]
            sum += B_ij[p - l - 1] * (ringrad /r_ij[p - l - 1]) ** 2 * operator(vec, k, l, x, p)
    return sum, orient_ij



def summand_matrix_entry_distance_2(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0
    orient = orientation()
    for l in range(number_of_spins - 1):
        B_ij = all_B_ij[l]
        r_ij = all_r_ij[l]
        for p in range(1, number_of_spins):
            orient_ij = orient[p - 1]
            sum += B_ij[p - l - 1] * (ringrad /r_ij[p - l - 1]) ** 2 * operator(vec, k, l, x, p)

    return sum, orient_ij



def summand_matrix_entry_distance_3(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0
    orient = orientation()
    for l in range(number_of_spins - 1):
        B_ij = all_B_ij[l]
        r_ij = all_r_ij[l]
        for p in range(1, number_of_spins):
            orient_ij = orient[p - 1]
            sum += B_ij[p - l - 1] * (ringrad /r_ij[p - l - 1]) ** 2 * operator(vec, k, l, x, p)
    return sum, orient_ij



def summand_matrix_entry_distance_4(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0
    orient = orientation()
    for l in range(number_of_spins - 1):
        B_ij = all_B_ij[l]
        r_ij = all_r_ij[l]
        for p in range(1, number_of_spins):
            orient_ij = orient[p - 1]
            sum += B_ij[p] * (ringrad /r_ij[p]) ** 2 * operator(vec, k, l, x, p)
    return sum, orient_ij

def summand_matrix_entry_distance_heis_1(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0
    for l in range(len(all_B_ij)):
        B_ij = all_B_ij[l]
        print("B_ij", B_ij)
        for p in range(len(B_ij)):
            sum += B_ij[p] * operator(vec, k, l, x, p)
    return sum
def summand_matrix_entry_distance_heis(operator, vec, k, x):
    all_B_ij, all_r_ij = summand_prefactor(number_of_spins)
    sum = 0
    for l in range(len(all_B_ij)):
        B_ij = all_B_ij[l]
        for p in range(len(B_ij)):
            sum += B_ij[p - 1] * operator(vec, k, l, x, p)
    return sum


def nhdw(): #nicht heisenbergsche Dipol-Wirkung, also der hintere Teil von dem dikken Hässlon
    matrix = - 3 * (dipol_1() + dipol_2() + dipol_3() + dipol_4())
    return matrix

def heisenberg_part():
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = summand_matrix_entry_distance_heis(s_plus_s_minus, vec, k, x)
            summand_2 = summand_matrix_entry_distance_heis(s_minus_s_plus, vec, k, x)
            summand_3 = summand_matrix_entry_distance_heis(sum_sz_sz_dipol, vec, k, x)

            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += 0.5 * (summand_1 + summand_2) + summand_3  # Hier wird die Berechnung in die Matrix geschrieben
    heisenberg_matrix = matrix
    return heisenberg_matrix
def dipol_dipol():
    H_dd = prefactor_sum_dd * (heisenberg_part() + nhdw())
    return H_dd

