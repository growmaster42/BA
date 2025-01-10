from intial_values import *
from S_plus_s_minus_s_z import *
import numpy as np
quarter = 0.25
quarter_i = - 0.25j
# function to create the orientation matrix that origins from the connection vectors between two spins
#returns 3x3 arrays within array with length of all possible interactions
def orientation():
    # initializing list for each matrix that occurs in the sum
    all_b_ij = []
    orientation_matrix_list = []
    all_r_ij = []
    #creating the angles that occur inbetween the spins by defining alpha_i and alpha_j
    for i in range(number_of_spins + 1):
        b_ij = []  # Liste für die aktuelle äußere Iteration
        r_ij = []
        alpha_i = 2 * np.pi * (i + 1) / number_of_spins
        for j in range(i + 2, number_of_spins + 1):
            alpha_j = 2 * np.pi * j / number_of_spins
            spin_angle = alpha_j - alpha_i
            r = 2 * ringrad * np.sin(spin_angle / 2)

            b = 1 / (r ** 3)

            frac = (ringrad / r) ** 2
            #abbreviate differences between cosines and sines of angles
            delta_cos = np.cos(alpha_j) - np.cos(alpha_i)
            delta_sin = np.sin(alpha_j) - np.sin(alpha_i)
            #defining the 3x3 matrix
            orientation_matrix = np.array([[frac * delta_cos ** 2,            frac * delta_cos * delta_sin, 0],
                               [frac * delta_sin * delta_cos,     frac * delta_sin ** 2,        0],
                               [0,                  0,                            0]], dtype=complex)
            orientation_matrix_list.append(orientation_matrix)
            r_ij.append(r)
            b_ij.append(b)
        if b_ij and r_ij:  # Stelle sicher, dass die Listen nicht leer sind
            all_b_ij.append(b_ij)
            all_r_ij.append(r_ij)

        # Winkel zur aktuellen Liste hinzufügen

    return orientation_matrix_list, all_b_ij, all_r_ij
#transforming the 3x3 matrix elements to +,-,z, page 35 j.s. script
def transform_matrix(matrix):

    m, a, b = matrix

    all_matrices = np.array([])
    for l in range(len(m)):

        j_xx = m[l][0, 0]
        j_xy = m[l][0, 1]
        j_xz = m[l][0, 2]
        j_yx = m[l][1, 0]
        j_yy = m[l][1, 1]
        j_yz = m[l][1, 2]
        j_zx = m[l][2, 0]
        j_zy = m[l][2, 1]
        j_zz = m[l][2, 2]
        #row 1
        m[l][0, 0] = quarter * j_xx + quarter_i * j_xy + quarter_i * j_yx - quarter * j_yy #transformation for s+s+
        m[l][0, 1] = quarter * j_xx - quarter_i * j_xy + quarter_i * j_yx + quarter * j_yy #transformation for s+s-
        m[l][0, 2] = 2 * quarter * j_xz + 2 * quarter_i * j_yz #transformation for s-s+
        #row 2
        m[l][1, 0] = quarter * j_xx + quarter_i * j_xy - quarter_i * j_yx + quarter * j_yy #transformation for s-s-
        m[l][1, 1] = quarter * j_xx - quarter_i * j_xy - quarter_i * j_yx - quarter * j_yy
        m[l][1, 2] = 2 * quarter * j_xz - 2 * quarter_i * j_yz
        #row 3
        m[l][2, 0] = 2 * quarter * j_zx + 2 * quarter_i * j_zy
        m[l][2, 1] = 2 * quarter * j_zx - 2 * quarter_i * j_zy
        m[l][2, 2] = j_zz
        all_matrices = np.append(all_matrices, [m])
    return m

def transform_matrix_unity(matrix):
    c, a, b = matrix
    m = [np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=complex)]

    all_matrices = np.array([])
    for l in range(len(c)):

        j_xx = m[l][0, 0]
        j_xy = m[l][0, 1]
        j_xz = m[l][0, 2]
        j_yx = m[l][1, 0]
        j_yy = m[l][1, 1]
        j_yz = m[l][1, 2]
        j_zx = m[l][2, 0]
        j_zy = m[l][2, 1]
        j_zz = m[l][2, 2]
        #row 1
        m[l][0, 0] = quarter * j_xx + quarter_i * j_xy + quarter_i * j_yx - quarter * j_yy #transformation for s+s+
        m[l][0, 1] = quarter * j_xx - quarter_i * j_xy + quarter_i * j_yx + quarter * j_yy #transformation for s+s-
        m[l][0, 2] = 2 * quarter * j_xz + 2 * quarter_i * j_yz #transformation for s-s+
        #row 2
        m[l][1, 0] = quarter * j_xx + quarter_i * j_xy - quarter_i * j_yx + quarter * j_yy #transformation for s-s-
        m[l][1, 1] = quarter * j_xx - quarter_i * j_xy - quarter_i * j_yx - quarter * j_yy
        m[l][1, 2] = 2 * quarter * j_xz - 2 * quarter_i * j_yz
        #row 3
        m[l][2, 0] = 2 * quarter * j_zx + 2 * quarter_i * j_zy
        m[l][2, 1] = 2 * quarter * j_zx - 2 * quarter_i * j_zy
        m[l][2, 2] = j_zz
        all_matrices = np.append(all_matrices, [m])
    return m

#function that creates the sum for a single matrix entry
def summand_matrix_entry_distance_1(operator, vec, k, x):
    m, all_b_ij, all_r_ij = orientation()
    orient = transform_matrix(orientation())
    total_sum = 0
    for l in range(number_of_spins - 1):
        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            total_sum += b_ij[p - l - 1] * operator(vec, k, l, x, p)

    return total_sum

def nhdw_matrix():
    orient = transform_matrix(orientation())
    matrix = np.zeros((dim, dim), dtype=complex)
    vec = basvec()
    for k in range(dim):
        for x in range(dim):
            s_p_s_p = orient[0][0, 0] * summand_matrix_entry_distance_1(s_plus_s_plus, vec, k, x)
            s_p_s_m = orient[0][0, 1] * summand_matrix_entry_distance_1(s_plus_s_minus, vec, k, x)
            s_m_s_p = orient[0][1, 0] * summand_matrix_entry_distance_1(s_minus_s_plus, vec, k, x)
            s_m_s_m = orient[0][1, 1] * summand_matrix_entry_distance_1(s_minus_s_minus, vec, k, x)
            s_j_s = - 3 * (s_p_s_p + s_p_s_m + s_m_s_p + s_m_s_m)
            matrix[k, x] += s_j_s
            if k == 1 and x == 2:
                print("s_p_s_p=", s_p_s_p, "s_p_s_m= ", s_p_s_m, "s_m_s_p= ",s_m_s_p, "s_m_s_m= ", s_m_s_m )
                print("spsp",summand_matrix_entry_distance_1(s_plus_s_minus, vec, k, x) )
    return matrix



#the functions dipol_1, dipol_2, dipol_3 and dipol_4 create various terms for the whole matrix
#k refers to row, x to column
def dipol_1():
    vec = basvec()
    m, all_b_ij, all_r_ij = orientation()
    orient = transform_matrix(orientation())
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = orient[0][0, 0] * summand_matrix_entry_distance_1(s_plus_s_plus, vec, k, x)
            summand_2 = orient[0][0, 1] * summand_matrix_entry_distance_1(s_plus_s_minus, vec, k, x)
            summand_3 = orient[0][1, 0] * summand_matrix_entry_distance_1(s_minus_s_plus, vec, k, x)
            summand_4 = orient[0][1, 1] * summand_matrix_entry_distance_1(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += (summand_1 + summand_2 + summand_3 + summand_4)
            if k == 3 and x == 0:
                print("dipol 1, matrix entry 3,0:" ,matrix[k, x])
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

def summand_matrix_entry_distance_2(operator, vec, k, x):
    m, all_b_ij, all_r_ij = orientation()
    orient = transform_matrix(orientation())
    sum = 0

    for l in range(number_of_spins - 1):
        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            sum += b_ij[p - l - 1] * operator(vec, k, l, x, p)

    return sum

def dipol_2():
    orient = transform_matrix(orientation())
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = orient[0][0, 0] * summand_matrix_entry_distance_2(s_plus_s_plus, vec, k, x)
            summand_2 = - orient[0][0, 1] * summand_matrix_entry_distance_2(s_plus_s_minus, vec, k, x)
            summand_3 = orient[0][1, 0] * summand_matrix_entry_distance_2(s_minus_s_plus, vec, k, x)
            summand_4 = - orient[0][1, 1] * summand_matrix_entry_distance_2(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += (summand_1 + summand_2 + summand_3 + summand_4)
            if k == 3 and x == 0:
                print("dipol 2, matrix entry 3,0:", matrix[k, x])
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix
def summand_matrix_entry_distance_3(operator, vec, k, x):
    m, all_b_ij, all_r_ij = orientation()
    orient = transform_matrix(orientation())
    sum = 0

    for l in range(number_of_spins - 1):
        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            sum += b_ij[p - l - 1] * operator(vec, k, l, x, p)

    return sum
def dipol_3():
    orient = transform_matrix(orientation())
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = orient[0][0, 0] * summand_matrix_entry_distance_3(s_plus_s_plus, vec, k, x)
            summand_2 = orient[0][0, 1] * summand_matrix_entry_distance_3(s_plus_s_minus, vec, k, x)
            summand_3 = - orient[0][1, 0] * summand_matrix_entry_distance_3(s_minus_s_plus, vec, k, x)
            summand_4 = - orient[0][1, 1] * summand_matrix_entry_distance_3(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += (summand_1 + summand_2 + summand_3 + summand_4)
            if k == 3 and x == 0:
                print("dipol 3, matrix entry 3,0:" ,matrix[k, x])
    # Hier wird die Berechnung in die Matrix geschrieben
    return matrix
def summand_matrix_entry_distance_4(operator, vec, k, x):
    m, all_b_ij, all_r_ij = orientation()
    orient = transform_matrix(orientation())
    sum = 0
    for l in range(number_of_spins - 1):
        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            sum += b_ij[p - l - 1] * operator(vec, k, l, x, p)
    return sum
def dipol_4():
    orient = transform_matrix(orientation())
    vec = basvec()
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = - orient[0][0, 0] * summand_matrix_entry_distance_4(s_plus_s_plus, vec, k, x)
            summand_2 = orient[0][0, 1] * summand_matrix_entry_distance_4(s_plus_s_minus, vec, k, x)
            summand_3 = orient[0][1, 0] * summand_matrix_entry_distance_4(s_minus_s_plus, vec, k, x)
            summand_4 = -  orient[0][1, 1] * summand_matrix_entry_distance_4(s_minus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += (summand_1 + summand_2 + summand_3 + summand_4)
            if k == 3 and x == 0:
                print("dipol 4, matrix entry 3,0:" ,matrix[k, x])
            # Hier wird die Berechnung in die Matrix geschrieben
    return matrix
def summand_matrix_entry_distance_heis(operator, vec, k, x):
    orient, all_b_ij, all_r_ij = orientation()
    sum = 0
    for l in range(number_of_spins - 1):
        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            sum += b_ij[p - l - 1] * operator(vec, k, l, x, p)
    return sum
def heisenberg_part():
    vec = basvec()
    orient = transform_matrix_unity(orientation())
    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            summand_1 = orient[0][0, 1] * summand_matrix_entry_distance_heis(s_plus_s_minus, vec, k, x) #s_i_plus * s_j_minus
            summand_2 = orient[0][1, 0] * summand_matrix_entry_distance_heis(s_minus_s_plus, vec, k, x)    #s_i_minus * s_j_plus
            summand_3 = orient[0][2, 2] * summand_matrix_entry_distance_heis(sz_sz_dipol, vec, k, x)   #s_i_z * s_j_z
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += (summand_1 + summand_2) + summand_3  # Hier wird die Berechnung in die Matrix geschrieben
    heisenberg_matrix = matrix
    return heisenberg_matrix
def nhdw(): #nicht heisenbergsche Dipol-Wirkung, also der hintere Teil von dem dikken Hässlon
    matrix = - 3 * (dipol_1() + dipol_2() + dipol_3() + dipol_4())
    return matrix
def dipol_dipol():
    H_dd = prefactor_sum_dd * (nhdw_matrix() +heisenberg_part())
    return H_dd