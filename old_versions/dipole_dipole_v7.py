from intial_values import *
from Masochismus.old_versions.S_plus_s_minus_s_z import *
import numpy as np
quarter = 0.25
quarter_i = - 0.25j

# function to create the orientation matrix that origins from the connection vectors between two spins
#returns 3x3 arrays within array with length of all possible interactions
def orientation():
    # initializing list for each matrix that occurs in the sum
    all_b_ij = []
    orientation_matrix_list = []
    unity_matrices = []
    all_r_ij = []
    unity = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=complex)
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
            unity_matrices.append(unity)
            r_ij.append(r)
            b_ij.append(b)
        if b_ij and r_ij:  # Stelle sicher, dass die Listen nicht leer sind
            all_b_ij.append(b_ij)
            all_r_ij.append(r_ij)

        # Winkel zur aktuellen Liste hinzufügen

    return orientation_matrix_list, unity_matrices, all_b_ij, all_r_ij
#transforming the 3x3 matrix elements to +,-,z, page 35 j.s. script
def transform_matrix(matrix):

    m, c, a, b = matrix

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
    c, m, a, b = orientation() #import of orientation function that returns 4 variables, second (m) is unity matrix

    all_matrices = np.array([])
    for l in range(len(c)):

        #row 1
        m[l][0, 0] = 0.5  #transformation for s+s+
        m[l][0, 1] = 0.5 #transformation for s+s-
        m[l][0, 2] = 0 #transformation for s-sz
        #row 2
        m[l][1, 0] = 0.5  #transformation for s-s-
        m[l][1, 1] = 0
        m[l][1, 2] = 0
        #row 3
        m[l][2, 0] = 0
        m[l][2, 1] = 0
        m[l][2, 2] = 1
        all_matrices = np.append(all_matrices, m[l])


    return m

#function that creates the sum for a single matrix entry
def summand_matrix_entry_distance_1(operator, vec, k, x, a, b):
    m, c, all_b_ij, all_r_ij = orientation()
    orient = transform_matrix(orientation())
    total_sum = 0
    for l in range(number_of_spins - 1):
        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            total_sum += orient[l][a, b] * b_ij[p - l - 1] * operator(vec, k, l, x, p)
    return total_sum

def nhdw_matrix():

    matrix = np.zeros((dim, dim), dtype=complex)
    vec = basvec()
    for k in range(dim):
        for x in range(dim):
            s_p_s_p = summand_matrix_entry_distance_1(s_plus_s_plus, vec, k, x, 0, 0)
            s_p_s_m = summand_matrix_entry_distance_1(s_plus_s_minus, vec, k, x, 0, 1)
            s_m_s_p = summand_matrix_entry_distance_1(s_minus_s_plus, vec, k, x, 1, 0)
            s_m_s_m = summand_matrix_entry_distance_1(s_minus_s_minus, vec, k, x, 1, 1)
            s_j_s = - 3 * (s_p_s_p + s_p_s_m + s_m_s_p + s_m_s_m)
            matrix[k, x] = s_j_s

    return matrix

def summand_matrix_entry_distance_heis(operator, vec, k, x, a, b):
    m, c,  all_b_ij, all_r_ij = orientation()
    orient = transform_matrix_unity(orientation())
    total_sum = 0
    for l in range(number_of_spins - 1):

        b_ij = all_b_ij[l]
        for p in range(l + 1, number_of_spins):
            total_sum += orient[l][a, b] * b_ij[p - l - 1] * operator(vec, k, l, x, p)

    return total_sum
def heisenberg_part():
    vec = basvec()

    matrix = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        for x in range(dim):
            spsm = summand_matrix_entry_distance_heis(s_plus_s_minus, vec, k, x, 0, 1) #s_i_plus * s_j_minus
            smsp = summand_matrix_entry_distance_heis(s_minus_s_plus, vec, k, x, 1, 0)    #s_i_minus * s_j_plus
            szsz = summand_matrix_entry_distance_heis(sz_sz_dipol, vec, k, x, 2, 2)   #s_i_z * s_j_z
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] = (spsm + smsp) + szsz  # Hier wird die Berechnung in die Matrix geschrieben
    heisenberg_matrix = matrix
    return heisenberg_matrix

def dipol_dipol():
    H_dd = prefactor_sum_dd * (nhdw_matrix() + heisenberg_part())
    return H_dd