from Masochismus.old_versions.S_plus_s_minus_s_z import *


def Heisenberg_slow():
    vec = basvec()
    matrix = - 2 * np.diag(sum_sz_i_sz_j(vec, number_of_spins))
    for k in range(dim):
        for x in range(dim):
            summand_1 = summand_matrix_entry(s_minus_s_plus, vec, k, x)
            summand_2 = summand_matrix_entry(s_plus_s_minus, vec, k, x)
            # Ersetze dies ggf. durch eine andere Berechnung
            matrix[k, x] += - (summand_1 + summand_2)  # Hier wird die Berechnung in die Matrix geschrieben
    return matrix

