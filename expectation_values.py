from basis_vectors import *
from operators import *
from interaction_models import *
from misc import load_spin_system

def expectation_value_sz_sz(spin, number_of_spins, eigenvector, spin_pair):
    num_spins = number_of_spins
    i, j = spin_pair
    basis_vectors = basvec(spin, number_of_spins)
    psi = np.zeros(len(eigenvector), dtype=complex)
    for k, eig_v in enumerate(eigenvector):
        vec_num = vector_number(basis_vectors[k], spin)
        # sz_sz
        sz = s_z(basis_vectors, spin, k, i) * s_z(basis_vectors, spin, k, j)
        psi[vec_num] += sz * eigenvector[vec_num]


    return np.dot(eigenvector.conj(), psi)


def expectation_value_s_plus_s_minus(spin, number_of_spins, eigenvector, spin_pair):
    num_spins = number_of_spins
    erw =0.0
    i, j = spin_pair
    basis_vectors = basvec(spin, number_of_spins)
    psi = np.zeros(len(eigenvector), dtype=complex)
    for k, eig_v in enumerate(eigenvector):
        vec_num = vector_number(basis_vectors[k], spin)
        for l, eig_w in enumerate(eigenvector):

            # s_plus_s_minus_operator
            if any(value > 2 * spin or value < 0 for value in s_plus_s_minus_ket_change(basis_vectors, k, i, j)):
                s_p_s_m = 0
            elif np.array_equal(basis_vectors[l], s_plus_s_minus_ket_change(basis_vectors, k, i, j)):

                s_p_s_m = s_plus(basis_vectors, spin, k, i) * s_minus(basis_vectors, spin, k, j)
            else:
                s_p_s_m = 0

            #psi[vec_num] += s_p_s_m * eigenvector[vec_num]
            erw+=s_p_s_m * eigenvector[l].conj()*eigenvector[k]

    #return np.dot(eigenvector.conj(), psi)
    return erw


def expectation_value_s_minus_s_plus(spin, number_of_spins, eigenvector, spin_pair):
    i, j = spin_pair
    basis_vectors = basvec(spin, number_of_spins)
    psi = np.zeros(len(eigenvector), dtype=complex)
    erw=0.0
    for k, eig_v in enumerate(eigenvector):
        vec_num = vector_number(basis_vectors[k], spin) #ist das nicht einfach k?
        for l, eig_w in enumerate(eigenvector):
            # s_plus_s_minus_operator
            if any(value > 2 * spin or value < 0 for value in s_minus_s_plus_ket_change(basis_vectors, k, i, j)):
                s_p_s_m = 0
            elif np.array_equal(basis_vectors[l], s_minus_s_plus_ket_change(basis_vectors, k, i, j)):

                s_p_s_m = s_minus(basis_vectors, spin, k, i) * s_plus(basis_vectors, spin, k, j)
            else:
                s_p_s_m = 0

            #psi[vec_num] += s_p_s_m * eigenvector[vec_num]
            erw+=s_p_s_m * eigenvector[l].conj()*eigenvector[k]

    #return np.dot(eigenvector.conj(), psi)
    return erw


def expectation_value_s_i_s_i_plus_one(spin, number_of_spins, eigenvector, spin_pair):
    i, j = spin_pair
    exp_value = expectation_value_sz_sz(spin, number_of_spins, eigenvector, spin_pair) \
                + 0.5 * expectation_value_s_minus_s_plus(spin, number_of_spins, eigenvector, spin_pair) \
                + 0.5 * expectation_value_s_plus_s_minus(spin, number_of_spins, eigenvector, spin_pair)
    return exp_value



def expectation_value_s_i_s_i_plus_one_spin_ring(spin, num_spins, spin_pair, dg_deg):
    eigenvector = load_spin_system(spin, num_spins)[(spin, num_spins, dg_deg)]
    szsz = expectation_value_sz_sz(spin, num_spins, eigenvector, spin_pair)
    spsm = expectation_value_s_plus_s_minus(spin, num_spins, eigenvector, spin_pair)
    smsp = expectation_value_s_minus_s_plus(spin, num_spins, eigenvector, spin_pair)
    exp_value =   szsz + 0.5 * spsm + 0.5 * smsp

    return np.real(exp_value)


def all_exp_values(spin, max_spins, dg_deg):
    values = {}
    for i in range(max_spins - 1):
        j = i + 1
        spin_pair = i, j
        values[i, i + 1] = expectation_value_s_i_s_i_plus_one_spin_ring(spin, max_spins, spin_pair, dg_deg)
    return values




















