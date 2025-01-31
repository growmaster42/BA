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

def expectation_value_sz(spin, number_of_spins, eigenvector, spin_pair):
    num_spins = number_of_spins
    i, j = spin_pair
    basis_vectors = basvec(spin, number_of_spins)
    psi = np.zeros(len(eigenvector), dtype=complex)
    exp_value = 0
    for k, eig_v in enumerate(eigenvector):
        vec_num = vector_number(basis_vectors[k], spin)
        # sz_sz
        sz = s_z(basis_vectors, spin, k, i)
        psi[vec_num] += sz * eigenvector[vec_num]
        exp_value +=sz * eig_v.conj() * eig_v
        print(exp_value)
    return exp_value

def quantum_number_i(spin, number_of_spins, eigenvector):
    m_i = 0
    for i in range(number_of_spins):
        m_i += expectation_value_sz(spin, number_of_spins, eigenvector, (i, 0))
    return m_i

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


import numpy as np


def compute_sz_expectation(state_vector):
    """
    Compute the expectation value of S^z operator for a given quantum state vector.

    Parameters:
    -----------
    state_vector : numpy.ndarray
        The quantum state vector (must be normalized)
        Shape should be (2^n,) where n is the number of qubits/spins

    Returns:
    --------
    float
        The expectation value <ψ|S^z|ψ>
    """
    # Check if input is a numpy array
    if not isinstance(state_vector, np.ndarray):
        state_vector = np.array(state_vector)

    # Get number of qubits/spins from state vector dimension
    n_qubits = int(np.log2(len(state_vector)))
    if 2 ** n_qubits != len(state_vector):
        raise ValueError("State vector length must be a power of 2")

    # Check normalization
    if not np.isclose(np.sum(np.abs(state_vector) ** 2), 1.0, atol=1e-6):
        raise ValueError("State vector must be normalized")

    # Initialize expectation value
    expectation = 0.0

    # Calculate expectation value for each qubit
    for i in range(n_qubits):
        # Create the S^z operator matrix for the i-th qubit
        sz_i = np.zeros((2 ** n_qubits, 2 ** n_qubits), dtype=complex)

        # Fill the diagonal elements
        for j in range(2 ** n_qubits):
            # Check if i-th bit is 0 or 1
            if (j >> i) & 1:
                sz_i[j, j] = -0.5  # Spin down
            else:
                sz_i[j, j] = 0.5  # Spin up

        # Add contribution from this qubit
        expectation += np.real(state_vector.conj() @ sz_i @ state_vector)

    return expectation


# Example usage:
if __name__ == "__main__":
    # Example: Single qubit in the state |0⟩
    psi_0 = np.array([1, 0], dtype=complex)
    print(f"Expectation value for |0⟩: {compute_sz_expectation(psi_0)}")

    # Example: Single qubit in the state |1⟩
    psi_1 = np.array([0, 1], dtype=complex)
    print(f"Expectation value for |1⟩: {compute_sz_expectation(psi_1)}")

    # Example: Two-qubit state (|00⟩ + |11⟩)/√2
    bell_state = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    print(f"Expectation value for Bell state: {compute_sz_expectation(bell_state)}")

















