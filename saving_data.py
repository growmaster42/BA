import numpy as np
import os
import math
import pickle

# Write a vector to a text file
def write_vector_to_txt(vector, directory, filename):
    os.makedirs(directory, exist_ok=True)  # Create directory if it doesn't exist
    full_path = os.path.join(directory, filename)
    with open(full_path, 'w') as f:
        vector_str = ' '.join(map(str, vector))
        f.write(vector_str + '\n')


# Read a vector from a text file in a custom directory
def read_vector_from_txt(directory, filename):
    full_path = os.path.join(directory, filename)
    with open(full_path, 'r') as f:
        line = f.readline()
        vector = np.array(list(map(float, line.split())))
        return vector


def export_ground_states(hamiltonian, spin, num_spins):
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    ground_state_1 = eigenvectors[:, 0]
    key = (spin, num_spins, 0)
    all_ground_states = {key: ground_state_1}
    for i, value in enumerate(eigenvalues):
        j = i + 1
        if j == len(eigenvalues):
            break
        elif math.isclose(eigenvalues[i], eigenvalues[j], abs_tol=1e-9):
            ground_state_n = eigenvectors[:, j]
            all_ground_states[(spin, num_spins, j)] = ground_state_n
        else:
            break
    return all_ground_states


def save_dictionary(dictionary, file_name, directory):
    # Ensure the directory exists
    os.makedirs(directory, exist_ok=True)

    # Create the full file path
    file_path = os.path.join(directory, file_name)
    with open(file_path, 'wb') as file:
        pickle.dump(dictionary, file)


def load_dict(file_name, directory):
    file_path = os.path.join(directory, file_name)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as file:
            return pickle.load(file)
    else:
        raise FileNotFoundError(f"No such file: {file_name}")


def store_systems(hamiltonian, min_spin, max_spin, spin):
    directory = "data/ground_state_vectors"
    if min_spin < 2:
        print("Minimum Spin must be => 2!")
    else:
        for num_spins in range(min_spin, max_spin + 1):
            dictionary = export_ground_states(hamiltonian(spin, num_spins), spin, num_spins)
            file_name = "spin=" + str(spin) + "_num_spins=" + str(num_spins)
            save_dictionary(dictionary, file_name, directory)
        print("Data has been stored successfully!!!!!")


def load_spin_system(spin, num_spins):
    directory = 'data/ground_state_vectors'
    file_name = "spin=" + str(spin) + "_num_spins=" + str(num_spins)
    spin_system = load_dict(file_name, directory)
    return spin_system


def calculate_eigenvalues(matrix):
    """
    Berechnet die Eigenwerte einer symmetrischen Matrix und speichert sie in einer Textdatei.

    Args:
        matrix: Eine symmetrische numpy-Array Matrix
    """
    # Überprüfe ob Matrix symmetrisch ist
    if not np.allclose(matrix, matrix.T):
        raise ValueError("Die eingegebene Matrix ist nicht symmetrisch!")

    # Berechne Eigenwerte
    eigenvalues = np.linalg.eigvalsh(matrix)

    # Speichere Eigenwerte in Textdatei
    with open('eigenvalues.txt', 'w') as f:
        f.write("Eigenwerte der Matrix:\n")
        for value in eigenvalues:
            f.write(f"{value}\n")



