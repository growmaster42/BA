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


def load_spin_system(spin, num_spins, directory):
    file_name = "spin=" + str(spin) + "_num_spins=" + str(num_spins)
    spin_system = load_dict(file_name, directory)
    return spin_system


import numpy as np
import json
from pathlib import Path
from datetime import datetime


def calculate_degeneracy(eigenvalues, tolerance=1e-10):
    """
    Calculate the degeneracy of the ground state.

    Args:
        eigenvalues: Sorted array of eigenvalues (ascending order)
        tolerance: Numerical tolerance for comparing eigenvalues

    Returns:
        int: Degree of degeneracy of the ground state
    """
    if len(eigenvalues) == 0:
        return 0

    degeneracy = 1
    ground_state = eigenvalues[0]

    for i in range(1, len(eigenvalues)):
        if abs(eigenvalues[i] - ground_state) <= tolerance:
            degeneracy += 1
        else:
            break

    return degeneracy


def process_and_store_eigenvalues(matrix, spin, num_spins, j_ij):
    """
    Calculate eigenvalues, determine ground state degeneracy, and store in JSON file.

    Args:
        matrix: Input matrix to calculate eigenvalues from
        spin: Spin value
        num_spins: Number of spins in the system
        j_ij: Coupling constant
    """
    # Create data directory if it doesn't exist
    data_dir = Path('data/eigenvalues')
    data_dir.mkdir(parents=True, exist_ok=True)

    # Create filename based on parameters
    filename = f"spin={spin}_num_spins={num_spins}_J_ij={j_ij}.json"
    filepath = data_dir / filename

    # Calculate eigenvalues and sort them
    eigenvalues = np.sort(np.linalg.eigvals(matrix).real)

    # Calculate ground state degeneracy
    dg_deg_ground = calculate_degeneracy(eigenvalues)

    # Create the dictionary key as a tuple
    dict_key = (float(spin), int(num_spins), float(j_ij), int(dg_deg_ground))

    # Convert eigenvalues to a list and round to avoid numerical artifacts
    eigenvalues_list = [float(round(ev, 10)) for ev in eigenvalues]

    # Prepare the data structure
    data = {}
    if filepath.exists():
        with open(filepath, 'r') as f:
            data = json.load(f)

    # Add metadata to the storage
    metadata = {
        "last_modified": datetime.utcnow().isoformat(),
        "last_modified_by": "growmaster42",
        "parameters": {
            "spin": spin,
            "num_spins": num_spins,
            "j_ij": j_ij
        }
    }

    # Store both metadata and eigenvalues
    key_str = str(dict_key)
    if "metadata" not in data:
        data["metadata"] = metadata
    else:
        data["metadata"].update(metadata)

    data[key_str] = eigenvalues_list

    # Save updated data with pretty printing
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=4)
    print("-----------------------------------------------------------"
              "\n Eigenvalues for spin system with", num_spins,
              "spins have been calculated and stored. \n "
              "-----------------------------------------------------------")





