import multiprocessing as mp
import numpy as np
from interaction_models import spin_ring
import pickle
import os
import json
from expectation_values import expectation_value_s_i_s_i_plus_one
import ast
def spin_system_to_dict(spin, num_spins, j_ij):
    """This function takes the spin, the number of spins and the interaction model and returns a dictionary with the
    eigenvalues, eigenvectors and degeneracies of the system. The key of the dictionary is a tuple of the function name,
    the spin, the number of spins and J_ij. The degeneracies are counted by rounding the eigenvalues to 9 decimal
    places and counting the unique values."""
    eigenvalues, eigenvectors = np.linalg.eigh(spin_ring(spin, num_spins, j_ij))
    # round eigenvalues for degeneracy counting
    rounded_eigenvalues = np.round(eigenvalues, 9)
    # count the degeneracies
    unique, counts = np.unique(rounded_eigenvalues, return_counts=True)
    # add counts to the dictionary, data structure of degeneracies is a list where the position of the degeneracy refers
    # to the actual energy level or state
    dictionary =  {'eigenvalues' : eigenvalues,
                         'eigenvectors' : eigenvectors,
                         'degeneracies' : counts}
    # testing if the degeneracies are correct by summing them up and comparing them to the len(eigenvalues)
    assert sum(counts) == len(eigenvalues) # assert means that the condition has to be true, otherwise the program stops
    return dictionary

def save_system(spin, num_spins, j_ij):
    system = spin_system_to_dict(spin, num_spins, j_ij)
    filename = f"spin_{spin}_num_spins_{num_spins}_j_ij_{j_ij}.pickle"
    path = os.getcwd()
    folder_path = os.path.join(path, "saved_systems")
    os.makedirs(folder_path, exist_ok=True)
    file_path = os.path.join(folder_path, filename)

    with open(file_path, 'wb') as f:
        pickle.dump(system, f)
    print(f"System saved as {filename}")

def load_system(spin, num_spins, j_ij):
    filename = f"spin_{spin}_num_spins_{num_spins}_j_ij_{j_ij}.pickle"
    path = os.getcwd()
    folder_path = os.path.join(path, "saved_systems")
    file_path = os.path.join(folder_path, filename)
    with open(file_path, 'rb') as f:
        system = pickle.load(f)
    return system

def save_1(min, max):
    for n in range(min, max + 1):
        for j in range(-1, 2):
            save_system(0.5, n, j)
            print(f"-------------------:{j}----------------------")
    print("OOOOOOOOOOOOOOspin 0.5 doneOOOOOOOOOOOOOO")

def save_2(min, max):
    for n in range(min, max + 1):
        for j in range(-1, 2):
            save_system(1, n, j)
            print(f"-------------------:{j}----------------------")
    print("OOOOOOOOOOOOOOspin 1 doneOOOOOOOOOOOOOO")


def save_3(min, max):
    for n in range(min, max + 1):
        for j in range(-1, 2):
            save_system(1.5, n, j)
            print(f"-------------------:{j}----------------------")
    print("OOOOOOOOOOOOOOspin 1.5 doneOOOOOOOOOOOOOO")

def save_4(min, max):
    for n in range(min, max + 1):
        for j in range(-1, 2):
            save_system(2, n, j)
            print(f"-------------------:{j}----------------------")
    print("OOOOOOOOOOOOOOspin 2 doneOOOOOOOOOOOOOO")

def save_5(min, max):
    for n in range(min, max + 1):
        for j in range(-1, 2):
            save_system(2.5, n, j)
            print(f"-------------------:{j}----------------------")
    print("OOOOOOOOOOOOOOspin fully 2.5 doneOOOOOOOOOOOOOO")


def run_parallel_saves_explicit(min_value):
    # Create process objects
    processes = [
        mp.Process(target=save_1, args=(min_value, 10)),
        mp.Process(target=save_2, args=(min_value, 7)),
        mp.Process(target=save_3, args=(min_value, 5)),
        mp.Process(target=save_4, args=(min_value, 5)),
        mp.Process(target=save_5, args=(min_value, 4))
    ]

    # Start all processes
    for p in processes:
        p.start()

    # Wait for all processes to complete
    for p in processes:
        p.join()

def state_degeneracy(spin, num_spins, j_ij, state):
    spin_system = load_system(spin, num_spins, j_ij)
    degeneracies = spin_system['degeneracies']
    if state >= len(degeneracies):
        return ValueError("State not in system")
    else:
        return degeneracies[state]

def collect_degeneracies(spin, j_ij, state):
    key_1 = (spin, j_ij)
    all_degeneracies = {key_1: {}}
    if spin == 0.5:
        num_spins = 10
    elif spin == 1:
        num_spins = 7
    elif spin == 1.5:
        num_spins = 5
    elif spin == 2:
        num_spins = 5
    elif spin == 2.5:
        num_spins = 4
    for n in range(2, num_spins + 1):
        key_2 = n
        degeneracy = state_degeneracy(spin, n, j_ij, state)
        all_degeneracies[key_1][key_2] = degeneracy
    return all_degeneracies


def make_exp_dict(spin, j_ij, spins_min, spins_max):
    exp_dict = {}
    for num_spins in range(spins_min, spins_max + 1):
        spin_pair=(0, 1)
        eigenvector =load_system(spin, num_spins, j_ij)['eigenvectors'][:, 0]
        exp_value = np.real(expectation_value_s_i_s_i_plus_one(spin, num_spins, eigenvector, spin_pair))
        exp_dict[num_spins] = exp_value
    return exp_dict


def expec_value_all_pairs_to_dict(spin, number_of_spins, j_ij, state):
    """Calculate the expectation values of all pairs of spins in a ring and return them as a dictionary."""
    # Pre-load eigenvectors to avoid recomputation in the loop
    eigenvector = load_system(spin, number_of_spins, j_ij)['eigenvectors'][:, state]

    # Create the dictionary with converted keys inline
    exp_value_dict = {
        str((i, (i + 1) % number_of_spins)):
            np.real(
                expectation_value_s_i_s_i_plus_one(spin, number_of_spins, eigenvector, (i, (i + 1) % number_of_spins)))
        for i in range(number_of_spins)
    }

    return exp_value_dict


def save_spin_rings_single_expect(data, spin, j_ij, spins_min , spins_max):
    # Ensure the directory exists
    directory = 'data/expectation_values_spin_rings'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Define the file path
    file_path = os.path.join(directory,
        f'spin={spin}_j_ij={j_ij}_spins_min={spins_min}_spins_max={spins_max}.json')

    # Write the dictionary to a JSON file
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

    print(f"Dictionary successfully saved to {file_path}")

def save_spin_rings_all_pairs_exp(data, spin, j_ij, spins_max, directory='data/expectation_values_spin_rings_all_pairs'):
    """Save the expectation values of all pairs of spins in a ring to a JSON file.
    Args:
        data (dict): Dictionary containing the expectation values of all pairs of spins,
                     function expectation_values_all_pairs can be used here
        spin (float): Spin of the system
        j_ij (float): Coupling constant
        spins_max (int): Maximum number of spins in the ring
        Returns:
                None """
    # Ensure the directory exists
    if not os.path.exists(directory):
        os.makedirs(directory)
    # Define the file path
    file_path = os.path.join(directory,
        f'spin={spin}_j_ij={j_ij}_spins_max={spins_max}.json')

    # Write the dictionary to a JSON file
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

    print(f"Dictionary successfully saved to {file_path}")


def load_json(path, name):
    # Construct the full file path
    file_path = os.path.join(path, name)

    # Open the file and load the JSON data
    with open(file_path, 'r') as file:
        data = json.load(file)

    return {float(k) if '.' in k else int(k): v for k, v in data.items()}


def convert_keys_to_tuple(data):
    """Convert string keys to tuples in a nested dictionary or list.
    Args:
        data (dict or list): Dictionary or list containing string keys
    Returns:
        dict or list: Dictionary or list with string keys converted to tuples"""
    if isinstance(data, dict):
        return {ast.literal_eval(k): convert_keys_to_tuple(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [convert_keys_to_tuple(item) for item in data]
    else:
        return data


def load_tuple_exp(file_path):
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)

    # Convert string keys back to tuples
    data_with_tuples = {ast.literal_eval(k): v for k, v in data.items()}

    return data_with_tuples


def create_latex_table(s, num, decimals=6):
    #print(f"Starting function with s={s}, num={num}, decimals={decimals}")  # Debug print

    # Dictionary to store data for different j_ij values
    data_dict = {}

    # Iterate over j_ij values
    for j_ij in range(-1, 2):
        filename = f"data/expectation_values_spin_rings_all_pairs/spin={s}_j_ij={j_ij}_spins_max={num}.json"
        #print(f"Attempting to open file: {filename}")  # Debug print

        try:
            with open(filename, 'r') as f:
                data = json.load(f)
                # Convert string representation of tuples back to actual tuples and round the values
                data = {eval(k): round(v, decimals) for k, v in data.items()}
                data_dict[j_ij] = data
                #print(f"Successfully loaded data for j_ij={j_ij}")  # Debug print
        except FileNotFoundError:
            print(f"Warning: File {filename} not found!")
            return

    # Get all unique pairs (keys) from all dictionaries
    all_pairs = set()
    for j_data in data_dict.values():
        all_pairs.update(j_data.keys())
    all_pairs = sorted(all_pairs)

    #print(f"Found {len(all_pairs)} unique pairs")  # Debug print

    latex_output = [
        "\\begin{table}[h]",
        "\\centering",
        "\\begin{tabular}{|c|c|c|c|}",
        "\\hline",
        "Pair & $J_{ij}=-1$ & $J_{ij}=0$ & $J_{ij}=1$ \\\\",
        "\\hline"
    ]

    for pair in all_pairs:
        row = [
            f"({pair[0]},{pair[1]})",
            f"{data_dict[-1].get(pair, '-'):.{decimals}f}",
            f"{data_dict[0].get(pair, '-'):.{decimals}f}",
            f"{data_dict[1].get(pair, '-'):.{decimals}f}"
        ]
        latex_output.append(" & ".join(row) + " \\\\")
        latex_output.append("\\hline")

    latex_output.extend([
        "\\end{tabular}",
        f"\\caption{{Spin ring with s={s} and \\#spins = {num}}}",
        "\\end{table}"
    ])

    #print("\nGenerated LaTeX table:")  # Debug print
    for line in latex_output:
        print(line)


# Now actually call the function with some parameters:
if __name__ == "__main__":
    print("Starting program")
    create_latex_table(s=0.5, num=4, decimals=3)
