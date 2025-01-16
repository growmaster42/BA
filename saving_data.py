import multiprocessing as mp
from numpy import linalg as lg
from interaction_models import *
from numpy import isclose
import pickle
import os

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




