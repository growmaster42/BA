import numpy as np
from concurrent.futures import ProcessPoolExecutor

from Masochismus.expectation_values import expectation_value_sz, quantum_number_i, compute_sz_expectation
from Masochismus.interaction_models import heisenberg_ring
from Masochismus.old_versions.saving_data import *
from interaction_models import *
from saving_data import *
from plotting import *
from make_table import *
from arrows import *

spin = 0.5
n = 5
j_ij = - np.pi
hamiltonian = heisenberg_ring(basvec(spin, n), spin, n, j_ij)
eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)

def save_array_to_file(array, filename="output.txt"):
    """Writes an array to a text file, each element on a new line."""
    with open(filename, "w") as file:
        for item in array:
            file.write(f"{item}\n")

save_array_to_file(eigenvalues, f"eigenvalues_of_heis_ring_s={spin}_n={n}_jij={j_ij}.txt")