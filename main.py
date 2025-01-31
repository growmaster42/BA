import numpy as np

from Masochismus.expectation_values import expectation_value_sz, quantum_number_i, compute_sz_expectation
from Masochismus.interaction_models import heisenberg_ring
from Masochismus.old_versions.saving_data import *
from interaction_models import *
from saving_data import *
from concurrent.futures import ProcessPoolExecutor
import time
from plotting import *
from make_table import *

def show_table():
    for j_ij in range(-1, 2):
        print(f"J_ij = {j_ij}")
        data = {}
        for num_spins in range(3, 8):
            eigenvalue = load_system(1, num_spins, j_ij)['eigenvalues'][0]
            data[str(num_spins)] = eigenvalue
        print(data)


eigenvector = load_system(0.5, 3, j_ij=-1)['eigenvectors'][0]
exp_val = quantum_number_i(0.5, 3, eigenvector)
print("exp val function", exp_val)

print(compute_sz_expectation(eigenvector))