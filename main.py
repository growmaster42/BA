import numpy as np

from Masochismus.expectation_values import expectation_value_sz
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

#show_table()
#dict1 = load_tuple_exp('data/expectation_values_spin_rings_all_pairs_second/spin=0.5_j_ij=-1_spins_max=3.json')
#plot_sys_exp_pairs(0.5, 3, dict1, -1)
for num in range(3, 8):
    create_latex_table(
        1, num, decimals=4)