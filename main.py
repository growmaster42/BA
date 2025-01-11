from Masochismus.time_estimation import estimate_time
from plotting import *
from interaction_models import *
from operators import *
from time_estimation import main
import time as tm
from basis_vectors import basvec
import json
from misc import *
from saving_data import *
from expectation_values import *
from distance_calculations import *

store_systems(spin_ring_hamiltonian, 2, 10, 0.5)
print("done")
store_systems(spin_ring_hamiltonian, 2, 7, 1)
print("done")
store_systems(spin_ring_hamiltonian, 2, 5, 1.5)
print("done")
store_systems(spin_ring_hamiltonian, 2, 5, 2)
print("done")
store_systems(spin_ring_hamiltonian, 2, 4, 2.5)
print("done")













