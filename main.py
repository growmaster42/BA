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
spin = 0.5
num_spins = 10
hamiltonian = spin_ring_hamiltonian(spin, num_spins)
eigenvalue, eigenvector = np.linalg.eigh(hamiltonian)
print(eigenvalue[:2])












