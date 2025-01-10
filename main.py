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

spin = 1.5
num_spins = 5
eigenvalue = np.linalg.eigvalsh(spin_ring_hamiltonian(spin, num_spins))
print(f"Eigenvalues of the Hamiltonian for {num_spins} spins with spin {spin}: \n", eigenvalue)

#degeneracy_plot()
#print((1/0.19245) * s_j_s_heis(basvec(0.5, 3), 0.5, 3))













