
from interaction_models import *
from operators import *
from time_estimation import main
import time as tm
from basis_vectors import basvec
import json
from misc import *
from saving_data import *
from expectation_values import *


spin = 1
num_spins = 4
dg_deg = 0
spin_pair = (1, 2)

eigenstate = load_spin_system(spin, num_spins)[(spin, num_spins, dg_deg)]

spin_pair_system, sz_sz, spsm, smsp = spin_pair_operator(spin, spin_pair, num_spins)
exp_value = np.real(np.dot(eigenstate.conj(), np.dot(spin_pair_system, eigenstate)))
exp_value_sz_sz = np.real(np.dot(eigenstate.conj(), np.dot(sz_sz, eigenstate)))
exp_value_spsm = np.real(np.dot(eigenstate.conj(), np.dot(spsm, eigenstate)))
exp_value_smsp = np.real(np.dot(eigenstate.conj(), np.dot(smsp, eigenstate)))

print("expectation value from matrix operation \U0001F913 : \n ", exp_value)
print("expectation value my function s_i_s_i+1 lüpt \U0001F44D : \n", expectation_value_s_i_s_i_plus_one_spin_ring(spin, num_spins, spin_pair, dg_deg))














