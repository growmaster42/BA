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

#print(load_system(0.5, 2, -1))

