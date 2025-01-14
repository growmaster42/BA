import matplotlib.pyplot as plt
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
import multiprocessing as mp
import logging
from functools import partial
matrix = dipole_dipole(basvec(0.5, 8), 0.5, 8)
process_and_store_eigenvalues(matrix, 0.5, 8, 0)










