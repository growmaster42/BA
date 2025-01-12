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

import multiprocessing as mp
from multiprocessing import Process

def worker(spin_range):
    """Worker function that processes a specific spin value and its range"""
    spin, start, end = spin_range
    for num_spins in range(start, end):
        matrix = spin_ring_hamiltonian(spin, num_spins)
        process_and_store_eigenvalues(matrix, spin, num_spins, J_ij)

def main():
    # Define the work for each core as (spin_value, start_range, end_range)
    spin_ranges = [
        (0.5, 2, 11),  # core 1
        (1.0, 2, 7),   # core 2
        (1.5, 2, 6),   # core 3
        (2.0, 2, 6),   # core 4
        (2.5, 2, 5),   # core 5
    ]

    # Create processes
    processes = []
    for spin_range in spin_ranges:
        p = Process(target=worker, args=(spin_range,))
        processes.append(p)

    # Start all processes
    for p in processes:
        p.start()

    # Wait for all processes to complete
    for p in processes:
        p.join()

if __name__ == '__main__':
    # This guard is important for multiprocessing on macOS
    main()






