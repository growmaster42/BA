import multiprocessing as mp
import logging
from interaction_models import spin_ring_hamiltonian
from Masochismus.old_versions.saving_data import store_systems

def parallel_store_systems(max_processes=None):
    # Define the parameters for each core
    params = [
        (spin_ring_hamiltonian, 2, 10, 0.5),
        (spin_ring_hamiltonian, 2, 7, 1.0),
        (spin_ring_hamiltonian, 2, 5, 1.5),
        (spin_ring_hamiltonian, 2, 5, 2.0),
        (spin_ring_hamiltonian, 2, 3, 2.5)
    ]

    # If max_processes is None, use the number of CPU cores available
    if max_processes is None:
        max_processes = mp.cpu_count()  # Get the number of CPU cores

    try:
        # Using context manager for proper resource cleanup
        with mp.Pool(processes=max_processes) as pool:
            # Execute the functions in parallel
            results = pool.starmap(store_systems, params)
            return results

    except Exception as e:
        logging.error(f"Error in parallel execution: {str(e)}")
        raise


if __name__ == '__main__':
    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Set the start method for macOS
    mp.set_start_method('fork')

    # Execute with error handling
    try:
        results = parallel_store_systems()
        logging.info("All calculations completed successfully")
    except Exception as e:
        logging.error(f"Failed to complete calculations: {str(e)}")