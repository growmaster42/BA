import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Dictionary to store data for different spin systems
spin_systems_data = {
    0.5: {
        "n_values": np.array([2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "times": np.array([0.005464076995849609, 0.010461091995239258, 0.06242799758911133, 0.30720996856689453, 1.8171741962432861, 10.644792795181274, 60.47462773323059, 338.1103768348694, 1726.8112230300903])
    },
    1: {
            "n_values": np.array([2, 3, 4]),
            "times": np.array([0.002875089645385742, 0.04745030403137207, 0.6411268711090088])  # Example data for spin 1
        },
    1.5: {
                "n_values": np.array([2, 3, 4]),
                "times": np.array([0.02449178695678711, 0.43862295150756836, 12.855217933654785])  # Example data for spin 1
        },
    2.0 : { "n_values": np.array([2, 3, 4, 5]),
                "times": np.array([0.02277088165283203, 0.9565587043762207, 42.07624816894531, 1682.0156149864197])},
    2.5: {
        "n_values": np.array([2, 3, 4]),
        "times": np.array([0.07218098640441895, 4.7457239627838135, 358.77376198768616])  # Example data for spin 1
    },

    # Add more spin systems as needed
}

def estimate_time(spin, num_spins):
    """
    Estimates computation time for a given spin system and number of spins.

    Parameters:
        spin (float): The spin system (e.g., 0.5, 1, 1.5).
        num_spins (int): The number of spins to estimate for.

    Returns:
        float: Estimated computation time in seconds.
    """
    if spin not in spin_systems_data:
        print(f"Spin system {spin} is not available. Please add data first.")
        return None

    # Get data for the selected spin system
    data = spin_systems_data[spin]
    n_values, times = data["n_values"], data["times"]

    # Define the exponential model
    def model(n, a, b):
        return a * b**n

    # Fit the model
    params, _ = curve_fit(model, n_values, times)
    a, b = params

    # Predict time for the given number of spins
    estimated_time = model(num_spins, a, b)

    # Plot the data and the fit
    n_fit = np.linspace(min(n_values), max(num_spins, max(n_values)), 100)
    time_fit = model(n_fit, a, b)

    plt.scatter(n_values, times, label="Data", color="blue")
    plt.plot(n_fit, time_fit, label="Fit", color="red")
    plt.title(f"Spin {spin} System Runtime Estimate")
    plt.xlabel("Number of Spins (n)")
    plt.ylabel("Computation Time (s)")
    plt.yscale("log")
    plt.legend()


    return estimated_time

def main():
    """
    Main function to interactively estimate computation time.
    """
    print("Available spin systems:", list(spin_systems_data.keys()))
    try:
        spin = float(input("Enter the spin system you want to estimate (e.g., 0.5, 1, 1.5): "))
        num_spins = int(input("Enter the number of spins: "))

        estimated_time = estimate_time(spin, num_spins)
        if estimated_time is not None:
            print(f"Estimated computation time for spin {spin} with {num_spins} spins: {estimated_time:.2f} seconds")
    except ValueError:
        print("Invalid input. Please enter numerical values.")

if __name__ == "__main__":
    main()
