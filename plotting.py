import os
from expectation_values import all_exp_values
from pathlib import Path
from saving_data import collect_degeneracies
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


def degeneracy_plot(spin, j_ij, state):
    # Define the data_jij_m within the function, the data_jij_m refers to the coupling constant J_ij which is set to -1
    data = collect_degeneracies(spin, j_ij, state)
    # Create figure and axis
    fig, ax = plt.subplots()

    # Extract x and y values from the dictionary
    x_values = list(data[(spin, j_ij)].keys())
    y_values = list(data[(spin, j_ij)].values())

    # Create the scatter plot
    plt.scatter(x_values, y_values, color='blue')

    # Set the x-axis range from 2 to 10
    plt.xlim(1, max(data[(spin, j_ij)].keys()) + 1)

    # Set the y-axis range from 0 to 3
    plt.ylim(0, 3)

    # Add grid for better readability
    plt.grid(True, linestyle='--', alpha=0.7)

    state_name = r"Ground State" if state == 0 else  str(state) + r". Excited State"
    # Add labels and title
    plt.xlabel('Number of Spins')
    plt.ylabel('Degree of Degeneracy')
    if spin == 0.5:
        plt.title(state_name + r" Degeneracy for $s=\frac{1}{2}$ and $J_{ij}=$" + f"{j_ij}")
    elif spin == 1:
        plt.title(state_name + r" Degeneracy for $s=1$ and $J_ij=$" + f"{j_ij}")
    elif spin == 1.5:
        plt.title(state_name + r" Degeneracy for $s=\frac{3}{2}$ and $J_{ij}=$" + f"{j_ij}")
    elif spin == 2:
        plt.title(state_name + r" Degeneracy for $s=2$ and $J_{ij}=$" + f"{j_ij}")
    elif spin == 2.5:
        plt.title(state_name +r" Degeneracy for $s=\frac{5}{2}$ and $J_{ij}=$" + f"{j_ij}")
    # make os.makedirs to create the directory if it doesn't exist
    os.makedirs('data/plots', exist_ok=True)

    # Save the plot
    plt.savefig(f"data/plots/degeneracy_for_s={spin}_and_j_ij={j_ij}.pdf", format="pdf")
    # Show the plot
    #plt.show()








def plot_dict(data, spin, num_spins, dg_deg):
    output_dir = Path("data/plots")
    output_dir.mkdir(parents=True, exist_ok=True)  # Create directory if it doesn't exist

    output_path = output_dir / f"spin={spin}, num_spins={num_spins} and dg_deg={dg_deg}.jpg"
    plt.savefig(output_path, format="jpg", dpi=300)
    # Separate keys and values
    keys = list(data.keys())
    values = list(data.values())

    # Convert tuple keys to string for plotting
    keys_str = [str(key) for key in keys]

    # Create the scatter plot
    fig, ax = plt.subplots()
    ax.scatter(keys_str, values, color='blue')
    ax.set_xlabel('Spin-Pairs')
    ax.set_ylabel('Ground State Expectation Value')
    ax.set_title('Spin-Ring Ground States Expectation Values')

    # Add grid on top of the scatter plot
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

    plt.tight_layout()
    plt.savefig(output_path, format="jpg", dpi=300)
    print("Plot saved")
    #plt.show()


def save_all_plots(spin, max_spins):
    for num_spins in range(2, max_spins + 1):
        dg_deg = 0
        plot_dict(all_exp_values(spin, num_spins, dg_deg), spin, num_spins, dg_deg)


import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np


def create_connected_dots(num_dots, figsize=(8, 8)):
    """
    Creates a visualization with dots distributed on a circle and connected by dashed lines.

    Parameters:
    num_dots (int): Number of dots to distribute on the circle
    figsize (tuple): Size of the figure (width, height)
    """
    # Create a new figure with specified size
    plt.figure(figsize=figsize)

    # Create a circle
    circle = plt.Circle((0, 0), 1, fill=False, color='black')

    # Calculate dot positions
    angles = np.linspace(0, 2 * np.pi, num_dots, endpoint=False)
    x_coords = np.cos(angles)
    y_coords = np.sin(angles)

    # Add the circle to the plot
    ax = plt.gca()
    ax.add_patch(circle)

    # Plot dots
    plt.scatter(x_coords, y_coords, color='blue', s=100, zorder=2)

    # Connect all dots with dashed red lines
    for i in range(num_dots):
        for j in range(i + 1, num_dots):
            plt.plot([x_coords[i], x_coords[j]],
                     [y_coords[i], y_coords[j]],
                     'r--',
                     alpha=0.5,
                     zorder=1)

    # Set equal aspect ratio and limits
    plt.axis('equal')
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)

    # Remove axes
    plt.axis('off')
    # make os.makedirs to create the directory if it doesn't exist

    os.makedirs('data/rings', exist_ok=True)
    plt.savefig(f"data/rings/spin_ring_with_{num_dots}.pdf", format="pdf")
    plt.show()

# Example usage:
# create_connected_dots(5)  # Creates a visualization with 5 dots
# create_connected_dots(8)  # Creates a visualization with 8 dots

# Example usage:
# create_connected_dots(5)  # Creates a visualization with 5 dots

# Example usage:
# create_connected_dots(5)  # Creates a visualization with 5 dots
# create_connected_dots(8)  # Creates a visualization with 8 dots