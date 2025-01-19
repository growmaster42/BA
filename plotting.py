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
    x_values = [x for x in data[(spin, j_ij)].keys() if x != 2]
    y_values = [data[(spin, j_ij)][x] for x in x_values]

    # Create the scatter plot
    plt.scatter(x_values, y_values, color='blue')

    # Set the x-axis range from 2 to 10
    plt.xlim(2, max(data[(spin, j_ij)].keys()) + 1)

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
        plt.title(state_name + r" Degeneracy for $s=1$ and $J_{ij}=$" + f"{j_ij}")
    elif spin == 1.5:
        plt.title(state_name + r" Degeneracy for $s=\frac{3}{2}$ and $J_{ij}=$" + f"{j_ij}")
    elif spin == 2:
        plt.title(state_name + r" Degeneracy for $s=2$ and $J_{ij}=$" + f"{j_ij}")
    elif spin == 2.5:
        plt.title(state_name +r" Degeneracy for $s=\frac{5}{2}$ and $J_{ij}=$" + f"{j_ij}")
    # make os.makedirs to create the directory if it doesn't exist
    os.makedirs('data/plots', exist_ok=True)

    # Save the plot
    plt.savefig(f"data/plots/degeneracy_for_s={spin}_and_j_ij={j_ij}_in_state={state}.pdf", format="pdf")
    # Show the plot
    #plt.show()


import matplotlib.pyplot as plt
import numpy as np


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


import matplotlib.pyplot as plt
import numpy as np

def plot_vectors():
    fig, ax = plt.subplots(figsize=(10, 6))

    # Parameters
    n_vectors = 30
    vector_length = 1
    spacing = 1

    # First group: 30 evenly distributed aligned vectors pointing upwards
    x_start_left = np.zeros(n_vectors)
    y_start_left = np.linspace(0, (n_vectors - 1) * spacing, n_vectors)
    u_left = np.zeros(n_vectors)
    v_left = np.ones(n_vectors) * vector_length

    # Second group: 30 evenly distributed antiparallel aligned vectors
    x_start_right = np.ones(n_vectors) * 5
    y_start_right = np.linspace(0, (n_vectors - 1) * spacing, n_vectors)
    u_right = np.zeros(n_vectors)
    v_right = np.ones(n_vectors) * vector_length * (-1) ** np.arange(n_vectors)

    # Plotting the first group
    for i in range(n_vectors):
        color = 'red' if v_left[i] > 0 else 'blue'
        ax.quiver(x_start_left[i], y_start_left[i], u_left[i], v_left[i], angles='xy', scale_units='xy', scale=1, color=color)

    # Plotting the second group
    for i in range(n_vectors):
        color = 'red' if v_right[i] > 0 else 'blue'
        ax.quiver(x_start_right[i], y_start_right[i], u_right[i], v_right[i], angles='xy', scale_units='xy', scale=1, color=color)

    # Setting limits and labels
    ax.set_xlim(-1, 6)
    ax.set_ylim(-1, (n_vectors) * spacing)
    ax.set_aspect('equal')
    ax.set_title('Vector Alignment')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')

    # Save the figure as PDF
    plt.savefig('vector_alignment.pdf')
    plt.close()

# Call the function to create and save the plot

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import matplotlib.pyplot as plt

def generate_coupling_graph(coupling_type):
    """
    Generate a PDF illustrating ferromagnetic or antiferromagnetic coupling.

    Args:
        coupling_type (str): "ferro" for ferromagnetic coupling or "antiferro" for antiferromagnetic coupling.
    """
    if coupling_type not in ["ferro", "antiferro"]:
        raise ValueError("Invalid coupling type. Choose 'ferro' or 'antiferro'.")

    rows, cols = 5, 20  # Number of rows and columns
    fig, ax = plt.subplots(figsize=(12, 4))

    arrow_length = 0.5
    for i in range(rows):
        for j in range(cols):
            if coupling_type == "ferro":
                direction = 1  # All arrows point up
                color = "red"
            elif coupling_type == "antiferro":
                if i == 0:
                    direction = 1  # First row all up
                    color = "red"
                else:
                    direction = 1 if j % 2 == 0 else -1  # Alternating directions row-wise
                    color = "red" if direction == 1 else "blue"

            # Draw arrows aligned perfectly
            start_y = -i + (arrow_length / 2 if direction == -1 else -arrow_length / 2)
            ax.arrow(j, start_y, 0, arrow_length * direction, head_width=0.2, head_length=0.2, fc=color, ec=color, lw=0.5)

    # Set axis limits and remove the axes for better visual appeal
    ax.set_xlim(-1, cols)
    ax.set_ylim(-rows, 1)
    ax.axis("off")

    # Save the figure
    filename = "ferromagnetic.pdf" if coupling_type == "ferro" else "antiferromagnetic.pdf"
    plt.savefig(filename, bbox_inches="tight")
    plt.close()
    print(f"PDF saved as {filename}")

# Example usage
# generate_coupling_graph("ferro")
# generate_coupling_graph("antiferro")


