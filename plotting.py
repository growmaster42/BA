import matplotlib.pyplot as plt
from expectation_values import all_exp_values
from pathlib import Path
import matplotlib.pyplot as plt
from intial_values import J_ij

import matplotlib.pyplot as plt

def degeneracy_plot(spin_system):
    # Define the data within the function
    data = {0.5 : {2: 1, 3: 2, 4: 1, 5: 2, 6: 1, 7: 2, 8: 1, 9: 2, 10: 1},
            1 : {2: 1, 3: 1, 4: 1, 5: 1},
            1.5: {2: 1, 3: 2, 4: 1, 5: 2},
            2: {2: 1, 3: 2, 4: 1, 5: 1},
            2.5 : {2: 1, 3: 2, 4: 1}}

    # Check if the input spin system exists in the data
    if spin_system not in data:
        print(f"Spin system {spin_system} not found in the data.")
        return

    # Create figure and axis
    fig, ax = plt.subplots()

    # Extract x and y values from the dictionary
    x_values = list(data[spin_system].keys())
    y_values = list(data[spin_system].values())

    # Create the scatter plot
    plt.scatter(x_values, y_values, color='blue')

    # Set the x-axis range from 2 to 10
    plt.xlim(1, 11)

    # Set the y-axis range from 0 to 3
    plt.ylim(0, 3)

    # Add grid for better readability
    plt.grid(True, linestyle='--', alpha=0.7)

    # Add labels and title
    plt.xlabel('Number of Spins')
    plt.ylabel('Degree of Degeneracy')
    plt.title(f'Degeneracy Plot for Spin System {spin_system}')

    # Show the plot
    plt.show()








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



