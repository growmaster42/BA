from pathlib import Path
from saving_data import *
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

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


def create_connected_dots(num_dots, figsize=(8, 8)):
    """
    Creates a visualization with dots distributed on a circle and connected by dashed lines.
    The dots are numbered starting from 0 at the northernmost point.

    Parameters:
    num_dots (int): Number of dots to distribute on the circle
    figsize (tuple): Size of the figure (width, height)
    """
    # Create a new figure with specified size
    fig, ax = plt.subplots(figsize=figsize)

    # Create a circle
    circle = plt.Circle((0, 0), 1, fill=False, color='black')

    # Calculate dot positions
    angles = np.linspace(0, 2 * np.pi, num_dots, endpoint=False)
    x_coords = np.cos(angles)
    y_coords = np.sin(angles)

    # Add the circle to the plot
    ax.add_patch(circle)

    # Plot dots and add numbers
    ax.scatter(x_coords, y_coords, color='blue', s=100, zorder=2, label='spin')

    # Add numbers next to dots
    for i in range(num_dots):
        # Calculate text position slightly outside the dot
        # Increase the radius by 0.15 to position the text outside the dot
        text_x = 1.15 * x_coords[i]
        text_y = 1.15 * y_coords[i]

        # Special case for the northernmost dot (i=0)
        if i == 0:
            ax.text(text_x, text_y, f'({i},{num_dots})',
                    ha='center', va='bottom', fontsize=22)
        else:
            ax.text(text_x, text_y, f'({i})',
                    ha='center', va='center', fontsize=22)

    # Connect all dots with dashed red lines
    for i in range(num_dots):
        for j in range(i + 1, num_dots):
            ax.plot([x_coords[i], x_coords[j]],
                    [y_coords[i], y_coords[j]],
                    'r--',
                    alpha=0.5,
                    zorder=1,
                    label='dipole dipole interaction' if (i == 0 and j == 1) else "")

    # Set equal aspect ratio and limits
    ax.set_aspect('equal')
    ax.set_xlim(-1.3, 1.3)  # Increased limits to accommodate the numbers
    ax.set_ylim(-1.3, 1.3)

    # Remove axes
    ax.axis('off')

    # Create custom legend handles
    legend_handles = [
        Line2D([0], [0], color='blue', marker='o', linestyle='None', markersize=10, label='spin'),
        Line2D([0], [0], color='red', linestyle='--', label='dipole dipole interaction'),
        Line2D([0], [0], color='black', linestyle='-', label='ring')
    ]

    # Add legend below the circle, aligned to the left
    ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=1, fontsize=22,
              frameon=True)

    # Adjust layout to reduce white space
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15)

    # Make directories if they don't exist
    os.makedirs('data/rings', exist_ok=True)

    # Save the plot
    plt.savefig(f"data/rings/spin_ring_with_{num_dots}.pdf", format="pdf", bbox_inches='tight')
    plt.show()

# Example usage




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



def plot_sys_expect_values(spin, min_spin, max_spin):
    """
    Plots three dictionaries with keys as x-values and values as y-values.
    Each dictionary gets a different color, and the points are connected by dashed splines.
    The plot is saved as a PDF in the 'data/plots' directory.
    """
    # Ensure the directory exists
    directory = "data/plots"
    os.makedirs(directory, exist_ok=True)  # Create the folder if it doesn't exist
    dict_1 = load_json('data/expectation_values_spin_rings',
                       f'spin={spin}_j_ij=-1_spins_min={min_spin}_spins_max={max_spin}.json')
    dict_2 = load_json('data/expectation_values_spin_rings',
                       f'spin={spin}_j_ij=0_spins_min={min_spin}_spins_max={max_spin}.json')
    dict_3 = load_json('data/expectation_values_spin_rings',
                       f'spin={spin}_j_ij=1_spins_min={min_spin}_spins_max={max_spin}.json')
    j_ij = [-1, 0, 1]
    # Construct the full file path
    file_path = os.path.join(directory, f"expectation_values_spin={spin}.pdf")

    # Assign unique colors for each dictionary
    colors = ['red', 'blue', 'green']
    dicts = [dict_1, dict_2, dict_3]

    # Create a figure and axis
    plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots(figsize=(10, 6))
    #create empty list for custom legend handles
    custom_legend_handles = []
    # Loop through each dictionary and plot
    for idx, data in enumerate(dicts):
        # Sort the keys for consistent plotting
        x = np.array(sorted(data.keys()))
        y = np.array([data[key] for key in x])

        # Create a spline for smooth curves
        x_new = np.linspace(x.min(), x.max(), 300)
        if len(dict_1) >= 4:
            spline = make_interp_spline(x, y, k=3) # Cubic spline
        else:
            spline = make_interp_spline(x, y, k=2)
        y_smooth = spline(x_new)

        # Plot with dashed lines and scatter points
        plt.plot(x_new, y_smooth, color=colors[idx], linestyle='--')
        plt.scatter(x, y, color=colors[idx])
        #add custom legend handles
        custom_legend_handles.append(
        Line2D([0], [0], color=colors[idx], linestyle='--', marker='o', markersize=5, label=r'$J_{ij}=$' +
                                                                                            f'{j_ij[idx]} K'))
    # Adjust axis and labels
    ax.set_xlabel(r'Number of Spins', fontsize=22)
    ax.set_ylabel(r'$\langle \hat{S}_i \cdot \hat{S}_{i+1} \rangle$', fontsize=22)

    if spin == 0.5:
        spin_number = r'$\frac{1}{2}$'
    elif spin == 1:
        spin_number = '$1$'
    elif spin == 1.5:
        spin_number = r'$\frac{3}{2}$'
    else:
        print("Invalid spin value")
    ax.set_title(r'Correlation Function for Spin Rings with $s=$' + spin_number, fontsize=30)

    # Set legend font size
    ax.legend(handles=custom_legend_handles, fontsize=30)

    # Customize tick labels (axis values) font size
    ax.tick_params(axis='both', labelsize=30)


    # Add grid
    ax.grid(True, linestyle='--', alpha=0.7)
    # Add a legend with custom handles

    plt.grid(True)
    plt.tight_layout()

    # Save the plot as PDF and show it
    plt.savefig(file_path, dpi=300, bbox_inches='tight', format='pdf', transparent=True)
    plt.show()


def plot_sys_exp_pairs(spin, num_spins, show_dict1=True, show_dict2=True, show_dict3=True):
    """
    Plot system experimental pairs with options to show/hide different dictionaries.

    Parameters:
    -----------
    spin : float
        Spin value for the plot title and filename
    num_spins : int
        Number of spins for the plot title and filename
    show_dict1 : bool, optional (default=True)
        Whether to show dictionary 1 (j_{ij}=-1)
    show_dict2 : bool, optional (default=True)
        Whether to show dictionary 2 (j_{ij}=0)
    show_dict3 : bool, optional (default=True)
        Whether to show dictionary 3 (j_{ij}=0)
    """
    # Initialize empty dictionaries
    dict_1 = load_tuple_exp(
        f"data/expectation_values_spin_rings_all_pairs/spin={spin}_j_ij=-1_spins_max={num_spins}.json")
    dict_2 = load_tuple_exp(
        f"data/expectation_values_spin_rings_all_pairs/spin={spin}_j_ij=0_spins_max={num_spins}.json")
    dict_3 = load_tuple_exp(
        f"data/expectation_values_spin_rings_all_pairs/spin={spin}_j_ij=1_spins_max={num_spins}.json")

    # Create 'data/plots' directory if it doesn't exist
    os.makedirs('data/plots', exist_ok=True)

    # Create figure and axis
    plt.figure(figsize=(10, 6))
    fig, ax = plt.subplots(figsize=(10, 6))
    # Plot each dictionary with different colors and dashed lines if show_dict is True
    if show_dict1:
        plt.plot([str(k) for k in dict_1.keys()], list(dict_1.values()),
                 'ro--', label=r"$J_{ij}=-1$ K", markersize=8)

    if show_dict2:
        plt.plot([str(k) for k in dict_2.keys()], list(dict_2.values()),
                 'bo--', label=r"$J_{ij}=0$ K", markersize=8)

    if show_dict3:
        plt.plot([str(k) for k in dict_3.keys()], list(dict_3.values()),
                 'go--', label=r"$J_{ij}=1$ K", markersize=8)

    # Customize the plot
    ax.set_xlabel(r'Spin Pairs $(i,j)$', fontsize=30)
    ax.set_ylabel(r'Correlation Value', fontsize=30)
    if spin == 0.5:
        spin_number = r'$ s= \frac{1}{2}$'
    elif spin == 1:
        spin_number = r'$s = 1$'
    ax.set_title(f'Correlators for spin pairs with ' + str(spin_number), fontsize=30)

    # Set legend font size
    ax.legend(fontsize=16)

    # Customize tick labels (axis values) font size
    ax.tick_params(axis='both', labelsize=18)

    # Rotate x-axis labels
    plt.xticks(rotation=45)

    # Add grid
    ax.grid(True, linestyle='--', alpha=0.7)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save the plot
    plt.savefig(f'data/plots/correlators_all_pairs_spin={spin}_num_spins={num_spins}.pdf')

    # Display the plot
    #plt.show()
def plot_eigenstate_evolution_one_linear(save_path='data/plots',
                                      timestamp='2025-01-24 16:00:27',
                                      user='growmaster42'):
    # Create directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)

    # Define the linear function for fitting
    def linear_func(x, m, b):
        return m * x + b

    # Add R-squared calculation function
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # Data
    data = {
        r'$J_{ij} = -1$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array([-6.275025852521177, -12.611985280603264, -14.402318738965866, -25.5146930659471,
                           -48.4529810455124])
        },
        r'$J_{ij} = 0$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array([-1.8430620905058719, -5.952907421739498, -14.711049150806287, -30.733810600289107,
                           -57.205131899154885])
        },
        r'$J_{ij} = 1$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array(
                [-7.1952301231976685, -11.046189449961695, -18.178957879591092, -36.89708294135723, -66.0228013819197])
        }
    }

    # Colors for each dataset
    colors = ['blue', 'green', 'red']

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Print fitting coefficients
    print("Fitting coefficients (m, b) for y = mx + b:")
    print("-" * 50)

    # First plot all fitted curves
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']

        # Fit the linear function
        popt, _ = curve_fit(linear_func, x, y)

        # Calculate R-squared
        y_pred = linear_func(x, *popt)
        r2 = r_squared(y, y_pred)

        # Print coefficients and R-squared
        print(f"{label}:")
        print(f"m = {popt[0]:.8f}")
        print(f"b = {popt[1]:.8f}")
        print(f"R² = {r2:.8f}")
        print("-" * 50)

        # Generate smooth curve for the fit
        x_fit = np.linspace(min(x), max(x), 200)
        y_fit = linear_func(x_fit, *popt)

        # Plot the fitted curve (dashed line)
        plt.plot(x_fit, y_fit, '--', color=color, label=r'Fitting: ' + label)

    # Then plot all data points
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']
        # Plot the original data points
        plt.plot(x, y, 'o', color=color, markerfacecolor=color,
                 markersize=8, markeredgecolor='white', markeredgewidth=1, label=label)

    # Customize the plot
    plt.xlabel('Number of Spins', fontsize=24)
    plt.ylabel(r'Eigenvalues', fontsize=24)
    plt.title(r'Eigenstate Evolution $s = 1$', fontsize=30)
    plt.legend(loc='best', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(save_path, 'eigenstate_evolution_s=1_linear.pdf'),
                bbox_inches='tight')

    # Show the plot
    plt.show()

def plot_eigenstate_evolution_one_quad(save_path='data/plots',
                                      timestamp='2025-01-24 16:00:27',
                                      user='growmaster42'):
    # Create directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)

    # Define the quadratic function for fitting
    def quad_func(x, a, b, c):
        return a * x ** 2 + b * x + c

    # Add R-squared calculation function
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # Data
    data = {
        r'$J_{ij} = -1$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array([-6.275025852521177, -12.611985280603264, -14.402318738965866, -25.5146930659471,
                           -48.4529810455124])
        },
        r'$J_{ij} = 0$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array([-1.8430620905058719, -5.952907421739498, -14.711049150806287, -30.733810600289107,
                           -57.205131899154885])
        },
        r'$J_{ij} = 1$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array(
                [-7.1952301231976685, -11.046189449961695, -18.178957879591092, -36.89708294135723, -66.0228013819197])
        }
    }

    # Colors for each dataset
    colors = ['blue', 'green', 'red']

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Print fitting coefficients
    print("Fitting coefficients (a, b, c) for y = ax² + bx + c:")
    print("-" * 50)

    # First plot all fitted curves
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']

        # Fit the quadratic function
        popt, _ = curve_fit(quad_func, x, y)

        # Calculate R-squared
        y_pred = quad_func(x, *popt)
        r2 = r_squared(y, y_pred)

        # Print coefficients and R-squared
        print(f"{label}:")
        print(f"a = {popt[0]:.8f}")
        print(f"b = {popt[1]:.8f}")
        print(f"c = {popt[2]:.8f}")
        print(f"R² = {r2:.8f}")
        print("-" * 50)

        # Generate smooth curve for the fit
        x_fit = np.linspace(min(x), max(x), 200)
        y_fit = quad_func(x_fit, *popt)

        # Plot the fitted curve (dashed line)
        plt.plot(x_fit, y_fit, '--', color=color, label=r'Fitting: ' + label)

    # Then plot all data points
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']
        # Plot the original data points
        plt.plot(x, y, 'o', color=color, markerfacecolor=color,
                 markersize=8, markeredgecolor='white', markeredgewidth=1, label=label)

    # Customize the plot
    plt.xlabel('Number of Spins', fontsize=24)
    plt.ylabel(r'Eigenvalues', fontsize=24)
    plt.title(r'Eigenstate Evolution $s = 1$', fontsize=30)
    plt.legend(loc='best', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(save_path, 'eigenstate_evolution_s=1_quad.pdf'),
                bbox_inches='tight')

    # Show the plot
    plt.show()


# Call the function
def plot_eigenstate_evolution_one_exp(save_path='data/plots',
                                      timestamp='2025-01-24 16:00:27',
                                      user='growmaster42'):
    # Create directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)

    # Define the exponential function for fitting
    def exp_func(x, a, b, c):
        return a * np.exp(b * x) + c

    # Add R-squared calculation function
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # Data
    data = {
        r'$J_{ij} = -1$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array([-6.275025852521177, -12.611985280603264, -14.402318738965866, -25.5146930659471,
                           -48.4529810455124])
        },
        r'$J_{ij} = 0$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array([-1.8430620905058719, -5.952907421739498, -14.711049150806287, -30.733810600289107,
                           -57.205131899154885])
        },
        r'$J_{ij} = 1$ K': {
            'x': np.array([3, 4, 5, 6, 7]),
            'y': np.array(
                [-7.1952301231976685, -11.046189449961695, -18.178957879591092, -36.89708294135723, -66.0228013819197])
        }
    }

    # Colors for each dataset
    colors = ['blue', 'green', 'red']

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Print fitting coefficients
    print("Fitting coefficients (a, b, c) for y = a * exp(b * x) + c:")
    print("-" * 50)

    # First plot all fitted curves
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']

        # Fit the exponential function
        popt, _ = curve_fit(exp_func, x, y)

        # Calculate R-squared
        y_pred = exp_func(x, *popt)
        r2 = r_squared(y, y_pred)

        # Print coefficients and R-squared
        print(f"{label}:")
        print(f"a = {popt[0]:.8f}")
        print(f"b = {popt[1]:.8f}")
        print(f"c = {popt[2]:.8f}")
        print(f"R² = {r2:.8f}")
        print("-" * 50)

        # Generate smooth curve for the fit
        x_fit = np.linspace(min(x), max(x), 200)
        y_fit = exp_func(x_fit, *popt)

        # Plot the fitted curve (dashed line)
        plt.plot(x_fit, y_fit, '--', color=color, label=r'Fitting: ' + label)

    # Then plot all data points
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']
        # Plot the original data points
        plt.plot(x, y, 'o', color=color, markerfacecolor=color,
                 markersize=8, markeredgecolor='white', markeredgewidth=1, label=label)

    # Customize the plot
    plt.xlabel('Number of Spins', fontsize=24)
    plt.ylabel(r'Eigenvalues', fontsize=24)
    plt.title(r'Eigenstate Evolution $s = 1$', fontsize=30)
    plt.legend(loc='best', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(save_path, 'eigenstate_evolution_s=1_exp.pdf'),
                bbox_inches='tight')

    # Show the plot
    plt.show()

def plot_eigenstate_evolution_half_linear(save_path='data/plots',
                                   timestamp='2025-02-05 12:04:13',
                                   user='growmaster42'):
    # Create directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)

    # Define the linear function for fitting
    def linear_func(x, m, b):
        return m * x + b

    # Function to calculate R-squared
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # Data
    data = {
        r'$J_{ij} = -1$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-1.6117863820993121, -4.032461143809934, -4.033541112369547, -6.569731929672131,
                           -12.118103416763434, -21.640729558870905, -35.837003990575, -55.92696708714688])
        },
        r'$J_{ij} = 0$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-0.4685528, -1.5254977, -3.6906414, -7.6981213,
                           -14.3146174, -24.4865779, -39.2991429, -59.9840949])
        },
        r'$J_{ij} = 1$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-1.8596573, -3.0596840, -4.6940646, -9.3198199,
                           -16.5563872, -27.3516828, -42.7705081, -64.0460788])
        }
    }

    # Colors for each dataset
    colors = ['blue', 'green', 'red']

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Print fitting coefficients
    print("Fitting coefficients (m, b) for y = mx + b:")
    print("-" * 50)

    # First plot all fitted curves
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']

        # Fit the linear function
        popt, _ = curve_fit(linear_func, x, y)

        # Calculate R-squared
        y_pred = linear_func(x, *popt)
        r2 = r_squared(y, y_pred)

        # Print coefficients and R-squared
        print(f"{label}:")
        print(f"m (slope) = {popt[0]:.8f}")
        print(f"b (intercept) = {popt[1]:.8f}")
        print(f"R² = {r2:.8f}")
        print("-" * 50)

        # Generate smooth curve for the fit
        x_fit = np.linspace(min(x), max(x), 200)
        y_fit = linear_func(x_fit, *popt)

        # Plot the fitted curve (dashed line)
        plt.plot(x_fit, y_fit, '--', color=color, label=r'Linear fit: ' + label)

    # Then plot all data points
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']
        # Plot the original data points
        plt.plot(x, y, 'o', color=color, markerfacecolor=color,
                 markersize=8, markeredgecolor='white', markeredgewidth=1, label=label)

    # Customize the plot
    plt.xlabel('Number of Spins', fontsize=22)
    plt.ylabel(r'Eigenvalues', fontsize=22)
    plt.title(r'Eigenstate Evolution $s = \frac{1}{2}$', fontsize=30)
    plt.legend(loc='best', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(save_path, 'eigenstate_evolution_s=0.5_linear.pdf'),
                bbox_inches='tight')

    # Show the plot
    plt.show()

def plot_eigenstate_evolution_half_quad(save_path='data/plots',
                              timestamp='2025-01-24 16:00:27',
                              user='growmaster42'):
    # Create directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)

    # Define the quadratic function for fitting
    def quad_func(x, a, b, c):
        return a * x**2 + b * x + c

    # Add R-squared calculation function
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # Data
    data = {
        r'$J_{ij} = -1$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-1.6117863820993121, -4.032461143809934, -4.033541112369547, -6.569731929672131,
                           -12.118103416763434, -21.640729558870905, -35.837003990575, -55.92696708714688])
        },
        r'$J_{ij} = 0$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-0.4685528, -1.5254977, -3.6906414, -7.6981213,
                           -14.3146174, -24.4865779, -39.2991429, -59.9840949])
        },
        r'$J_{ij} = 1$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-1.8596573, -3.0596840, -4.6940646, -9.3198199,
                           -16.5563872, -27.3516828, -42.7705081, -64.0460788])
        }
    }

    # Colors for each dataset
    colors = ['blue', 'green', 'red']

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Print fitting coefficients
    print("Fitting coefficients (a, b, c) for y = ax² + bx + c:")
    print("-" * 50)

    # First plot all fitted curves
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']

        # Fit the quadratic function
        popt, _ = curve_fit(quad_func, x, y)

        # Calculate R-squared
        y_pred = quad_func(x, *popt)
        r2 = r_squared(y, y_pred)

        # Print coefficients and R-squared
        print(f"{label}:")
        print(f"a = {popt[0]:.8f}")
        print(f"b = {popt[1]:.8f}")
        print(f"c = {popt[2]:.8f}")
        print(f"R² = {r2:.8f}")
        print("-" * 50)

        # Generate smooth curve for the fit
        x_fit = np.linspace(min(x), max(x), 200)
        y_fit = quad_func(x_fit, *popt)

        # Plot the fitted curve (dashed line)
        plt.plot(x_fit, y_fit, '--', color=color, label=r'Fitting: ' + label)

    # Then plot all data points
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']
        # Plot the original data points
        plt.plot(x, y, 'o', color=color, markerfacecolor=color,
                markersize=8, markeredgecolor='white', markeredgewidth=1, label=label)

    # Customize the plot
    plt.xlabel('Number of Spins', fontsize=22)
    plt.ylabel(r'Eigenvalues', fontsize=22)
    plt.title(r'Eigenstate Evolution $s = \frac{1}{2}$', fontsize=30)
    plt.legend(loc='best', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(save_path, 'eigenstate_evolution_s=0.5_quad.pdf'),
                bbox_inches='tight')

    # Show the plot
    plt.show()

def plot_eigenstate_evolution_half_exp(save_path='data/plots',
                              timestamp='2025-01-24 16:00:27',
                              user='growmaster42'):
    # Create directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)

    # Define the exponential function for fitting
    def exp_func(x, a, b, c):
        return a * np.exp(b * x) + c

    # Add R-squared calculation function
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # Data
    data = {
        r'$J_{ij} = -1$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-1.6117863820993121, -4.032461143809934, -4.033541112369547, -6.569731929672131,
                           -12.118103416763434, -21.640729558870905, -35.837003990575, -55.92696708714688])
        },
        r'$J_{ij} = 0$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-0.4685528, -1.5254977, -3.6906414, -7.6981213,
                           -14.3146174, -24.4865779, -39.2991429, -59.9840949])
        },
        r'$J_{ij} = 1$ K': {
            'x': np.array([3, 4, 5, 6, 7, 8, 9, 10]),
            'y': np.array([-1.8596573, -3.0596840, -4.6940646, -9.3198199,
                           -16.5563872, -27.3516828, -42.7705081, -64.0460788])
        }
    }

    # Colors for each dataset
    colors = ['blue', 'green', 'red']

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Print fitting coefficients
    print("Fitting coefficients (a, b, c) for y = a * exp(b * x) + c:")
    print("-" * 50)

    # First plot all fitted curves
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']

        # Fit the exponential function
        popt, _ = curve_fit(exp_func, x, y)

        # Calculate R-squared
        y_pred = exp_func(x, *popt)
        r2 = r_squared(y, y_pred)

        # Print coefficients
        # Print coefficients and R-squared
        print(f"{label}:")
        print(f"a = {popt[0]:.8f}")
        print(f"b = {popt[1]:.8f}")
        print(f"c = {popt[2]:.8f}")
        print(f"R² = {r2:.8f}")
        print("-" * 50)

        # Generate smooth curve for the fit
        x_fit = np.linspace(min(x), max(x), 200)
        y_fit = exp_func(x_fit, *popt)

        # Plot the fitted curve (dashed line)
        plt.plot(x_fit, y_fit, '--', color=color, label=r'Fitting: ' + label)

    # Then plot all data points
    for (label, dataset), color in zip(data.items(), colors):
        x = dataset['x']
        y = dataset['y']
        # Plot the original data points
        plt.plot(x, y, 'o', color=color, markerfacecolor=color,
                markersize=8, markeredgecolor='white', markeredgewidth=1, label=label)

    # Customize the plot
    plt.xlabel('Number of Spins', fontsize=22)
    plt.ylabel(r'Eigenvalues', fontsize=22)
    plt.title(r'Eigenstate Evolution $s = \frac{1}{2}$', fontsize=30)
    plt.legend(loc='best', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(save_path, 'eigenstate_evolution_s=0.5_exp.pdf'),
                bbox_inches='tight')

    # Show the plot
    plt.show()