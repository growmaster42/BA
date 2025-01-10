import matplotlib.pyplot as plt
from expectation_values import all_exp_values
from pathlib import Path
def plot_degeneracy():
    # Data dictionaries
    spin_one_half = {(0.5, 2): 2, (0.5, 3): 2, (0.5, 5): 2, (0.5, 7): 2}
    spin_one = {(1, 2): 2}
    spin_one_half = {(1.5, 2): 2, (1.5, 3): 2, (1.5, 5): 2}
    spin_two = {(2, 2): 2}
    spin_two_half = {(2.5, 2): 2, (2.5, 3): 2}

    # Extracting points for plotting
    x1, y1 = zip(*[(x[1], y) for x, y in spin_one_half.items()])
    x2, y2 = zip(*[(x[1], y) for x, y in spin_one.items()])
    x3, y3 = zip(*[(x[1], y) for x, y in spin_one_half.items()])
    x4, y4 = zip(*[(x[1], y) for x, y in spin_two.items()])
    x5, y5 = zip(*[(x[1], y) for x, y in spin_two_half.items()])

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.scatter(x1, y1, color='red', label='Spin 0.5', marker='o')
    plt.scatter(x2, y2, s= 200,  color='blue', label='Spin 1', marker='x')
    plt.scatter(x3, y3, color='green', label='Spin 1.5', marker='*')
    plt.scatter(x4, y4, color='purple', label='Spin 2', marker='+')
    plt.scatter(x5, y5, color='orange', label='Spin 2.5', marker='.')

    plt.xlabel('Number of Spins')
    plt.ylabel('Integer Value')
    plt.ylim(0, 3)
    plt.xlim(0, 8)
    plt.legend()
    plt.title('Spin Values in XY-Plane')
    plt.show()

#plot_degeneracy()



import matplotlib.pyplot as plt

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
        dg_deg = 1
        plot_dict(all_exp_values(spin, num_spins, dg_deg), spin, num_spins, dg_deg)


spin = 1
num_spins = 2
dg_deg = 0
plot_dict(all_exp_values(spin, num_spins, dg_deg), spin, num_spins, dg_deg)