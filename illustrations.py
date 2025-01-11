import numpy as np
import matplotlib.pyplot as plt


# Function to generate magnetic dipoles within the magnet
def generate_dipoles(magnet_size, dipole_density):
    dipoles = []
    for _ in range(dipole_density):
        # Randomly place dipoles within the magnet
        x = np.random.uniform(0, magnet_size[0])
        y = np.random.uniform(0, magnet_size[1])
        # Randomize the angle of the dipole (0 to 2*pi radians)
        angle = np.random.uniform(0, 2 * np.pi)
        dipoles.append((x, y, angle))
    return dipoles


# Function to plot the magnetic dipoles
def plot_magnetic_dipoles(magnet_size, dipoles):
    fig, ax = plt.subplots()
    ax.set_xlim(0, magnet_size[0])
    ax.set_ylim(0, magnet_size[1])

    # Plot the larger magnet as a rectangle
    ax.add_patch(plt.Rectangle((0, 0), magnet_size[0], magnet_size[1], fill=True, color='lightgray', alpha=0.5))

    # Plot each dipole as an arrow
    for x, y, angle in dipoles:
        dx = np.cos(angle)  # x-component of dipole vector
        dy = np.sin(angle)  # y-component of dipole vector
        ax.arrow(x, y, dx, dy, head_width=0.2, head_length=0.3, fc='blue', ec='blue')

    ax.set_title("Magnetic Dipole Alignment within a Larger Magnet")
    ax.set_aspect('equal')
    plt.show()


# Parameters
magnet_size = (10, 10)  # Size of the larger magnet (width, height)
dipole_density = 50  # Number of magnetic dipoles

# Generate dipoles and plot them
dipoles = generate_dipoles(magnet_size, dipole_density)
plot_magnetic_dipoles(magnet_size, dipoles)
