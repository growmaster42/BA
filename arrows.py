import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def generate_arrow_pdfs(output_file1="alternating_arrows.pdf", output_file2="all_up_arrows.pdf"):
    """
    Generate two compact PDF files with centered arrow patterns:
    1. A pattern of 4 rows with 5 arrows each, alternating up (red) and down (blue)
    2. A pattern of 4 rows with 5 arrows each, all pointing up (red)

    The arrow patterns are precisely centered in each PDF with minimal margins.

    Parameters:
    -----------
    output_file1 : str
        Filename for the alternating arrow pattern PDF
    output_file2 : str
        Filename for the all-up arrow pattern PDF

    Returns:
    --------
    tuple
        Paths to the two generated PDF files
    """

    # Function to create arrow vertices
    def create_arrow_up(x, y, size=1, width_factor=0.25):
        """Create an up arrow with narrow width"""
        width = size * width_factor
        return np.array([
            [x - width, y],  # base left
            [x - width, y + size * 0.6],  # shaft left
            [x - size * 0.5, y + size * 0.6],  # left wing bottom
            [x, y + size * 1.5],  # tip
            [x + size * 0.5, y + size * 0.6],  # right wing bottom
            [x + width, y + size * 0.6],  # shaft right
            [x + width, y]  # base right
        ])

    def create_arrow_down(x, y, size=1, width_factor=0.25):
        """Create a down arrow with narrow width"""
        width = size * width_factor
        return np.array([
            [x - width, y],  # base left
            [x - width, y - size * 0.6],  # shaft left
            [x - size * 0.5, y - size * 0.6],  # left wing top
            [x, y - size * 1.5],  # tip
            [x + size * 0.5, y - size * 0.6],  # right wing top
            [x + width, y - size * 0.6],  # shaft right
            [x + width, y]  # base right
        ])

    # Common settings
    arrow_size = 0.5
    x_spacing = 1.0
    y_spacing = 1.2

    # Calculate dimensions for proper centering
    pattern_width = 5 * x_spacing
    pattern_height = 3 * y_spacing + arrow_size * 1.5

    # Generate first PDF - Alternating arrows
    with PdfPages(output_file1) as pdf:
        # Create figure with minimal margins
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_axes([0, 0, 1, 1])  # Use full figure area
        ax.set_aspect('equal')

        # Calculate center positions
        center_x = 2.5  # Center of the 5-inch figure
        center_y = 1.5  # Center of the 3-inch figure

        # Starting positions to center the pattern
        start_x = center_x - pattern_width / 2 + x_spacing / 2
        start_y = center_y + pattern_height / 2 - arrow_size * 1.5

        # Draw 4 rows of 5 arrows, alternating up and down
        for row in range(4):
            base_y = start_y - row * y_spacing

            for col in range(5):
                x_pos = start_x + col * x_spacing

                if col % 2 == 0:  # Even columns: up arrows (red)
                    arrow = create_arrow_up(x_pos, base_y, arrow_size)
                    ax.add_patch(Polygon(arrow, closed=True, facecolor='red', edgecolor='black', linewidth=0.5))
                else:  # Odd columns: down arrows (blue)
                    arrow = create_arrow_down(x_pos, base_y + arrow_size * 1.5, arrow_size)
                    ax.add_patch(Polygon(arrow, closed=True, facecolor='blue', edgecolor='black', linewidth=0.5))

        # Set limits to properly frame the centered content
        margin = 0.3  # Small margin around the pattern
        ax.set_xlim(center_x - pattern_width / 2 - margin, center_x + pattern_width / 2 + margin)
        ax.set_ylim(center_y - pattern_height / 2 - margin, center_y + pattern_height / 2 + margin)
        ax.axis('off')  # No axes

        pdf.savefig(fig, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)

    # Generate second PDF - All up arrows
    with PdfPages(output_file2) as pdf:
        # Create figure with minimal margins
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_axes([0, 0, 1, 1])  # Use full figure area
        ax.set_aspect('equal')

        # Draw 4 rows of 5 arrows, all pointing up
        for row in range(4):
            base_y = start_y - row * y_spacing

            for col in range(5):
                x_pos = start_x + col * x_spacing
                arrow = create_arrow_up(x_pos, base_y, arrow_size)
                ax.add_patch(Polygon(arrow, closed=True, facecolor='red', edgecolor='black', linewidth=0.5))

        # Set limits to properly frame the centered content
        ax.set_xlim(center_x - pattern_width / 2 - margin, center_x + pattern_width / 2 + margin)
        ax.set_ylim(center_y - pattern_height / 2 - margin, center_y + pattern_height / 2 + margin)
        ax.axis('off')  # No axes

        pdf.savefig(fig, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)

    return output_file1, output_file2


# Example usage
if __name__ == "__main__":
    file1, file2 = generate_arrow_pdfs()
    print(f"Generated PDFs: {file1} and {file2}")