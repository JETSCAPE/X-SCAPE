import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import math
import matplotlib.colors as mcolors

# Define plot styling
mpl.rcParams['figure.figsize'] = [6., 4.5]
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['xtick.minor.width'] = 0.8
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['ytick.minor.width'] = 0.8
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['font.size'] = 15

# Define the contour levels and color map
levels = np.linspace(0.0, 100.0, 1000)
colors1 = np.array([[1, 1, 1, 1]])
colors2 = plt.cm.jet(np.linspace(0., 1, 10))
colors = np.vstack((colors1, colors2))
my_cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)

# Define the x and y grids with step 0.1
x = np.arange(-10, 10.1, 0.1)  # Step of 0.1 in x
y = np.arange(-10, 10.1, 0.1)  # Step of 0.1 in y
X, Y = np.meshgrid(x, y)

# Function to read data from the file
def read_data(filename):
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) == 0:
                continue
            T = float(parts[0])  # Time (first column)
            X_val = float(parts[1])  # X coordinate (second column)
            Y_val = float(parts[2])  # Y coordinate (third column)
            #t_val = float(parts[4])
            vx = float(parts[4])  # Temperature value (fifth column)
            vy = float(parts[5])
            s_val = float(parts[7])
            #if math.isnan(t_val):# or math.isnan(vy):
            #    print("NaN value detected in vx or vy ",vx,vy,line)       
            # Store the temperature data by time (T)
            if T not in data:
                data[T] = {'X': [], 'Y': [], 'vx': [], 'vy':[], 't':[]}
            data[T]['X'].append(X_val)
            data[T]['Y'].append(Y_val)
            data[T]['vx'].append(vx)
            data[T]['vy'].append(vy)
            data[T]['t'].append(s_val)

    return data

# Function to generate the contour plot for each time snapshot with arrows
def plot_temperature_snapshots(data, output_folder):
    for T in sorted(data.keys()):
        X_vals = np.array(data[T]['X'])
        Y_vals = np.array(data[T]['Y'])
        vx_vals = np.array(data[T]['vx'])
        vy_vals = np.array(data[T]['vy'])
        t_vals = np.array(data[T]['t'])

        # Create a grid for the data (initialized with NaNs)
        grid_data = np.full((201, 201), np.nan)  # Corrected for 201x201 grid

        # Create a grid for the velocity vectors
        vx_grid = np.full((201, 201), np.nan)
        vy_grid = np.full((201, 201), np.nan)

        # Populate the grid with magnitude and velocity components
        for i in range(len(X_vals)):
            x_idx = int(round((X_vals[i] + 10.0) * 10))  # Adjusted for 0.1 step
            y_idx = int(round((Y_vals[i] + 10.0) * 10))  # Adjusted for 0.1 step
            #grid_data[x_idx, y_idx] = math.sqrt(vx_vals[i]**2 + vy_vals[i]**2)
            grid_data[x_idx, y_idx] = t_vals[i]
            vx_grid[x_idx, y_idx] = vx_vals[i]
            vy_grid[x_idx, y_idx] = vy_vals[i]

        # Create a contour plot
        fig, ax = plt.subplots(figsize=(6, 4.5))
        cont = ax.contourf(X, Y, grid_data.transpose(), levels, cmap=my_cmap, extend='both')
        cbar = fig.colorbar(cont)
        cbar.set_label(r"Entropy Density", fontsize=15)
        time_text = ax.text(-12, 11, r"$\tau = {0:4.2f}$ fm/c".format(T), fontsize=15)

        # Add quiver plot for the vector field
        skip = 10  # Skip every 10 grid points to avoid clutter
        X_quiver, Y_quiver = X[::skip, ::skip], Y[::skip, ::skip]
        vx_quiver = vx_grid[::skip, ::skip].transpose()
        vy_quiver = vy_grid[::skip, ::skip].transpose()

        ax.quiver(
            X_quiver, Y_quiver, vx_quiver, vy_quiver, 
            angles='xy', scale_units='xy', scale=1, 
            color='black', alpha=0.7
        )

        # Set axis labels and limits
        ax.set_xlabel(r"$x$ (fm)", fontsize=15)
        ax.set_ylabel(r"$y$ (fm)", fontsize=15)
        ax.set_xlim([-10, 10])
        ax.set_ylim([-10, 10])
        
        plt.tight_layout()
        fig.suptitle('OnFly Profile 50-60 % TRENTO', fontsize=10)

        # Save the figure
        output_filename = f"{output_folder}/OnFly_Contour_XY_{int(T*10)}.png"  # Multiply T by 10 for better filename
        plt.savefig(output_filename)
        plt.close(fig)
        print(f"Saved snapshot for time {T:.2f} fm/c")

# Main function to execute the plotting 
def main():
    input_file = '/Users/ritobandatta/Desktop/X-SCAPE/ResultFolder/XScapeMusicProfileConcurrentrun0.txt'  # Input file containing the data
    output_folder = '/Users/ritobandatta/Desktop/X-SCAPE/ResultFolder'  # Folder to save the plots

    # Read the data from the filef
    data = read_data(input_file)

    # Ensure the output folder exists
    import os
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Generate and save the plots
    plot_temperature_snapshots(data, output_folder)

if __name__ == "__main__":
    main()