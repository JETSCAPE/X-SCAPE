import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import math
import matplotlib.colors as mcolors

# Define plot styling
mpl.rcParams['figure.figsize'] = [10.0, 8.0]  # Updated to larger plot size
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
levels = np.linspace(0, 10.0, 100)
colors1 = np.array([[1, 1, 1, 1]])
colors2 = plt.cm.jet(np.linspace(0., 1, 10))
colors = np.vstack((colors1, colors2))
my_cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)

# Define the x and y grids with step 0.1
x = np.arange(-15.0, 15.0, 0.3)  # Step of 0.1 in x
y = np.arange(-15.0, 15.0, 0.3)  # Step of 0.1 in y
X, Y = np.meshgrid(x, y)

# Function to read data from the file
def read_data(filename):
    data={}
    data[0] = {'X': [], 'Y': [], 'S': []}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            #print(parts)
            #if len(parts) == 0 or len(parts)>3:
            #    continue
            X = float(parts[0])  # Time (first column)
            Y = float(parts[1])  # X coordinate (second column)
            #Z = float(parts[2])  # Y coordinate (third column)
            #t_val = float(parts[4])
            S = float(parts[3])  # Temperature value (fifth column)
            #vy = float(parts[5])
            #if math.isnan(t_val):# or math.isnan(vy):
            #    print("NaN value detected in vx or vy ",vx,vy,line)       
            # Store the temperature data by time (T)
            T = 0
            #if Z !=0 or X>15 or X<-5 or Y>5 or Y<-5:
            #    continue

            data[T]['X'].append(X)
            data[T]['Y'].append(Y)
            data[T]['S'].append(S)
            #data[T]['t'].append(t_val)

    return data

# Function to generate the contour plot for each time snapshot
def plot_temperature_snapshots(data, output_folder):
    for T in sorted(data.keys()):
        X_vals = np.array(data[T]['X'])
        Y_vals = np.array(data[T]['Y'])
        #t_vals = np.array(data[T]['t'])
        S_vals = np.array(data[T]['S'])
        #vy_vals = np.array(data[T]['vy'])
        #print(X_vals)
        #print(Y_vals)
        # Create a grid for the data (initialized with NaNs)
        grid_data = np.full((100, 100), 0)  # Corrected for 201x201 grid

        # Populate the grid with temperature values
        for i in range(len(X_vals)):
            x_idx = round((X_vals[i] + 15.0)/0.3)  # Adjusted for 0.1 step
            y_idx = round((Y_vals[i] + 15.0)/0.3)  # Adjusted for 0.1 step
            #grid_data[x_idx,y_idx] = t_vals[i]
            #print(X_vals[i],x_idx,Y_vals[i],y_idx,S_vals[i])
            grid_data[x_idx, y_idx] = S_vals[i]
            #if X_vals[i] == -3.0:
                #print(X_vals[i], Y_vals[i], x_idx, y_idx, vx_vals[i], vy_vals[i],grid_data[x_idx, y_idx])

        # Create a contour plot
        #print(grid_data)
        #norm = mcolors.LogNorm(vmin=max(np.nanmin(grid_data), 1e-3), vmax=np.nanmax(grid_data))
        fig, ax = plt.subplots(figsize=(8, 5.5))
        cont = ax.contourf(X, Y, grid_data.transpose(), levels, cmap=my_cmap, extend='both')
        cbar = fig.colorbar(cont)
        cbar.set_label(r"Entropy Density ", fontsize=15)
        time_text = ax.text(-3, 7.5 , r"$\tau = {0:4.2f}$ fm/c".format(T), fontsize=12)



        # Set axis labels
        ax.set_xlabel(r"$x$ (fm)", fontsize=15)
        ax.set_ylabel(r"$y$ (fm)", fontsize=15)
        #ax.set_xlim([-15, 15])
        #ax.set_ylim([-15, 15])
        
        plt.tight_layout()
        fig.suptitle('On the Fly Profile 50-60% TRENTO', fontsize=10)
        #fig.subplots_adjust(top=0.002) 
        # Save the figure
        #output_filename = f"{output_folder}/InitialStage_Contour_XY_{int(T*10)}.png"  # Multiply T by 10 for better filename
        plt.savefig("/Users/ritobandatta/Desktop/X-SCAPE/ResultFolder/OnFly_TR.png")
        plt.close(fig)
        print(f"Saved snapshot for time {T:.2f} fm/c")

# Main function to execute the plotting
def main():
    input_file = '/Users/ritobandatta/Desktop/X-SCAPE/ResultFolder/XScapeTrentoProfileConcurrentrun0.txt'  # Input file containing the data
    output_folder = 'readfromtrento'  # Folder to save the plots

    # Read the data from the files
    data = read_data(input_file)

    # Ensure the output folder exists
    import os
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Generate and save the plots
    plot_temperature_snapshots(data, output_folder)

if __name__ == "__main__":
    main()
