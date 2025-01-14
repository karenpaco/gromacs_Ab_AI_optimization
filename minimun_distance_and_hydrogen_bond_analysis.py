import numpy as np
import matplotlib.pyplot as plt

def parse_xvg(filename):
    data = []
    with open(filename, "r") as file:
        for line in file:
            # Skip comment and header lines
            if line.startswith(("#", "@")):
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)

# Load the mindrdist.xvg file
distances = parse_xvg("Downloads/mindistnew.xvg")

# Check the shape of the data
if distances.shape[1] != 2:
    raise ValueError(f"Unexpected number of columns in mindist.xvg: {distances.shape[1]}")

# Extract time and aggregated distances
time = distances[:, 0]
distances = distances[:, 1]

# Calculate the average temperature
def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')
    
# Define a window size for the running average
window_size = 10 #adjust this to your convenience. according to gromacs 10 or 11 is a good approach
distances_avg = moving_average(distances, window_size)
time_avg = time[:len(distances_avg)]

# Plot aggregated distances over time
plt.figure(figsize=(10, 6))
plt.plot(time, distances, label="Minimum Distance")
plt.plot(time_avg, distances_avg, label=f"Average (Window: {window_size})", color="red", linewidth=2)
plt.xlabel("Time (ps)")
plt.ylabel("Distance (Å)")
plt.title("Minimum Residue-Residue Distance Over Time")
plt.legend()
plt.grid(color = 'grey', linestyle = '-', linewidth = 0.5)
plt.tight_layout()
#grid configuration
ax = plt.gca()
ax.grid(which = "both")
ax.grid(which = "major", linewidth = 1)
ax.grid(which = "minor", linewidth = 0.2)
ax.minorticks_on()

# Save the plot
plt.savefig("minimum_distance_plot.png")
plt.show()


import numpy as np
import matplotlib.pyplot as plt

def parse_xvg(filename):
    data = []
    with open(filename, "r") as file:
        for line in file:
            # Skip comment and header lines
            if line.startswith(("#", "@")):
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)

# Load the mindrdist.xvg file
hbond_distances = parse_xvg("Downloads/hbondsnew.xvg")

# Check the shape of the data
if hbond_distances.shape[1] != 2:
    raise ValueError(f"Unexpected number of columns in hbond.xvg: {hbond_distances.shape[1]}")

# Extract time and aggregated distances
time = hbond_distances[:, 0]
hbond_distances = hbond_distances[:, 1]

# Calculate the average temperature
def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')
    
# Define a window size for the running average
window_size = 10 #adjust this to your convenience. according to gromacs 10 or 11 is a good approach
hbond_avg = moving_average(hbond_distances, window_size)
time_avg = time[:len(hbond_avg)]

# Plot aggregated distances over time
plt.figure(figsize=(10, 6))
plt.plot(time, hbond_distances, label="Minimum Distance")
plt.plot(time_avg, hbond_avg, label=f"Average (Window: {window_size})", color="red", linewidth=2)
plt.xlabel("Time (ps)")
plt.ylabel("Number of Hidrogen Bonds")
plt.title("Hidrogen Bond Distance Over Time")
plt.legend()
plt.grid(color = 'grey', linestyle = '-', linewidth = 0.5)
plt.tight_layout()
#grid configuration
ax = plt.gca()
ax.grid(which = "both")
ax.grid(which = "major", linewidth = 1)
ax.grid(which = "minor", linewidth = 0.2)
ax.minorticks_on()

# Save the plot
plt.savefig("minimum_distance_plot.png")
plt.show()
