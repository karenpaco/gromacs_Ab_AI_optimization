import numpy as np

# Load the .dat file
data = np.loadtxt('Downloads/ss.dat', skiprows=1, dtype=str)  # Adjust skiprows for header

with open('Downloads/ss.dat', 'r') as file:
    lines = file.readlines()[1:]  # Skip header

# Convert each row into a list of characters
data = np.array([list(line.strip()) for line in lines])

# Access rows and columns
frame_0 = data[0, :]  # First row (frame 0)
residue_0 = data[:, 0]  # First column (residue 0 across all frames)

# Analyze secondary structure frequencies
unique, counts = np.unique(data, return_counts=True)
structure_counts = dict(zip(unique, counts))
print("Secondary Structure Counts:", structure_counts)

import matplotlib.pyplot as plt

# Example: Plot secondary structure occurrence
frames = np.arange(data.shape[0])
helix_counts = (data == 'H').sum(axis=1) #H is for alpha helix. It can be changed accordingly

plt.plot(frames, helix_counts, label='Alpha Helices')
plt.xlabel('Frame')
plt.ylabel('Residue Count')
plt.title('Secondary Structure Evolution')
plt.legend()
plt.show()


#Secondary structure evolution map considering all the structures

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load the .dat file (adjust skiprows if header exists)
with open("Downloads/ss.dat", "r") as f:
    lines = f.readlines()[1:]  # Skip the header if present

# Convert secondary structure assignments to a 2D list
data = np.array([list(line.strip()) for line in lines])

# Map secondary structure codes to colors
structure_map = {
    'H': 0,  # Alpha helix
    'E': 1,  # Beta sheet
    'C': 2,  # Coil
    'T': 3,  # Turn
    'S': 4,  # Bend
    'G': 5,  # 3-10 helix
    'I': 6,  # Pi helix
    ' ': 7,  # No assignment
    '~': 8,  # Coil
    'P': 9,  # Polyproline II helix or unknown
    'B': 10, # Beta-bridge
    '=': 11  # Chain Separator
}

# Create a numerical array using the structure_map
numeric_data = np.array([[structure_map[ss] for ss in frame] for frame in data])

# Create a colormap for visualization
cmap = mcolors.ListedColormap(['black', 'gray', 'white', 'blue', 'green', 'red', 'yellow', 'pink', 'purple', 'orange', 'cyan', 'magenta'])
bounds = np.arange(len(structure_map) + 1) - 0.5
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Plot the heatmap
plt.figure(figsize=(10, 6))
plt.imshow(numeric_data.T, aspect='auto', cmap=cmap, norm=norm)

# Add labels and title
plt.xlabel('Time (ps)')
plt.ylabel('Residue Index')
plt.title('Secondary Structure Evolution')

# Define custom labels for the legend
custom_labels = [
    'Alpha Helix',      # For 'H'
    'Beta Sheet',       # For 'E'
    'Coil',             # For 'C'
    'Turn',             # For 'T'
    'Bend',             # For 'S'
    '3-10 Helix',       # For 'G'
    'Pi Helix',         # For 'I'
    'No Assignment',    # For ' '
    'Coil',        # For '~'
    'Polyproline II',   # For 'P'
    'Beta Bridge',      # For 'B'
    'Chain Separator'         # For '='
]

# Create custom legend handles
handles = [plt.Line2D([0], [0], marker='o', color='w', label=label, 
                       markerfacecolor=cmap(i), markersize=10) for i, label in enumerate(custom_labels)]

# Add the legend to the plot
plt.legend(edgecolor="Black", handles=handles, title='Secondary Structure', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    ticks=range(len(structure_map)),
    label='Secondary Structure'
).ax.set_yticklabels(list(structure_map.keys()))

# Save or display the plot
plt.tight_layout()
plt.show()


#Secondary structure evolution map considering just “β-sheet”, “coil”, and the rest as other secondary structure elements

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load the .dat file (adjust skiprows if header exists)
with open("Downloads/ss.dat", "r") as f:
    lines = f.readlines()[1:]  # Skip the header if present

# Convert secondary structure assignments to a 2D list
data = np.array([list(line.strip()) for line in lines])

# Map secondary structure codes to colors
structure_map = {
    'H': 0,  # Alpha helix
    'E': 1,  # Beta sheet
    'C': 2,  # Coil
    'T': 3,  # Turn
    'S': 4,  # Bend
    'G': 5,  # 3-10 helix
    'I': 6,  # Pi helix
    ' ': 7,  # No assignment
    '~': 8,  # Undefined
    'P': 9,  # Polyproline II helix or unknown
    'B': 10, # Beta-bridge
    '=': 11  # Potentially another undefined or system-specific state
}

# Create a numerical array using the structure_map
numeric_data = np.array([[structure_map[ss] for ss in frame] for frame in data])

# Create a colormap for visualization
cmap = mcolors.ListedColormap(['black', 'gray', 'red', 'white', 'white', 'white', 'white', 'white', 'white', 'white', 'white', 'white'])
bounds = np.arange(len(structure_map) + 1) - 0.5
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Plot the heatmap
plt.figure(figsize=(10, 6))
plt.imshow(numeric_data.T, aspect='auto', cmap=cmap, norm=norm)

# Add labels and title
plt.xlabel('Time (ps)')
plt.ylabel('Residue Index')
plt.title('Secondary Structure Evolution')

# Define custom labels for the legend (grouping some as 'Others')
custom_labels = [
    'Alpha Helix',      # For 'H'
    'Beta Sheet',       # For 'E'
    'Coil',             # For 'C'
    'Others'            # Grouping T, S, G, I, ' ', ~, P, B, =
]

# Create custom legend handles
handles = [
    plt.Line2D([0], [0], marker='o', color='w', label=custom_labels[0], 
               markerfacecolor=cmap(0), markersize=10),  # 'H'
    plt.Line2D([0], [0], marker='o', color='w', label=custom_labels[1],
               markerfacecolor=cmap(1), markersize=10),  # 'E'
    plt.Line2D([0], [0], marker='o', color='w', label=custom_labels[2],
               markerfacecolor=cmap(2), markersize=10),  # 'C'
    plt.Line2D([0], [0], marker='o', color='w', label=custom_labels[3], 
               markerfacecolor=cmap(3), markersize=10),  # 'Others'
]

# Add the legend to the plot
plt.legend(edgecolor="Black", handles=handles, title='Secondary Structure', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    ticks=range(len(structure_map)),
    label='Secondary Structure'
).ax.set_yticklabels(list(structure_map.keys()))

# Save or display the plot
plt.tight_layout()
plt.show()
