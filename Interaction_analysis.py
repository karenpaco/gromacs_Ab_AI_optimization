#Interaction .txt file that saves contacts but no considering COGs

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

# Load trajectory
u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")

# Select residues from two proteins
protein1 = u.select_atoms("resid 1:223 and protein")
protein2 = u.select_atoms("resid 224:366 and protein")

# Define cutoff distance (in Å)
cutoff = 4.0

# Store results
interactions = []

# Analyze frame-by-frame
for ts in u.trajectory:
    distances = distance_array(protein1.positions, protein2.positions)
    close_pairs = np.where(distances < cutoff)

    for p1, p2 in zip(*close_pairs):
        res1 = protein1[p1].resname + str(protein1[p1].resid)
        res2 = protein2[p2].resname + str(protein2[p2].resid)
        interactions.append((ts.frame, res1, res2))

# Save results to a file
with open("interacting_residues_noCOGs.txt", "w") as f:
    f.write("Frame\tResidue1\tResidue2\n")
    for frame, res1, res2 in interactions:
        f.write(f"{frame}\t{res1}\t{res2}\n")



#Generating a index file to be used in gromacs considering from the previous interaction aminoacids
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np

# Load trajectory
u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")

# Select residues from two proteins
protein1 = u.select_atoms("resid 1:223 and protein")
protein2 = u.select_atoms("resid 224:366 and protein")

# Define cutoff distance (in Å)
cutoff = 4.0

# Initialize sets to store interacting residues
interacting_residues_protein1 = set()
interacting_residues_protein2 = set()

# Analyze frame-by-frame
for ts in u.trajectory:
    distances = distance_array(protein1.positions, protein2.positions)
    close_pairs = np.where(distances < cutoff)

    # Collect interacting residues
    for p1, p2 in zip(*close_pairs):
        interacting_residues_protein1.add(protein1[p1].resid)
        interacting_residues_protein2.add(protein2[p2].resid)

# Write the NDX file
with open("interacting_residues.ndx", "w") as f:
    # Group for interacting residues in protein 1
    f.write("[ Interacting_Protein1 ]\n")
    for resid in sorted(interacting_residues_protein1):
        atoms = protein1.select_atoms(f"resid {resid}").indices + 1  # Convert to GROMACS 1-based indexing
        f.write(" ".join(map(str, atoms)) + "\n")
    
    # Group for interacting residues in protein 2
    f.write("\n[ Interacting_Protein2 ]\n")
    for resid in sorted(interacting_residues_protein2):
        atoms = protein2.select_atoms(f"resid {resid}").indices + 1  # Convert to GROMACS 1-based indexing
        f.write(" ".join(map(str, atoms)) + "\n")

print("NDX file generated: interacting_residues.ndx") #this index file will help you to analyze only the aminoacids that are interacting between each other

#Interaction .txt file that saves positions and residues considering COGs

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import distance_array

# Load the trajectory
u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")

# Select residues from two proteins
protein1 = u.select_atoms("resid 1:223 and protein")
protein2 = u.select_atoms("resid 224:366 and protein")

# Define the number of residues and cutoff distance
n_res1 = len(protein1.residues)
n_res2 = len(protein2.residues)
cutoff = 4.0  # Å

# Initialize contact map
contact_map = np.zeros((n_res1, n_res2))

# List to store interacting residues
interacting_residues = []

# Loop through trajectory frames
for ts in u.trajectory:
    # Calculate COGs for each residue in the current frame
    cog1 = np.array([residue.atoms.center_of_geometry() for residue in protein1.residues])
    cog2 = np.array([residue.atoms.center_of_geometry() for residue in protein2.residues])

    # Compute distance matrix between the two sets of COGs
    distances = distance_array(cog1, cog2)

    # Update contact map and store interacting residues for distances within the cutoff
    for i in range(n_res1):
        for j in range(n_res2):
            if distances[i, j] < cutoff:
                contact_map[i, j] += 1
                interacting_residues.append((ts.frame, protein1.residues[i].resname, protein1.residues[i].resid,
                                             protein2.residues[j].resname, protein2.residues[j].resid))

with open("interacting_residues_no_frequencies.txt", "w") as f:
    f.write("Frame\tResidue1_Name\tResidue1_ID\tResidue2_Name\tResidue2_ID\n")
    for frame, res1_name, res1_id, res2_name, res2_id in interacting_residues:
        f.write(f"{frame}\t{res1_name}\t{res1_id}\t{res2_name}\t{res2_id}\n")
    #this script does the same that the previous but it is considering COGs (centers of geometry or rach aminoacid, does not consider the atomes per each aminoacid
    #it is up to you if you choose COGs or no COGs (this one considers each atom from each aminoacid)


#Interaction .txt file that saves positions, residues, frequencies of interaction and a heatmap plot with frequencies considering COGs

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import distance_array

# Load trajectory
u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")

# Select residues from two proteins
protein1 = u.select_atoms("resid 1:223 and protein")
protein2 = u.select_atoms("resid 224:366 and protein")

# Define cutoff distance (in Å)
cutoff = 4.0

# Initialize variables
n_res1 = len(protein1.residues)
n_res2 = len(protein2.residues)
interaction_counts = np.zeros((n_res1, n_res2))

# Analyze frame-by-frame to calculate interaction frequencies
for ts in u.trajectory:
    # Calculate distance matrix for the centers of geometry (COGs) of residues
    cog1 = np.array([res.atoms.center_of_geometry() for res in protein1.residues])
    cog2 = np.array([res.atoms.center_of_geometry() for res in protein2.residues])
    
    distances = distance_array(cog1, cog2)
    contact_map = (distances < cutoff).astype(int)
    
    # Accumulate interactions
    interaction_counts += contact_map

# Normalize by the number of frames to get the frequency
interaction_frequencies = interaction_counts / len(u.trajectory)

# Save interaction frequencies to a .txt file
with open("interaction_frequencies:printed_frequencies.txt", "w") as f:
    f.write("Residue1\tResidue2\tFrequency\n")
    for i, res1 in enumerate(protein1.residues):
        for j, res2 in enumerate(protein2.residues):
            if interaction_frequencies[i, j] > 0:
                f.write(f"{res1.resname}{res1.resid}\t{res2.resname}{res2.resid}\t{interaction_frequencies[i, j]:.4f}\n")

# Print residues that interact at least once
print("Residues interacting at least once:")
with open("interacting_residues_once.txt", "w") as f:
    f.write("Residue1\tResidue2\n")
    for i, res1 in enumerate(protein1.residues):
        for j, res2 in enumerate(protein2.residues):
            if interaction_frequencies[i, j] > 0:
                print(f"{res1.resname}{res1.resid} interacts with {res2.resname}{res2.resid}")
                f.write(f"{res1.resname}{res1.resid}\t{res2.resname}{res2.resid}\n")

# Plot the interaction frequency map
plt.figure(figsize=(10, 8))
plt.imshow(interaction_frequencies, cmap="viridis", origin="lower", aspect="auto")
plt.colorbar(label="Frequency of Interaction")
plt.xlabel("Residue Index (Protein 2)")
plt.ylabel("Residue Index (Protein 1)")
plt.title("Residue Interaction Frequency Map")
plt.tight_layout()
plt.savefig("residue_interaction_frequency_map.png")
plt.show()

#Generating a index file to be used in gromacs considering COGs

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np

# Load trajectory
u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")

# Select residues from two proteins
protein1 = u.select_atoms("resid 1:223 and protein")
protein2 = u.select_atoms("resid 224:366 and protein")

# Define cutoff distance (in Å)
cutoff = 4.0

# Initialize sets to store interacting residues
interacting_residues_protein1 = set()
interacting_residues_protein2 = set()

# Analyze frame-by-frame
for ts in u.trajectory:
    # Calculate COGs for each residue
    protein1_cogs = np.array([res.atoms.center_of_geometry() for res in protein1.residues])
    protein2_cogs = np.array([res.atoms.center_of_geometry() for res in protein2.residues])

    # Compute distance matrix between COGs
    distances = distance_array(protein1_cogs, protein2_cogs)
    close_pairs = np.where(distances < cutoff)

    # Collect interacting residues
    for p1, p2 in zip(*close_pairs):
        interacting_residues_protein1.add(protein1.residues[p1].resid)
        interacting_residues_protein2.add(protein2.residues[p2].resid)

# Write the NDX file
with open("interacting_residues_cogs.ndx", "w") as f:
    # Protein 1 interacting residues group
    f.write("[ Interacting_Protein1 ]\n")
    for resid in sorted(interacting_residues_protein1):
        atoms = protein1.select_atoms(f"resid {resid}").indices + 1  # GROMACS indices start at 1
        f.write(" ".join(map(str, atoms)) + "\n")

    # Protein 2 interacting residues group
    f.write("\n[ Interacting_Protein2 ]\n")
    for resid in sorted(interacting_residues_protein2):
        atoms = protein2.select_atoms(f"resid {resid}").indices + 1  # GROMACS indices start at 1
        f.write(" ".join(map(str, atoms)) + "\n")

print("NDX file generated: interacting_residues_cogs.ndx")


#Salt birdges analysis

import MDAnalysis as mda
import numpy as np

# Load trajectory and topology
try:
    u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")
except ValueError as e:
    
    # Load full topology
    u = mda.Universe("Downloads/md_0_1.tpr")

    # Select only the atoms present in the trajectory
    subset = u.select_atoms("protein")  # Adjust the selection as needed
    u_subset = mda.Merge(subset)

    # Save the stripped topology
    stripped_topology = "subset.pdb"
    u_subset.write(stripped_topology)

    print(f"Stripped topology saved as {stripped_topology}.")
    u = mda.Universe(stripped_topology, "Downloads/md_0_1.xtc")

# Select acidic and basic residues
acidic = u.select_atoms("resname ASP GLU")
basic = u.select_atoms("resname LYS ARG HIS")

# Analyze salt bridges
cutoff = 4.0  # Å
results = []
for ts in u.trajectory:
    for acid in acidic:
        for base in basic:
            dist = np.linalg.norm(acid.position - base.position)
            if dist < cutoff:
                results.append((ts.frame, acid.resname, acid.resid, base.resname, base.resid, dist))

# Print results
print("Salt Bridges Detected:")
for frame, acid_resname, acid_resid, base_resname, base_resid, dist in results:
    print(f"Frame {frame}: {acid_resname}{acid_resid} - {base_resname}{base_resid} ({dist:.2f} Å)")


#GIF geneator of the previous heatmap plot per all the feames in the molecular dynamics

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import distance_array
import imageio
from io import BytesIO
from PIL import Image

# Load trajectory
u = mda.Universe("Downloads/md_0_1.tpr", "Downloads/md_0_1.xtc")

# Select residues from two proteins
protein1 = u.select_atoms("resid 1:223 and protein")
protein2 = u.select_atoms("resid 224:366 and protein")

# Define cutoff distance (in Å)
cutoff = 4.0

# Initialize variables
n_res1 = len(protein1.residues)
n_res2 = len(protein2.residues)
gif_frames = []  # Store frames for the GIF

# Analyze frame-by-frame and generate interaction maps
for ts in u.trajectory:
    # Calculate distance matrix for the centers of geometry (COGs) of residues
    cog1 = np.array([res.atoms.center_of_geometry() for res in protein1.residues])
    cog2 = np.array([res.atoms.center_of_geometry() for res in protein2.residues])
    
    distances = distance_array(cog1, cog2)
    
    # Create a contact map for the current frame
    contact_map = (distances < cutoff).astype(int)

    # Plot the interaction map for the current frame
    plt.figure(figsize=(10, 8))
    plt.imshow(contact_map, cmap="viridis", origin="lower", aspect="auto")
    plt.colorbar(label="Interaction (1=Yes, 0=No)")
    plt.xlabel("Residue Index (Protein 2)")
    plt.ylabel("Residue Index (Protein 1)")
    plt.title(f"Residue Interaction Map - Frame {ts.frame}")
    plt.tight_layout()
    
    # Save the plot as an image in memory
    buffer = BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    gif_frames.append(Image.open(buffer))
    plt.close()

# Create a GIF directly from memory
gif_frames[0].save(
    "residue_interaction_map.gif",
    save_all=True,
    append_images=gif_frames[1:],
    duration=100,
    loop=0
)

print("GIF saved as 'residue_interaction_map.gif'.")
#I couldnt make this work but if you can fix it it theoretically give you a heat map per frame so you wil see the evolution of aminoacid interaction per frame in the MD
