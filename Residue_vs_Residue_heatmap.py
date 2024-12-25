import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("Downloads/positions.xvg", comments=["#", "@"])
n_residues = 366 #add the number of residues (aminoacids) in your protein
n_frames = data.shape[0]
data = data[:, 1:]
data = data.reshape(n_frames, n_residues, 3)
residue_displacements = np.mean(data, axis=2)
correlation_matrix = np.corrcoef(residue_displacements.T)
plt.imshow(correlation_matrix, cmap="coolwarm", vmin=-1, vmax=1)
plt.colorbar(label="Correlation Coefficient")
plt.title("Residue Correlation Map")
plt.xlabel("Residue Index")
plt.ylabel("Residue Index")
plt.show()
