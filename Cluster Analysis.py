import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/rmsd_dist.xvg", comments=["#", "@"])

# Extract
counts = data[:, 0]  # First column: time (ns)
rms = data[:, 1]  # Second column: RMSD (nm)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(counts, rms, label="RMS (nm)", color="black", linewidth=1)

# Potential Plot
plt.xlabel("RMS (nm)")
plt.ylabel("Counts (x10^6)")
plt.title("Cluster Analysis using GROMOS method")
plt.legend()
plt.grid()
