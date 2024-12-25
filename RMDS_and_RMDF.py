import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/rmsd.xvg", comments=["#", "@"])

# Extract
time = data[:, 0]  # First column: time (ns)
rmsd = data[:, 1]  # Second column: RMSD (nm)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(time, rmsd, label="rmsd (nm)", color="black", linewidth=1)

# Potential Plot
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("RMSD of the Backnone to Backbone")
plt.legend()
plt.grid()

######################################################################

import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/rmsf.xvg", comments=["#", "@"])

# Extract
atom = data[:, 0]  # First column: atom
rmsf = data[:, 1]  # Second column: nm

# Plot
plt.figure(figsize=(10, 8))
plt.plot(atom, rmsf, label="nm", color="black", linewidth=1)

# Potential Plot
plt.xlabel("Atom")
plt.ylabel("nm")
plt.title("RMSF Fluctuation")
plt.legend()
plt.grid()
