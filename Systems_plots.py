import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/temperature.xvg", comments=["#", "@"])

# Extract
time = data[:, 0]  # First column: time (ns)
temperature = data[:, 1]  # Second column: temperature (K)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(time, temperature, label="Temperature (K)", color="red", linewidth=1)

#Check average and standard deviation
mean = np.mean(temperature)
std = np.std(temperature)
plt.text(
    0.05, 0.96,  # Position (x, y) in axes fraction for the mean
    f"Mean: {mean:.2f} K",  # Text to display
    transform=plt.gca().transAxes,  # Position relative to plot axes
    fontsize=12, color="blue", bbox=dict(facecolor="white", alpha=0.5)
)
plt.text(
    0.21, 0.96,  # Position (x, y) in axes fraction for the standard deviation
    f"Std Dev: {std:.2f} K",  # Text to display
    transform=plt.gca().transAxes,  # Position relative to plot axes
    fontsize=12, color="green", bbox=dict(facecolor="white", alpha=0.5)
)

# Temperature Plot
plt.xlabel("Time (ns)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Over Time (NVT Simulation)")
plt.legend()
plt.grid()

# Show the plot
plt.tight_layout()
plt.show()

######################################################################################

import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/pressure.xvg", comments=["#", "@"])

# Extract
time = data[:, 0]  # First column: time (ns)
pressure = data[:, 1]  # Second column: pressure (bar)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(time, pressure, label="Pressure (bar)", color="red", linewidth=1)

#Check average and standard deviation
mean = np.mean(pressure)
std = np.std(pressure)
plt.text(
    0.05, 0.96,  # Position (x, y) in axes fraction for the mean
    f"Mean: {mean:.2f} bar",  # Text to display
    transform=plt.gca().transAxes,  # Position relative to plot axes
    fontsize=12, color="blue", bbox=dict(facecolor="white", alpha=0.5)
)
plt.text(
    0.21, 0.96,  # Position (x, y) in axes fraction for the standard deviation
    f"Std Dev: {std:.2f} bar",  # Text to display
    transform=plt.gca().transAxes,  # Position relative to plot axes
    fontsize=12, color="green", bbox=dict(facecolor="white", alpha=0.5)
)

# Temperature Plot
plt.xlabel("Time (ns)")
plt.ylabel("Pressure (bar)")
plt.title("Pressure Over Time (NPT Simulation)")
plt.legend()
plt.grid()

# Show the plot
plt.tight_layout()
plt.show()

######################################################################################

import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/density.xvg", comments=["#", "@"])

# Extract
time = data[:, 0]  # First column: time (ns)
density = data[:, 1]  # Second column: density (Kg/m^3)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(time, density, label="Density (Kg/m^3)", color="red", linewidth=1)

#Check average and standard deviation
mean = np.mean(density)
std = np.std(density)
plt.text(
    0.05, 0.96,  # Position (x, y) in axes fraction for the mean
    f"Mean: {mean:.2f} Kg/m^3",  # Text to display
    transform=plt.gca().transAxes,  # Position relative to plot axes
    fontsize=12, color="blue", bbox=dict(facecolor="white", alpha=0.5)
)
plt.text(
    0.28, 0.96,  # Position (x, y) in axes fraction for the standard deviation
    f"Std Dev: {std:.2f} Kg/m^3",  # Text to display
    transform=plt.gca().transAxes,  # Position relative to plot axes
    fontsize=12, color="green", bbox=dict(facecolor="white", alpha=0.5)
)

# Temperature Plot
plt.xlabel("Time (ns)")
plt.ylabel("Density (Kg/m^3)")
plt.title("Density Over Time (NPT Simulation)")
plt.legend()
plt.grid()

# Show the plot
plt.tight_layout()
plt.show()

######################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline

# Load data, skipping comments and metadata
data = np.loadtxt("Downloads/em_potential.xvg", comments=["#", "@"])

# Extract
time = data[:, 0]  # First column: time (ns)
potential = data[:, 1]  # Second column: potential (Kj/mol)

#Smoothing lines (optional)
xnew = np.linspace(time.min(), time.max(), 100)  
spl = make_interp_spline(time, potential, k=5) #the higher the K the smother, go by odd numbers
smooth = spl(xnew)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(xnew, smooth, label="Potential (KJ/mol)", color="red", linewidth=1)
#plt.plot(time, potential, label="Potential (KJ/mol)", color="red", linewidth=1) use this if you dont want the smooth algorithm

# Potential Plot
plt.xlabel("Time (ns)")
plt.ylabel("Potential (KJ/mol)")
plt.title("Energy Minimization")
plt.legend()
plt.grid()
