#  MATLAB FDTD Simulations

##  Introduction
This repository contains two distinct **Finite-Difference Time-Domain (FDTD)** simulations written in MATLAB. These scripts simulate electromagnetic wave interactions, demonstrating both the **Doppler Effect** with moving objects and **Electromagnetic Scattering** by static objects.


### 1. `static_sphere.m`: Realistic Doppler Effect
This script simulates a dielectric object moving back and forth within a 2D grid to demonstrate the Doppler frequency shift.


* **Key Features:**
    * **Moving Object:** A dielectric block ($\epsilon_r = 6.0$) that bounces between boundaries.
    * **Source:** A 1 GHz Gaussian pulse source.
    * **Probes:** Measures electric fields in front (approaching) and behind (receding) the object.
    * **Analysis:** Performs an FFT (Fast Fourier Transform) at the end to calculate and display the precise frequency shift in GHz .
* **Output:**
    * Real-time animation of the $E_z$ field and Magnetic Intensity $|H|$.
    * Final plots showing Time-domain signals and Frequency-domain shifts.

### 2. `static_FDTD_2D.m`: Scattering by Steel Cylinder
This script models the scattering of a plane wave by a highly conductive steel cylinder using advanced FDTD techniques.


* **Key Features:**
    * **TF/SF Formalism:** Uses Total-Field/Scattered-Field formulation to inject plane waves cleanly.
    * **UPML Boundaries:** Implements Uniaxial Perfectly Matched Layers (UPML) to absorb outgoing waves and prevent reflections.
    * **Material:** Simulates a steel cylinder with high conductivity ($\sigma = 10^7$).
* **Output:**
    * Real-time heatmaps of $E_z$ (Electric Field) and $H_z$ (Magnetic Field) .
    * Intensity vs. Time graphs at four specific probe points.

##  Requirements
* **Software:** MATLAB (No specific toolboxes required, uses standard math functions).
