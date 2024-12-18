README TODO:
Look at Cote Lab Ex

<h1 align = 'center'> Photo-Acoustic Simulations of Hb and HbO Saturations in Tissue </h1>
Picture: 

## Summary: 
This folder contains all the code necessary to run spectral unmixing for Oxygenated and Deoxygenated Hemoglobin for arbitrary tissue shape. Starting by generating arbitrary tissue shape in 2D(width, and depth), it simulates the acoustic pressure waves generated from the a laser given an input set of wavelengths. It then reconstructs the tissue from the simulated sensor data. All wave simulation is generated using k-wave. Spectral Unmixing is performed on the recovered pressures to produce maps of the Hb, HbO concentration and saturation. Additional inputs can be added such as noise and multiple iterations. This will generate error bars to characterize the performance of the reconstruction and spectral unmixing for different wavelengths. Note: these simulations can take many minutes to hours depending on how many iterations.

## Notable Features:
- Customizeable Tissue Shape with options for set shapes like circles
- Generates sensor data and reconstruction with K-Wave(faster than DAS)
- Spectral Unmixing for many wavelengths at once
- Characterize wavelength performance for different SNR
- Organizes results into easy to access data files
## Table of Contents:

## Installations:
# Requirements
- Matlab 2024b
- K-Wave 
## Getting Started:

## Features:

## Examples:

## Acknowledgements


