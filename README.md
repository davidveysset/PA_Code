README TODO:
Look at Cote Lab Ex

<h1 align = 'center'> Photo-Acoustic Simulations of Hb and HbO Saturations in Tissue </h1>
<p align="center">
 <img src="/README_images/collage.png" width="600" height="400">
</p>


## Summary: 
This folder contains all the code necessary to run spectral unmixing for Oxygenated and Deoxygenated Hemoglobin for arbitrary tissue shape. Starting by generating arbitrary tissue shape in 2D(width, and depth), it simulates the acoustic pressure waves generated from the a laser given an input set of wavelengths. It then reconstructs the tissue from the simulated sensor data. All wave simulation is generated using k-wave. Spectral Unmixing is performed on the recovered pressures to produce maps of the Hb, HbO concentration and saturation. Additional inputs can be added such as noise and multiple iterations. This will generate error bars to characterize the performance of the reconstruction and spectral unmixing for different wavelengths. Note: these simulations can take many minutes to hours depending on how many iterations.

## Notable Features:
- Customizeable Tissue Shape with options for set shapes like circles
- Generates sensor data and reconstruction with K-Wave(faster than DAS)
- Spectral Unmixing for many wavelengths at once
- Generalizes for any wavelength, molecule type
- Characterize wavelength performance for different SNR
- Can simulate different noise levels in sensor
- Simulate arbitrary iterations of different concentrations
- Organizes results into easy to access data files
## Table of Contents:

## Installations:
### Requirements
- Matlab 2024b
- K-Wave 1.4
Need to go to K-Wave website and make an account to download. Make sure that the path to the K-Wave(k-wave-toolbox-version-1.4) is added to Matlab
- Go to Home on Matlab
- Find 'SetPath' button in toolbar
- Add the Path to the K-Wave toolbox and Save
## Getting Started:
### Define tissue dimensions
```Matlab
Ny = 200; % Number of Pixels in X-Dir
Nx = 200; % Numbers of Pixels in Y-Dir
dx = 5e-6; % Size of x-pixel
dy = 5e-6; % Size of y-pixel
```
### Create a Mask/Or Use a downloaded Mask
```Matlab
 mask = zeros(Nx,Ny);
 mask = mask + makeDisc(Nx,Ny, Nx/2, Ny/2, 4);
```

or create a tissue using:
```Matlab
mask = create_tissue_mask_v2(epsilon, circle_radius, hbo_C, c_std, fill_pct, full_Nx, full_Ny, Nx, Ny)
```
### Build a Concentration Set
First pick a set of unique concentrations you are interested in and a num of iterations. create_concentrations will tile the unique concentration sets
to allow the algorithm to calculate error bars and average over multiple runs. 
```Matlab
unique_concentrations = [0.75, 0.25; 0.5, 0.5; 0.25 0.75;];
num_iter = 10;
% Generate sets of unique concentrations with num_iter repeating sets.
concentrations = create_concentrations(unique_concentrations, num_iter);

```

### Set Noise levels
Specify a maximum noise strength and then percentages of that max strength you are interested in.
```Matlab
 noise_levels = [0,0.05,0.5];
noise_strength = 70;
```
### Run Full Simulation
Note: There are few more flags for plotting and saving but they are reasonably straightforward and included in the code documentation.
This will simulate and perform the spectral unmixing on the simulated tissue. 
```Matlab
 [specunmix_c_data, ~] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);

```

### Run Analysis
This will perform analysis by calculating and plotting the mean, standard deviation(errorbars), weighted average, and the difference from expected saturations.
```Matlab
[full_specunmix_analysis_data] = specunmix_analysis(full_specunmix_data, mask, wavelength_sets, concentrations, noise_simple, type_names, Nx ,Ny, ...
    a_save_flag, num_iter, folder_path, a_plot_flag_c, a_plot_flag_s, a_plot_flag_w_avg, a_plot_flag_diff);

```



## Examples:
Putting it all together. Here is an example of running circles of Hb and HbO over a 1mm by 1mm plane with a tissues of radius 40e-6m in the middle. 
I'm interested in seeing how the 750-850 vs 770-780 wavelengths affect the Saturation of HbO over the concentrations [0.25,0.5,0.75]. I want to generate 10 iterations for each noise level at strength of 0, 3.5, 35.

```Matlab
 %%% MASTER PARAMS
% These parameters must be fixed for all the simulations.

% Dimensions of the simulated tissue area
Ny = 200; % Number of Pixels in X-Dir
Nx = 200; % Numbers of Pixels in Y-Dir
dx = 5e-6; % Size of x-pixel
dy = 5e-6; % Size of y-pixel

%Type Names needs to be cell array
type_names = {'Hb', 'HbO'};

% Define Unique Concentrations where each column represents a type
% IMPORTANT: the name of types and order of concentrations must match!

% Specify number of iterations to measure each set of concentrations
% Number of iterations will determine how many runs will be repeated to
% calculate error bars.
unique_concentrations = [0.75, 0.25; 0.5, 0.5; 0.25 0.75;];
num_iter = 10;
% Generate sets of unique concentrations with num_iter repeating sets.
concentrations = create_concentrations(unique_concentrations, num_iter);

% Need to make your own or pick a preset mask to simulate tissue. 
% There are examples in the folder, which circles, and "tissue" that you
% can import.
 mask = zeros(200,200);
 mask = mask + makeDisc(200,200, 100, 100, 4);

%Specify the noise strength and noise levels which allow you to pick a
%percentage of the noise level. 
noise_levels = [0,0.05,0.5];
noise_strength = 70;

% IMPORTANT: Plots are not considering num_iter so only use for testing!
% Creates figures with subplots of the pressure mask for each
% concentration iteration per wavelength.
plot_flag_pressures = false;
% Creates figures of reconstructed pressures with noise, for each
% concentration and wavelength.
plot_flag_recon = false;
% Create a figure with suplots of the unmixed concentrations by type for each
% nosie level
plot_flag_concentrations = false;
%Create a figure with suplots of the unmixe saturations by type for each
% noise level
plot_flag_saturations = true;

%saves specunmix_c_data
% structure of specunmix_c_data:
% cells{num_noise, num_concentrations, 4}
% for each noise and concentration it stores:
% {sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type};)
save_flag = true;
%folder to Save the Data
folder_path = '/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/';



% Specific PARAMS
%Specify the wavelengths.
%IMPORTANT: The epsilon values, must match with the order of the
%wavelengths and the types. 


epsilon_770_780 = [1361 636; 1075 710]; 
epsilon_750_850 = [1405,518;691,1050];
epsilon_780_1030 = [1075 710; 1024 206];

%EX: Full Run for 750-850
wavelengths = [750,850];
epsilon = epsilon_750_850;
[specunmix_c_data, ~] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
    type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);

specunmix_c_data_750_850 = specunmix_c_data;

%EX: Full Run for 770-780
wavelengths = [750,850];
epsilon = epsilon_770_780;
[specunmix_c_data, ~] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
    type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);

specunmix_c_data_770_780 = specunmix_c_data;

% Combine Them:
full_specunmix_data = {specunmix_c_data_750_850, specunmix_c_data_770_780};

% IMPORTANT Need the wavelength sets to be combined in the same order.

%Analysis PARAMS

% Need to concat the different wavelengths you are analyzing together
wavelength_sets = [750,850; 770, 780];

%This is just a convenient book-keeping
noise_simple = noise_strength*noise_levels;

%NOTE: All the analysis flags are denoted with an 'a_' at the front
% This will place it in the same folder, though you could specify a
% different one, by replacing it in the input parameters.
a_save_flag = true;

%Plot flags
% For all of these flags it will plot a new figure for each type, and unique concentration 
% The plots for concentration(c),saturation(s),weighted_average(w_avg) are line plots with error bars calculated using
% num_iter. The x-axis is noise-levels and the y-axis is the measure
% quantity. 
% The plots for the difference is a bar plot. 
a_plot_flag_c = true;
a_plot_flag_s = true;
a_plot_flag_w_avg = true;
a_plot_flag_diff = true;

% Analysis function
[full_specunmix_analysis_data] = specunmix_analysis(full_specunmix_data, mask, wavelength_sets, concentrations, noise_simple, type_names, Nx ,Ny, ...
    a_save_flag, num_iter, folder_path, a_plot_flag_c, a_plot_flag_s, a_plot_flag_w_avg, a_plot_flag_diff);
```

## Acknowledgements
Special thanks to David Veysset for mentoring me through this process. Thank you to Bhaskara Rao Chintada as well for help with the delay and sum algorithm.



