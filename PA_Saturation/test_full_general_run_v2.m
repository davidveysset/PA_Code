

%%% MASTER PARAMS
% Dimensions of the simulated tissue area
Ny = 200;
Nx = 200;
dx = 5e-6;
dy = 5e-6;

%Concentration Names needs to be cell array
type_names = {'Hb', 'HbO'};
% Concentratiosn where each column specifies the type
% and each row species an iteration of set concentrations to run
concentrations = [0.75, 0.25; 0.5, 0.5; 0.25 0.75;];
% Need to make your own or pick a preset mask to simulate tissue. 
 mask = zeros(200,200);
 mask = mask + makeDisc(200,200, 100, 100, 4);

%Epsilon Matrix, can also for Hb, HbO use the look up table
epsilon = [1361 636; 1075 710];

%Specify the noise strength and noise levels which allow you to pick a
%percentage of the noise level
noise_levels = [0,0.05,0.5];
noise_strength = 70;

% Will create a figure with subplots of the pressure mask for each
% concentration iteration per wavelength.
plot_flag_pressures = true;

plot_flag_recon = true;
% Create a figure with suplots of the unmixed concentrations by type for each
% nosie level
plot_flag_concentrations = true;
%Create a figure with suplots of the unmixe saturations by type for each
% noise level
plot_flag_saturations = true;

save_flag = true;

%Folder to Save the Data
folder_path = '/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/';

wavelengths = [770,780];

[specunmix_c_data, noisy_sensor_data_holder] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
    type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);
% analysis functions:

% input takes a set of specunmix_c_data from the different wavelengths = 
%{specunmix_c_data_770_780, specunmix_c_data_750_850}
% also takes a set of wavelengths = [wavelengths_770_780, wavelengths_750_850]
% also takes concentrations as the expected values.
% also takes noises to calculate the Signal Strength

% each specunmix_c_data is
% cell array of (num_noise, num_concentration(iter))
% each cell has an array of {sum_C, saturations_by_type,
% concentrations_by_type, w_avg_by_type}
% each of them is (num_type, Nx, Ny)




% Calculate Error bars for saturation of each noise level over the mask.
% Calculate Error bars for weighted average of each noise level.
% Plot the Unmixed on the mask to check saturation

