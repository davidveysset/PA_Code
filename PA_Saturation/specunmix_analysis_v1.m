function [] = specunmix_analysis_v1(specunmix_data, mask, wavelength_set, concentrations, noise_levels, noise_strength)

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


% for each concentration(/num_iter) plot the saturation and then calculate the
% error bar for the different noise levels and for all the wavelengths.
% The error bar must be calculated using num_iterations.
% The other thing is that you must take the saturation over the mask.

% w_avg is just a number so slightly easier. same thing...




% Calculate Error bars for saturation of each noise level over the mask.
% Calculate Error bars for weighted average of each noise level.
% Plot the Unmixed on the mask to check saturation



end