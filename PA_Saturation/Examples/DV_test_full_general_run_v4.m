%%% MASTER PARAMS
% These parameters must be fixed for all the simulations.

% Dimensions of the simulated tissue area
Ny = 200; % Number of Pixels in X-Dir
Nx = 200; % Numbers of Pixels in Y-Dir
dx = 5e-6; % Size ofx-pixel
dy = 5e-6; % Size of y-pixel

%Type Names needs to be cell array
type_names = {'Hb', 'HbO'};

% Define Unique Concentrations where each column represents a type
% IMPORTANT: the name of types and order of concentrations must match!

% Specify number of iterations to measure each set of concentrations
% Number of iterations will determine how many runs will be repeated to
% calculate error bars.f 
unique_concentrations = [0.35 0.65];
num_iter = 30; 
% Generate sets of unique concentrations with num_iter repeating sets.
concentrations = create_concentrations(unique_concentrations, num_iter);

% Need to make your own or pick a preset mask to simulate tissue. 
% There are examples in the folder, which circles, and "tissue" that you
% can import.
 mask = zeros(200,200);
 mask = mask + makeDisc(200,200, 100, 100, 4);

%Specify the noise strength and noise levels which allow you to pick a
%percentage of the noise level. 
noise_levels = logspace(-1,0,6);%[0,0.5,1]; %always have first level at 0 to compute SNR
noise_levels = [0, noise_levels];
noise_strength = 80;

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
plot_flag_saturations = false;

%saves specunmix_c_data
% structure of specunmix_c_data:
% cells{num_noise, num_concentrations, 4}
% for each noise and concentration it stores:
% {sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type};)
save_flag = true;
%folder to Save the Data
%folder_path = '/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/';
folder_path = '/Users/davidveysset/Library/CloudStorage/OneDrive-MassGeneralBrigham/A_PROJECTS/03_PAM/Matlab/PAT/Simulations results/';


% Specific PARAMS
%Specify the wavelengths.
%IMPORTANT: The epsilon values, must match with the order of the
%wavelengths and the types. 

epsilon_770_780 = [1361 636; 1075 710]; 
epsilon_750_850 = [1405,518; 691,1050];
epsilon_780_1030 = [1075 710; 1024 206];

%EX: Full Run for 750-850
wavelengths = [750,850];
epsilon = epsilon_750_850;
[specunmix_c_data, noisy_sensor_data] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
    type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);

specunmix_c_data_750_850 = specunmix_c_data;
noisy_sensor_data_750_850 = noisy_sensor_data;
snr_array_750_850 = snr_compute(noisy_sensor_data, wavelengths, Nx, noise_levels, noise_strength);


% %EX: Full Run for 770-780
% wavelengths = [750,850];
% epsilon = epsilon_770_780;
% [specunmix_c_data, ~] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
%     type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);
% 
% specunmix_c_data_770_780 = specunmix_c_data;

%EX: Full Run for 780-1030
wavelengths = [780,1030];
epsilon = epsilon_780_1030;
[specunmix_c_data, noisy_sensor_data] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
    type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag,folder_path);

specunmix_c_data_780_1030 = specunmix_c_data;
noisy_sensor_data_780_1030 = noisy_sensor_data;
snr_array_780_1030 = snr_compute(noisy_sensor_data, wavelengths, Nx, noise_levels, noise_strength);

% Combine Them:
%full_specunmix_data = {specunmix_c_data_750_850, specunmix_c_data_770_780};
full_specunmix_data = {specunmix_c_data_750_850, specunmix_c_data_780_1030};
full_snr_data = {snr_array_750_850,snr_array_780_1030};

% IMPORTANT Need the wavelength sets to be combined in the same order.

%Analysis PARAMS

% Need to concat the different wavelengths you are analyzing together
wavelength_sets = [750, 850; 780, 1030];

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
a_plot_flag_c = false;
a_plot_flag_s = true;
a_plot_flag_w_avg = false;
a_plot_flag_diff = false;

% Analysis function
 mask_analysis = zeros(200,200);
 mask_analysis = mask_analysis + makeDisc(200,200, 100, 100, 3);
[full_specunmix_analysis_data, saturation_mean_holder, saturation_std_holder] = specunmix_analysis_saturation_specific(full_specunmix_data, mask_analysis, wavelength_sets, concentrations, noise_simple, type_names, Nx ,Ny, ...
    a_save_flag, num_iter, folder_path, a_plot_flag_c, a_plot_flag_s, a_plot_flag_w_avg, a_plot_flag_diff);



num_wavelength_sets = length(full_specunmix_data);
num_type = size(concentrations,2);
num_concentration = size(concentrations,1); 
num_c_unique = num_concentration/num_iter;
colors = lines(num_wavelength_sets); % Generate distinct colors for each wavelength set
num_noise = length(noise_levels);


if a_plot_flag_s
for t = 1:num_type
    for c = 1:num_c_unique
        figure;
        hold on;

        for w = 1:num_wavelength_sets
            sat_mean = squeeze(saturation_mean_holder(w,:,c,t));
            sat_std = squeeze(saturation_std_holder(w,:,c,t));
            %h = plot( sat_mean, '-', 'LineWidth', 2, 'Color', colors(w, :));
            snr_x_axis = full_snr_data{w};
            errorbar(snr_x_axis,sat_mean, sat_std, 'o-', 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', colors(w, :), 'DisplayName', num2str(wavelength_sets(w, :)));
        end

        %title(sprintf('%s Saturation Averages for c = %d', type_names{t}, c), 'Interpreter', 'none');
        legend('show', 'Location', 'best', 'FontSize', 14);
        xlabel('SNR');
        ylabel('Saturation Mean');
        %xlim([0.8, 1200])
        %xticks(logspace(0, log10(1200), 10));
        set(gca, 'FontSize',14, 'XScale', 'log')
        %xticks(1:num_noise);
        %xticklabels(arrayfun(@num2str, noise_simple, 'UniformOutput', false));
        grid on;
    end
end
end
