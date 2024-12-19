% Parameters for Full Run
num_x_pixels = 200;
num_y_pixels = 200;
dx = 1e-6;
dy = 1e-6;

%Epsilon Absorption Coefficient Matrices
epsilon_770_780 = [1361 636; 1075 710]; 
epsilon_750_850 = [1405,518;691,1050];
epsilon_780_1030 = [1075 710; 1024 206];
epsilon_all = [1405,518;1361, 636;1075, 710;691,1050];

% Sets of Concentrations, Wavelengths, Noise Levels, Types
concentrations = [0.5; 0.5];
noise_levels = [1];
types = [1,2];

% FINAL SET PARAMS:

% MASTER PARAMS
concentrations = [0.5; 0.5];
expected_values = [0.5];
noise_levels = [1];
types = [1,2];


% Flag to plot the Hb plots for all the noise levels
plot_hbo = true;
%Flag to plot the HbO plots for all the noise levels
plot_hb = true;
% Flag to plot the analysis plots
plot_analysis = true;



% Returns a set of the all the reconstructed images with structure          
% (noise level, wavelength, type, concentration, x ,y)
% Also returns all the raw data with noise added with structure
% (noise level, wavelength, type, concentration, x ,t)

%{
% Run for 780_1030
E = epsilon_780_1030;
wavelengths = [780,1030];

[recon_noise_holder, noisy_sensor_data_holder] = build_pressures_w_noise(E, concentrations, wavelengths, noise_levels, types, num_x_pixels, num_y_pixels, dx,dy);
data = recon_noise_holder;
[C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error] = calc_error_saturation(data, E, wavelengths, expected_values, plot_hbo, plot_hb, plot_analysis);
data_780_1030 = {C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error};
%Returns the Concentrations calculated with non-negative least squares
% C_nnls = (noise level, concentration,type, x, y)
% Concentration of hbO: hbo_C = (noise level, concentration, x, y)
% also returns the concentration of hb similarly
% Sum of Concentration: sum_C = (noise level, concentration, x,y)
% hbo Saturation: hbo_S = (noise, level, concentration, x,y)
% Average value for circle mask: circle_concentrations = (noise_level, concentration)
% saturation error for each circle at a noise level and concentration:
% saturation_error = (noise_level, concentration)
%} 

% Run for 770_780
E = epsilon_770_780;
wavelengths = [770,780];

[recon_noise_holder, noisy_sensor_data_holder] = build_pressures_w_noise(E, concentrations, wavelengths, noise_levels, types, num_x_pixels, num_y_pixels, dx,dy);
data = recon_noise_holder;
[C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error] = calc_error_saturation(data, E, wavelengths, expected_values, plot_hbo, plot_hb, plot_analysis);
data_770_780 = {C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error};

% Run for 750_85E = epsilon_750_850;
wavelengths = [750,850];

[recon_noise_holder, noisy_sensor_data_holder] = build_pressures_w_noise(E, concentrations, wavelengths, noise_levels, types, num_x_pixels, num_y_pixels, dx,dy);
data = recon_noise_holder;
[C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error] = calc_error_saturation(data, E, wavelengths, expected_values, plot_hbo, plot_hb, plot_analysis);
data_750_850 = {C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error};


% Run for 750_770_780_850
E = epsilon_all;
wavelengths = [750,770,780,850];

[recon_noise_holder, noisy_sensor_data_holder] = build_pressures_w_noise(E, concentrations, wavelengths, noise_levels, types, num_x_pixels, num_y_pixels, dx,dy);
data = recon_noise_holder;
[C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error] = calc_error_saturation(data, E, wavelengths, expected_values, plot_hbo, plot_hb, plot_analysis);
data_all = {C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error};


% Combine the Data
total_data = {data_750_850, data_770_780, data_all};
wavelength_names = {'750-850', '770-780', '750-770-780-850'};

% All Wavelength Analysis
multi_wavelength_analysis(total_data,wavelength_names);