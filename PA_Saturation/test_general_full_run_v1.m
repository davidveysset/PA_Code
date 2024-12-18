% Structure
% two wavelengths, two types hb, hbo
% construct the maps of the 9 points for hb and inverse for hbo
% Input also needs to consider the epsilon values, concentration_map,
% num_x_steps, num_y_steps


% The strucutre has to be for each epsilon (N wavelength), all
% concentrations, all types, all_noise generate a set cell array:
% {saturation, concentrations}  
% - saturation = size(num_noise, num_wavelength, num_concentration, Nx,Ny)
% - concentration = size(num_noise, num_wavelength, num_concentration, Nx,Ny)

% Then you run it for all the epsilons you are interested in ... multiple
% Es
% Then you run the cell array of cell the 2 cell arrays that does the
% analysis

epsilon_770_780 = [1361 636; 1075 710];
epsilon_750_850 = [1405,518; 691,1050];
epsilon_780_1030 = [1075 710; 1024 206];
epsilon_all = [1405,518;1361, 636;1075, 710;691,1050];
wavelengths = {epsilon_770_780};
num_wavelength = length(wavelengths);


noises = [0,0.1,0.2,0.3,0.4];

%concentration is an num_concentrations(the iterations) x num_type
% ex: []
concentrations = [a b c ; d e f]

noise_strength = 100;
num_noise = length(noises);
types = [1,2];
num_types = length(types);

Nx = 400;
Ny = 200;
dx = 5e-6;
dy = 5e-6;



hold on;

mask = zeros(Nx,Ny);

%Pressure maks is of size(num_concentrations, num_wavelength, # x pixels, # y pixels)
pressure_mask = convert_mask_to_pressure(mask,epsilon, concentrations, Nx, Ny);

sensor_data_holder = zeros(num_wavelength, Nx, Ny);

noisy_sensor_data_holder = zeros(num_noise,num_wavelength, Nx, Ny);


for n= 1:num_noise
    for w = 1:num_wavelength
    
    p_mask = squeeze(pressure_mask(w, c, :));

    sensor_data = calculate_sensor_data_for_mask(p_mask, Nx, dx, Ny, dy);
    if n== 1
        sensor_data_holder(w,:,:) = sensor_data;
    end

    noise = noises(n);
    noise_map = noise*randn(Nx,Ny);
   
    noisy_sensor_data = sensor_data + noise_map*noise_strength;
    noisy_sensor_data_holder(n,w,:,:) = noisy_sensor_data;

    end
end


all_noise_data = cell(num_noise);

for n = 1:num_noise
    recon_data_w = squeeze(noisy_sensor_data_holder(n,:,:,:));
    [sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type] = spectral_unmixing(epsilon, recon_data_w, type_names, plot_flag);
    spectral_unmixing_data = {sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type};
    all_noise_data(n) = spectral_unmixing_data;
    if save_flag

        wavelength_name = '';
        for w  =1:num_wavelength
            wavelength_name = wavelength_name + '-' + wavelength_names(w);
        end
        % add the datetime to str
        sat_file_name = sprintf('saturation_plots_%s.mat',wavelength_name);
        c_file_name = sprintf('concentration_plots_%s.mat',wavelength_name);
        w_avg_file_name = sprintf('weighted_avg_plots_%s.mat', wavelength_name);
        su_file_name = sprintf('spectral_unmixign_data_%s.mat',wavelength_name);


        save(su_file_name,'spectral_unmxing_data');
        %save(sat_file_name, 'saturations_by_type');
        %save(c_file_name, 'concentrations_by_type');
        %save(w_avg_file_name,'w_avg_by_type');
    end
end









% Plot results
figure;
hold on;
colors = lines(num_wavelength); % Generate distinct colors for each wavelength
wavelength_names = {'770-780', '750-850', '780-1030', '750-770-780-850'}; % Updated namesy
for w = 1:num_wavelength
    plot(noises, w_avg_holder(w,:, 1), '-o', 'Color', colors(w, :), ...
        'DisplayName', wavelength_names{w});
end

xlabel('Noise Level');
ylabel('HBO Weighted Average');
title('HBO Weighted Average vs. Noise Level for Different Wavelengths');
legend('show');
grid on;
hold off;