


function [specunmix_noise_c_data, noisy_sensor_data_holder] = spectral_unmixing_simulation(mask, epsilon, wavelengths, concentrations, ...
  type_names, noise_levels, noise_strength, Nx, Ny, dx,dy, plot_flag_pressures, plot_flag_recon, plot_flag_concentrations, plot_flag_saturations, save_flag, folder_path) 
% Structure
% two wavelengths, two types hb, hbo
% construct the maps of the 9 points for hb and inverse for hbo
% Input also needs to consider the epsilon values, concentration_map,
% num_x_steps, num_y_steps


% The structure has to be for each epsilon (N wavelength), all
% concentrations, all types, all_noise generate a set cell array:
% {saturation, concentrations}  
% - saturation = size(num_noise, num_wavelength, num_concentration, Nx,Ny)
% - concentration = size(num_noise, num_wavelength, num_concentration, Nx,Ny)

% Then you run it for all the epsilons you are interested in ... multiple
% Es
% Then you run the cell array of cell the 2 cell arrays that does the
% analysis

if nargin < 1
    mask = zeros(200,200);
    mask = mask + makeDisc(200,200, 100, 100, 4);
end
if nargin <2
    epsilon = [1361 636; 1075 710]; %770-780
end
if nargin <3
    wavelengths = [770,780];
end
if nargin <4 
    concentrations = [0.25,0.75; 0.5, 0.5;];
end
if nargin <5
    type_names = {'Hb', 'HbO'};
end
if nargin < 6
    noise_levels = [0,0.05];
end
if nargin < 7
    noise_strength = 50;
end
if nargin < 8
    Nx = 200;
end
if nargin < 9
    Ny = 200;
end
if nargin < 10
    dx = 1e-6;
end
if nargin < 11
    dy = 1e-6;
end
if nargin < 12 
     plot_flag_pressures = false;
end
if nargin < 13
     plot_flag_concentrations = true;
end
if nargin < 14
     plot_flag_saturations = true;
end
if nargin < 15
    save_flag = true;
end
if nargin <16
    folder_path = '/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/';
end

num_wavelength = length(wavelengths);
num_noise = length(noise_levels);
num_concentrations = size(concentrations, 1);
num_type = size(concentrations,2);

noise_simple = noise_strength*noise_levels;

%Pressure mask is of size(num_concentrations, num_wavelength, # x pixels, # y pixels)
pressure_mask = convert_mask_to_pressure(mask, epsilon, concentrations, type_names, wavelengths, Nx, Ny, plot_flag_pressures);

%get size of sensor_data
p_mask = squeeze(pressure_mask(1, 1, :,:));
test_sensor_data = calculate_sensor_data_for_mask(p_mask, Nx, dx, Ny, dy);
num_transducers = size(test_sensor_data,1);
num_time = size(test_sensor_data,2);

sensor_data_holder = zeros(num_wavelength, num_concentrations, num_transducers , num_time );
noisy_sensor_data_holder = zeros(num_noise,num_wavelength, num_concentrations, num_transducers, num_time);
recon_image_holder = zeros(num_noise,num_wavelength, num_concentrations, Nx, Ny);


for w = 1:num_wavelength
    for c=1:num_concentrations
        p_mask = squeeze(pressure_mask(w, c, :,:));
        sensor_data = calculate_sensor_data_for_mask(p_mask, Nx, dx, Ny, dy);
        sensor_data_holder(w,c,:,:) = sensor_data;

    end
end



for n= 1:num_noise
    for w = 1:num_wavelength
        for c=1:num_concentrations
            sensor_data = squeeze(sensor_data_holder(w,c,:,:));
            noise = noise_levels(n);
            noise_map = noise_strength*noise*randn(num_transducers,num_time);
            noisy_sensor_data = sensor_data + noise_map;
            noisy_sensor_data_holder(n,w,c,:,:) = noisy_sensor_data;
            recon_image_holder(n,w,c,:,:) = reconstruct_image_k_wave(noisy_sensor_data, Nx, Ny, dx, dy);

        end
    end

end


if plot_flag_recon
   
   for n= 1:num_noise
        figure;
        count = 1;
        max_pressure = max(recon_image_holder(n,:,:,:,:),[],'all');
        for w = 1:num_wavelength
            for c=1:num_concentrations
            subplot(num_wavelength, num_concentrations,count)
            imagesc(squeeze(recon_image_holder(n,w,c,:,:)));
            colormap;
            colorbar;
            clim([0,max_pressure]);
            title(sprintf('Reconstructed at N= %d, wv = %d, c= %d', noise_simple(n), wavelengths(w), c));
            count = count + 1;
           
            end
        end
    end
end

specunmix_noise_c_data = cell(num_noise, num_concentrations,4);

for n = 1:num_noise
    for c=1:num_concentrations
        recon_data_w = squeeze(recon_image_holder(n,:,c,:,:));
        [sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type] = spectral_unmixing(epsilon, recon_data_w,type_names, Nx, Ny);
        su_concentration_data = {sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type};
        specunmix_noise_c_data(n,c,:) = su_concentration_data;
    end
end

if plot_flag_concentrations


for n = 1:num_noise
    for c=1:num_concentrations
        figure;
        su_concentration_data = specunmix_noise_c_data(n,c,:);
        concentrations_by_type = su_concentration_data{3};
        max_val = max(concentrations_by_type,[],'all');
        for t =1:num_type
            subplot(1,num_type, t);
            c_map = squeeze(concentrations_by_type(t,:,:));
            imagesc(c_map);
            title(sprintf('%s Concentration, c = %d, n = %d ', type_names{t}, c, noise_simple(n)));
            colormap;
            colorbar;
            clim([0,max_val]);
            xlabel('X-Axis');
            ylabel('Y-Axis');

        end
    end
end

end 

if plot_flag_saturations
for n = 1:num_noise
    for c=1:num_concentrations
        figure;
        su_concentration_data = specunmix_noise_c_data(n,c,:);
        saturations_by_type = su_concentration_data{2};
        max_val = max(saturations_by_type,[],'all');
        for t =1:num_type
            subplot(1,num_type, t);
            c_map = squeeze(saturations_by_type(t,:,:));
            imagesc(c_map);
            title(sprintf('%s Saturation, c = %d, n = %d ', type_names{t}, c, noise_simple(n)));
            colormap;
            colorbar;
            clim([0,max_val]);
            xlabel('X-Axis');
            ylabel('Y-Axis');

        end
    end
end

end


if save_flag
    rand_id = num2str(rand(0,1000));
    timestamp = datestr(datetime('now'), 'yyyymmdd');
    timestamp = [timestamp,'_', rand_id];
    wavelength_name = '';
    for w =1:num_wavelength
        wavelength_name = [wavelength_name,num2str(wavelengths(w))] ;
    end
    % add the datetime to str
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
    su_file_name = [folder_path,'spectral_unmixing_data_',wavelength_name, '_', timestamp, '.mat'];
    save(su_file_name,'specunmix_noise_c_data');

end






