function [] = specunmix_analysis(specunmix_data, mask, wavelength_sets, concentrations, noise_simple, type_names, Nx ,Ny)
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

if nargin <1 
   data_770_780 = load('/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/spectral_unmixing_data_770780_20241219_.mat');
   data_750_850 = load('/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/spectral_unmixing_data_750850_20241219_.mat');
   
   specunmix_770_780 = data_770_780.specunmix_noise_c_data;
   specunmix_750_850 = data_750_850.specunmix_noise_c_data;

   specunmix_data = {specunmix_770_780,specunmix_750_850};
end

if nargin <2 
    mask = zeros(200,200);
    mask = mask + makeDisc(200,200, 100, 100, 4);
end

if nargin < 3
    wavelength_sets = [770, 780; 750,850];
   
end

if nargin < 4
    concentrations = [0.75, 0.25; 0.5, 0.5; 0.25 0.75;];
end

if nargin < 5
    noise_simple = 70*[0,0.05,0.5];
end
if nargin < 6
    type_names = {'Hb', 'HbO'};
end
if nargin < 7
    Nx = 200;
end

if nargin < 8
    Ny = 200;
end



num_wavelength_sets = length(specunmix_data);
num_type = size(concentrations,2);
num_concentration = size(concentrations,1);
num_noise = length(noise_simple);

mask_num = 0;
for x = 1:Nx
    for y = 1:Ny
        if mask(x,y) == 1
            mask_num = mask_num + 1;
        end
    end
end

saturation_avg_holder = zeros(num_wavelength_sets, num_type,num_noise, num_concentration, num_type);
concentration_avg_holder = zeros(num_wavelength_sets, num_type,num_noise, num_concentration, num_type);


for w = 1:num_wavelength_sets
    for n = 1:num_noise
        for c= 1:num_concentration
            nc_cell_arr = specunmix_data{1,w};
            saturation_data = nc_cell_arr{n,c,2};
            for t = 1:num_type
                saturation_type_data = squeeze(saturation_data(t,:,:));
                saturation_type_mask = saturation_type_data .* mask;
                sat_avg = sum(saturation_type_mask, 'all')/(mask_num);
                saturation_avg_holder(w,n,c,t) = sat_avg;
            end
        end
    end
end


for w = 1:num_wavelength_sets
    for n = 1:num_noise
        for c= 1:num_concentration
            nc_cell_arr = specunmix_data{w};
            concentration_data = nc_cell_arr{n,c,3};
            for t = 1:num_type
                concentration_type_data = concentration_data(t,:,:);
                concentration_type_mask = concentration_type_data .* mask;
                c_avg = sum(concentration_type_mask, 'all')/(mask_num);
                concentration_avg_holder(w,n,c,t) = c_avg;
            end
        end
    end
end

num_iter = 1;
num_c_unique = num_concentration/num_iter;

saturation_mean_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);
saturation_std_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);

concentration_mean_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);
concentration_std_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);

for w = 1:num_wavelength_sets
    for n = 1:num_noise
        for c= 1:num_c_unique
            for t = 1:num_type
               saturation_mean_holder(w,n,c,t) = mean(squeeze(saturation_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));
               saturation_std_holder(w,n,c,t) = std(squeeze(saturation_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));
               concentration_mean_holder(w,n,c,t) = mean(squeeze(concentration_avg_holder(w,n, (c-1)*num_iter +1 : c*num_iter, t)));
               concentration_std_holder(w,n,c,t) = std(squeeze(concentration_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));  
            end
        end
    end
end


%Plot the error_bars



% Plotting
offset = 0.1; % Offset value
colors = lines(num_wavelength_sets); % Generate distinct colors for each wavelength set



for t = 1:num_type
    for c = 1:num_c_unique
        figure;
        hold on;
        
        for w = 1:num_wavelength_sets
            sat_mean = squeeze(saturation_mean_holder(w,:,c,t));
            sat_std = squeeze(saturation_std_holder(w,:,c,t));
            plot(noise_simple, sat_mean,'-', 'LineWidth', 2);
            errorbar(noise_simple, sat_mean, sat_std, 'LineWidth', 2, 'Color', colors(w, :));
            
        end
        wv_title = num2str(wavelength_sets(w,:));
        title(sprintf('%s : %s Saturation Averages for c = %d ', wv_title, type_names{t}, c));
    end
end


for t = 1:num_type
    for c = 1:num_c_unique
        figure;
        hold on;
        for w = 1:num_wavelength_sets
            c_mean = squeeze(concentration_mean_holder(w,:,c,t));
            c_std = squeeze(concentration_std_holder(w,:,c,t));
            plot(noise_simple, c_mean,'o', 'LineWidth', 2);
            errorbar(noise_simple, c_mean, c_std, 'LineWidth', 2, 'Color', colors(w, :))
            wv_title = num2str(wavelength_sets(w,:));
            title(sprintf('%s : %s Concentration Averages for c = %d ', wv_title, type_names{t}, c));
        end
    end
end








% Calculate Error bars for saturation of each noise level over the mask.
% Calculate Error bars for weighted average of each noise level.
% Plot the Unmixed on the mask to check saturation



end