% Structure
% two wavelengths, two types hb, hbo
% construct the maps of the 9 points for hb and inverse for hbo
% Input also needs to consider the epsilon values, concentration_map,
% num_x_steps, num_y_steps

epsilon_770_780 = [1361 636; 1075 710];
epsilon_750_850 = [1405,518; 691,1050];
epsilon_780_1030 = [1075 710; 1024 206];
epsilon_all = [1405,518;1361, 636;1075, 710;691,1050];
wavelengths = {epsilon_770_780};
num_wavelength = length(wavelengths);
noises = 0;
num_noise = length(noises);
types = [1,2];
num_types = length(types);

Nx = 400;
Ny = 200;
dx = 5e-6;
dy = 5e-6;



hold on;


w_avg_holder = zeros(num_wavelength,num_noise,num_types);

[plot_map_og, plot_map_hb_og, plot_map_hbo_og] = create_tissue_mask_v2();

for w = 1:num_wavelength
    epsilon = wavelengths{w};
    for n = 1:num_noise
    noise = noises(n);
    noise_map = noise*randn(size(plot_map_og));
    [sum_C, hb_C, hbo_C, hb_S, hbo_S, w_avg_hbo_S, w_avg_hb_S] = run_tissue_recon_w_analysis(epsilon, noise, ...
        plot_map_og, plot_map_hb_og, plot_map_hbo_og, types, Nx, Ny, dx,dy);
    
    w_avg_holder(w,n,:) = [w_avg_hbo_S, w_avg_hb_S];
    
    end
end

figure;
subplot(1,2,1)a
imagesc(hb_S)
title('Hb Saturation')
colormap;
colorbar
subplot(1,2,2)
imagesc(hbo_S)
title('HbO Saturation')
colormap;
colorbar



save('weighted_average_w_noise_wavelength.mat','w_avg_holder');

%{
w_avg_holder = load('weighted_average_w_noise_wavelength.mat');
w_avg_holder = w_avg_holder.w_avg_holder;
%}


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