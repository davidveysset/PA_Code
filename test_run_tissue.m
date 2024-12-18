epsilon_770_780 = [1361 636; 1075 710];
epsilon_750_850 = [1405,518; 691,1050];
epsilon_780_1030 = [1075 710; 1024 206];
epsilon_all = [1405,518;1361, 636;1075, 710;691,1050];
wavelengths = {epsilon_770_780, epsilon_750_850, epsilon_780_1030, epsilon_all};
noises = [0,0.05,0.5];
num_noise = length(noises);
types = [1,2];
num_types = length(types);

Nx = 400;
Ny = 200;
dx = 5e-6;
dy = 5e-6;

[plot_map_og, plot_map_hb_og, plot_map_hbo_og] = create_tissue_mask_v2();

num_wavelength = length(wavelengths);
w_avg_holder = zeros(num_wavelength, num_noise, num_types);

% Initialize the wait bar
total_iterations = num_wavelength * num_noise;
h = waitbar(0, 'Processing, please wait...');

iteration = 0;

for w = 1:num_wavelength
    epsilon = wavelengths{w};
    for n = 1:num_noise
        noise = noises(n);
        % Run the tissue reconstruction and analysis
        [sum_C, hb_C, hbo_C, hb_S, hbo_S, w_avg_hbo_S, w_avg_hb_S] = run_tissue_recon_w_analysis(epsilon, noise, ...
            plot_map_og, plot_map_hb_og, plot_map_hbo_og, types, Nx, Ny, dx, dy);

        % Store results
        w_avg_holder(w, n, :) = [w_avg_hbo_S, w_avg_hb_S];

        % Update the wait bar
        iteration = iteration + 1;
        waitbar(iteration / total_iterations, h, ...
            sprintf('Processing: Wavelength %d/%d, Noise %d/%d', w, num_wavelength, n, num_noise));
    end
end

% Close the wait bar
close(h);

save('weighted_average_w_noise_wavelength.mat', 'w_avg_holder');

%{
w_avg_holder = load('weighted_average_w_noise_wavelength.mat');
w_avg_holder = w_avg_holder.w_avg_holder;
%}

% Plot results
figure;
hold on;
colors = lines(num_wavelength); % Generate distinct colors for each wavelength
wavelength_names = {'770-780', '750-850', '780-1030', '750-770-780-850'}; % Updated names
for w = 1:num_wavelength
    plot(noises, w_avg_holder(w, :, 1), '-o', 'Color', colors(w, :), ...
        'DisplayName', wavelength_names{w});
end

xlabel('Noise Level');
ylabel('HBO Weighted Average');
title('HBO Weighted Average vs. Noise Level for Different Wavelengths');
legend('show');
grid on;
hold off;
