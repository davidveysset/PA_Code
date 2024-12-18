function plot_error_bars(data)

% Structure of all_circle_data = (num_circle, num_wavelengths, num_noise)
% Dimensions of the data
[num_circles, num_wavelength_sets, num_noise_levels] = size(data);

% Preallocate arrays for mean and standard deviation
mean_values = zeros(num_wavelength_sets, num_noise_levels);
std_values = zeros(num_wavelength_sets, num_noise_levels);

% Calculate mean and standard deviation for each wavelength and noise level
for w = 1:num_wavelength_sets
    for n = 1:num_noise_levels
        % Extract values for all circles at this wavelength and noise level
        circle_values = data(:, w, n);
        mean_values(w, n) = mean(circle_values);
        std_values(w, n) = std(circle_values);
    end
end
% Plotting
offset = 0.01; % Offset value
figure;
hold on;
for w = 1:num_wavelength_sets
    x_values = (1:num_noise_levels) + (w - (num_wavelength_sets+1)/2) * offset; % Offset x-values
    errorbar(x_values, mean_values(w, :), std_values(w, :), ...
        'LineWidth', 2, 'DisplayName', sprintf('Wavelength Set %d', w));
end
xlabel('Noise Level Index');
ylabel('Value');
title('Error Bars for Each Wavelength Set');
legend('show');
grid on;
hold off;
end