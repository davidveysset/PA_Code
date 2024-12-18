function plot_error_bars(data, type)

if nargin < 1
   data = load('all_circle_data_hbo.mat');
   data = data.all_circle_data_hbo;
end

if nargin < 2
    type = 'HBO'
end
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
        mean_values(w, n) = 100 * mean(circle_values);
        std_values(w, n) = 100 * std(circle_values);
    end
end

% Plotting
offset = 0.1; % Offset value
wavelength_names = {'750-850', '770-780', '780-1030', '750-770-780-850'}; % Updated names
colors = lines(num_wavelength_sets); % Generate distinct colors for each wavelength set

figure;
hold on;

for w = 1:num_wavelength_sets
    x_values = (1:num_noise_levels) + (w - (num_wavelength_sets + 1) / 2) * offset; % Offset x-values
    % Plot error bars
    errorbar(x_values, mean_values(w, :), std_values(w, :), ...
        'LineWidth', 2, 'Color', colors(w, :), 'DisplayName', wavelength_names{w});
    % Plot points without adding to legend
    plot(x_values, mean_values(w, :), 'o', 'LineWidth', 2, 'Color', colors(w, :), 'HandleVisibility', 'off');
end

% Adjust x-ticks and labels
xticks(1:num_noise_levels);
xticklabels({'None', '24.3', '12.2', '6.1', '2.4'});
xlabel('Noise Strength', 'FontSize', 20);
ylabel(sprintf('Saturation Value(%%) %s ' , type), 'FontSize', 20);


ax = gca;
ax.FontSize = 20;

% Set legend to top-left corner
legend('show', 'Location', 'northwest', 'FontSize', 20);

grid on;
hold off;

end

