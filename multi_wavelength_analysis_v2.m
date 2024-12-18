function [all_circle_data_hbo, all_circle_data_hb] = multi_wavelength_analysis_v2(data, wavelengths)
% Multi-wavelength analysis visualization function
%
% Parameters:
%   - data: Cell array containing wavelength sets. Each set contains cells
%           with the matrix returned from calc_error_saturation.
%   - wavelengths: Array of wavelengths to label in the plots.
%
% Data structure:
%   - data{w} contains cells, where index 6 holds the saturation error matrix.
%   - Saturation error matrix (noise levels x circles).

for n=1:2
if n ==1 
    sat_err_idx = 8;  % Index for saturation error in the data structure
end
if n==2
    sat_err_idx = 9;
end
reference_data = data{1};
reference_set = reference_data{sat_err_idx};

% Prepare data dimensions
num_concentration = size(reference_set, 2);  % Number of concentrations (circles)
num_noise_levels = size(reference_set, 1);   % Number of noise levels
num_wavelength_sets = length(wavelengths);   % Number of wavelength sets

% Create a figure for bar plots
figure;
sgtitle('Saturation Error Across Noise Levels for Each Concentration');

for c = 1:num_concentration
    % Data for the current concentration circle
    data_for_circle = zeros(num_wavelength_sets, num_noise_levels);

    % Extract data for each wavelength set
    for w = 1:num_wavelength_sets
        data_list = data{w};
        saturation_error = data_list{sat_err_idx};
        data_for_circle(w, :) = saturation_error(:, c);
    end

    % Create a subplot for the current concentration
    subplot(ceil(sqrt(num_concentration)), ceil(sqrt(num_concentration)), c);
    bar(data_for_circle', 'grouped');

    % Format subplot
    title(sprintf('Concentration %d', c));
    xlabel('Noise Level');
    ylabel('Saturation Error');
    xticks(1:num_noise_levels);
    xticklabels(arrayfun(@(n) sprintf('N%d', n), 1:num_noise_levels, 'UniformOutput', false));
%    ylim([-min(data_for_circle(:)), max(data_for_circle(:)) * 1.1]);
    legend(wavelengths, 'Location', 'northwest');
end

% Create a figure for line plots
figure;
sgtitle('Saturation Error Across Noise Levels for Each Circle by Wavelength');

if n== 1
all_circle_data = zeros(num_concentration, num_wavelength_sets, num_noise_levels);
end
if n==2 
    all_circle_data_hbo = all_circle_data;
    all_circle_data = zeros(num_concentration, num_wavelength_sets, num_noise_levels);
end

for c = 1:num_concentration
    % Data for the current concentration circle
    data_for_circle = zeros(num_wavelength_sets, num_noise_levels);

    % Extract data for each wavelength set
    for w = 1:num_wavelength_sets
        data_list = data{w};
        saturation_error = data_list{sat_err_idx};
        data_for_circle(w, :) = saturation_error(:, c);
    end
    all_circle_data(c, :,:) = data_for_circle;

    % Create a subplot for the current concentration
    subplot(ceil(sqrt(num_concentration)), ceil(sqrt(num_concentration)), c);
    hold on;

    % Plot line graphs for each wavelength set
    for w = 1:num_wavelength_sets
        plot(1:num_noise_levels, squeeze(data_for_circle(w, :)), '-o', 'LineWidth', 2, ...
            'DisplayName', wavelengths{w});
    end
    hold off;

    % Format subplot
    title(sprintf('Concentration %d', c));
    xlabel('Noise Level');
    ylabel('Saturation Error');
    xticks(1:num_noise_levels);
    xticklabels(arrayfun(@(n) sprintf('N%d', n), 1:num_noise_levels, 'UniformOutput', false));
    legend('show');
    grid on;
    
    if c == 10
        % Noise levels and their inverses
        noise_levels = [0, 0.05, 0.1, 0.2, 0.5, 1];
        inverse_noise_levels = [Inf, 1 ./ noise_levels(2:end)]; % Handle division by zero for the first element

        figure;
        hold on;

        % Prepare data for plotting (scaled by 100)
        data_for_circle_c3 = 100 * data_for_circle(:, 1:end-1);

        % Plot line graphs for each wavelength set
        for w = 1:num_wavelength_sets
            plot(1:num_noise_levels-1, data_for_circle_c3(w, :), '-o', 'LineWidth', 2, ...
                'DisplayName', wavelengths{w});
        end

        % Format plot
        %title(sprintf('Concentration %d Oxygenated Hemoglobin', 50));
        xlabel('Signal to Noise Ratio','FontSize', 24);
        ylabel('Saturation Error(%)', 'FontSize',24);

        % Update xticks and labels
        xticks(1:num_noise_levels - 1); % Exclude the last tick since data_for_circle ends at num_noise_levels-1
        xticklabels(["None", arrayfun(@(n) sprintf('%.1f', n), inverse_noise_levels(2:end), 'UniformOutput', false)]);
        ax = gca; % Get current axes handle
        ax.FontSize = 20; % Increase tick label font size

        legend('show','FontSize',14);
        grid on;

        figure;
        
       
    end
    if n==2
        all_circle_data_hb = all_circle_data;
    end

end
end

end