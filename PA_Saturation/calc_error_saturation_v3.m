function [C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration, saturation_error,saturation_error_hb, w_avg] = calc_error_saturation_v3(data, E, wavelengths, expected_values, plot_hbo, plot_hb, plot_analysis)

    % Create Title for Plots Based on Wavelength
    w_title = strjoin(string(wavelengths), ' ');

    % Rename because I'm dumb and this is what I called it
    pg_noise_set = data;
    % Sum along the types to combine the oxy and deoxy data

    % Get parameters for the function
    L_x = size(pg_noise_set, 4);
    L_y = size(pg_noise_set, 5);
    num_concentration = size(pg_noise_set, 3);
    num_wavelength = size(pg_noise_set, 2);
    num_noise = size(pg_noise_set, 1);
    num_types = size(E, 2);

    % Create Pressure Matrix as a function of (wavelength, x ,y)
    P = zeros(num_wavelength, L_x, L_y);

    % C_nnls is concentrations from non-negative least squares as a function of
    % (type, x, y)
    C_nnls = zeros(num_noise, num_concentration, num_types, L_x, L_y);

    % Change the convergence since it was causing error
    opts = optimset('lsqnonneg');
    opts.TolX = 1e-6;

    % Total iterations for progress bar
    total_iterations = num_noise * num_concentration * L_x * L_y;
    current_iteration = 0;

    % Create a progress bar
    h = waitbar(0, 'Calculating...');

    % Loop through each noise level, concentration, and optimize Arg_Min{P - EC}
    % and Store in C_nnls Matrix
    for n = 1:num_noise
        for c = 1:num_concentration
            P_nc = squeeze(pg_noise_set(n, :, c, :, :));
            for x = 1:L_x
                for y = 1:L_y
                    P_point = squeeze(P_nc(:, x, y));
                    C_point = lsqnonneg(E, P_point, opts);
                    C_nnls(n, c, :, x, y) = C_point;

                    % Update progress bar
                    current_iteration = current_iteration + 1;
                    if mod(current_iteration,1000) == 0
                    waitbar(current_iteration / total_iterations, h, ...
                        sprintf('Processing Noise %d, Concentration %d, Pixel (%d, %d)', n, c, x, y));
                    end
                end
            end
        end
    end

    % Close the progress bar
    close(h);

    C_nnls(C_nnls < 0) = 1e-9;

    % Calculate the total concentration
    sum_C = squeeze(sum(C_nnls, 3));

    % Calculate the hbO and hb Concentration
    hbo_C = squeeze(C_nnls(:, :, 1, :, :));
    hb_C = squeeze(C_nnls(:, :, 2, :, :));

    % Get rid of zeros so that in the saturation calculation no dividing by zero
    sum_C(sum_C == 0) = 1e-9;

    % Calculate the saturation
    hbo_S = hbo_C ./ sum_C;
    hb_S = hb_C ./ sum_C;

% Loop through each noise level, and concetration and plot the saturation
% for HbO
if plot_hbo
    for n = 1:num_noise
        figure; % Create a new figure for each noise level
        max_val =max(hbo_S(n, :, :, :),[],'all');
        
        for c = 1:num_concentration

            subplot(ceil(sqrt(num_concentration)), ceil(sqrt(num_concentration)), c); % Subplots for concentrations
            imagesc(squeeze(hbo_S(n, c, :, :))); % Show HbO concentrations as an image
            colormap turbo;
            colorbar; % Add colorbar for reference
            title(sprintf('Concentration %d for Noise Level %d  %s', c, n, w_title));
            xlabel('X');
            ylabel('Y');
            clim([0,max_val]);
        end
        sgtitle(sprintf('HbO Saturation for Noise Level %d  %s', n, w_title)); % Super-title for the figure
    end
end

if plot_hb
    for n = 1:num_noise
        figure; % Create a new figure for each noise level
        [max_val,max_idx] =max(hb_S(n, :, :, :),[],'all');
        for c = 1:num_concentration

            subplot(ceil(sqrt(num_concentration)), ceil(sqrt(num_concentration)), c); % Subplots for concentrations
            imagesc(squeeze(hb_S(n, c, :, :))); % Show HbO concentrations as an image
            colormap turbo;
            colorbar; % Add colorbar for reference
            title(sprintf('Concentration %d for Noise Level %d  %s', c, n, w_title));
            xlabel('X');
            ylabel('Y');
            clim([0,max_val]);
        end
        sgtitle(sprintf('Hb Saturation for Noise Level %d  %s', n, w_title)); % Super-title for the figure
    end
end

% create an array to store the circle concetnrations
circle_concentration = zeros(num_noise, num_concentration);
%create the saturation error matrix  to store vals
saturation_error = zeros(size(circle_concentration));
saturation_error_hb = zeros(size(circle_concentration));
% create the weighted_average_matrix 
w_avg = zeros(num_noise,num_concentration);

% circle parameters
radius = 4;
center_x = floor(L_x / 2);
center_y = floor(L_y / 2);

% generate indices for the circular mask
[X, Y] = meshgrid(1:L_x, 1:L_y);
distance_from_center = sqrt((X - center_x).^2 + (Y - center_y).^2);
circle_mask = distance_from_center <= radius;
%Remove bottom point
circle_mask(104,100) = 0;


% loop through each noise levels and concentration
for n = 1:num_noise
    for c = 1:num_concentration
        % Extract the 2D grid for the current noise level and concentration
        position_grid = squeeze(hbo_S(n, c, :, :));
        position_grid_hb = squeeze(hb_S(n,c,:,:));
        % Apply the circular mask and average the values
        circle_values = position_grid(circle_mask);
        circle_values_hb = position_grid_hb(circle_mask);
        % Store the average value
        avg_val = mean(circle_values, 'omitnan');
        avg_val_hb = mean(circle_values_hb,'omitnan');
        circle_concentration(n, c) = avg_val;
        %Store the value of the error in saturation error
        saturation_error(n,c) =  avg_val;
        saturation_error_hb(n,c) = avg_val_hb;
        w_avg_calc = sum((squeeze(sum_C(n,c,:,:)) .* position_grid),'all')/(sum(sum_C(n,c,:,:),'all'));
        w_avg(n,c) = w_avg_calc;
        %saturation_error(n,c) = w_avg_calc;

    end
end

if plot_analysis

    % Create a figure for bar plots

    figure;

    for n = 1:num_noise
        subplot(ceil(sqrt(num_noise)), ceil(sqrt(num_noise)), n); % Arrange subplots in a grid

        % Combine data for current noise level with expected values
        num_expected = length(expected_values);
        combined_data = [circle_concentration(n, :); expected_values]';

        % Bar plot
        bar(combined_data, 'grouped'); % Grouped bar plot

        % Formatting the subplot
        title(sprintf('Noise Level %d   %s', n, w_title));
        xlabel('Concentration Index');
        ylabel('HbO Average & Expected');
        xticks(1:num_concentration);
        xticklabels(arrayfun(@num2str, 1:num_concentration, 'UniformOutput', false));
        ylim([0 max([circle_concentration(:); expected_values(:)]) * 1.1]); % Uniform Y-axis limits

        % Add legend
        legend({'Measured', 'Expected'}, 'Location', 'northwest');
    end

    % Add a super-title for the entire figure
    sgtitle(sprintf('Average HbO Concentrations vs Expected Values for Each Noise Level %s', w_title));


    figure;

    for n = 1:num_noise
        subplot(ceil(sqrt(num_noise)), ceil(sqrt(num_noise)), n); % Arrange subplots in a grid

        % Bar plot for the current noise level
        bar(saturation_error(n, :));
        % Formatting the subplot
        title(sprintf('Noise Level %d %s', n, w_title));
        xlabel('Concentration Index');
        ylabel('Average HbO');
        xticks(1:num_concentration); % Set x-axis ticks to match concentration indices
        ylim([min(saturation_error(:)) max(saturation_error(:)) ]); % Set uniform y-axis limits for comparison
    end

    % Add a super-title for the entire figure
    sgtitle(sprintf('Error In Saturation %s', w_title));
end
end

