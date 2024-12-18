function [C_nnls, hbo_C, sum_C, hbo_S, circle_concentration,saturation_error] = calc_error_saturation(data, E, wavelengths)

w_title = '';
w_title = strjoin(string(wavelengths), '_');
pg_noise_set = data;
pg_noise_set = squeeze(sum(pg_noise_set,3));

L_x = size(pg_noise_set,4);
L_y = size(pg_noise_set,5);
num_concentration = size(pg_noise_set,3);
num_wavelength = size(pg_noise_set,2);
num_noise = size(pg_noise_set,1);

P = zeros(num_wavelength,L_x, L_y);
num_types = size(E,2);

C_nnls = zeros(num_noise, num_concentration,num_types,L_x, L_y);
opts = optimset('lsqnonneg');
opts.TolX = 1e-6;
for n = 1:num_noise
    for c = 1:num_concentration
        P_nc = squeeze(pg_noise_set(n,:,c,:,:));
        for x = 1:L_x
            for y = 1:L_y
                P_point = squeeze(P_nc(:,x, y));
                C_point = lsqnonneg(E, P_point, opts);
                C_nnls(n,c,:, x, y) = C_point;
            end
        end
    end

end


%data = load('unmixed_k_recon_concentrations.mat','C_nnls');
%C_nnls = data.C_nnls;

sum_C = squeeze(sum(C_nnls, 3));
hbo_C = squeeze(C_nnls(:,:,1,:,:));
sum_C(sum_C == 0) = 1e-9;
hbo_S = hbo_C./sum_C;




for n = 1:num_noise
    figure; % Create a new figure for each noise level
    [max_val,max_idx] =max(hbo_S(n, :, :, :),[],'all');
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


circle_concentration = zeros(num_noise, num_concentration);

% Circle parameters
radius = 4;
center_x = floor(L_x / 2);
center_y = floor(L_y / 2);

% Generate indices for the circular mask
[X, Y] = meshgrid(1:L_x, 1:L_y);
distance_from_center = sqrt((X - center_x).^2 + (Y - center_y).^2);
circle_mask = distance_from_center <= radius;

saturation_error = zeros(size(circle_concentration));
expected_values = [0, 0.25, 0.5, 0.75, 1];

% Loop through noise levels and concentrations
for n = 1:num_noise
    for c = 1:num_concentration
        % Extract the 2D grid for the current noise level and concentration
        current_grid = squeeze(hbo_S(n, c, :, :));
        
        % Apply the circular mask and average the values
        circle_values = current_grid(circle_mask);
        avg_val = mean(circle_values, 'omitnan'); % Store the average value
        circle_concentration(n, c) = avg_val;
        saturation_error(n,c) = expected_values(c) - avg_val;
    end
end


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




% load the data in ...
% data is organized recon_holder_noise = (num_noise, num_wavelengths, num_concentrations, x, y)


wavelengths = [750, 850];
%data = load('k_recon_concentrations_w_noise_750_850.mat','recon_holder_noise');
%test = data.recon_holder_noise;
%data_750_850 = test;
%E_750_850 = [1405,518;691,1050]; % For 750 and 850


data = load('k_recon_concentrations_w_noise.mat','recon_holder_noise');
test = data.recon_holder_noise;
data_770_780 = test;
E_770_780 = [1361 636; 1075 710]; % For 770 and 780



data = load('k_recon_concentrations_w_noise_all.mat','recon_holder_noise');
test = data.recon_holder_noise;
data_all= test;
E_all = [1405,518; 1361 ,636; 1075,710; 691,1050];

%err_cbn_oxy_750_850 = calc_error_saturation_v2(data_750_850, E_750_850);


[C_nnls, hbo_C, hbo_S, sum_C, circle_concentration,saturation_error] = calc_error_saturation_k(data_770_780, E_770_780);
unmixing_results_770_780 = {C_nnls, hbo_C, hbo_S, sum_C, circle_concentration,saturation_error};
save('unmixing_results_770_780_v1.mat','unmixing_results_770_780');

%{
[C_nnls, hbo_C, hbo_S, sum_C, circle_concentration,saturation_error] = calc_error_saturation_k(data_750_850, E_750_850);
unmixing_results_750_850 = {C_nnls, hbo_C, hbo_S, sum_C, circle_concentration,saturation_error};
save('unmixing_results_750_850_v1.mat','unmixing_results_750_850');
%err_cbn_oxy_all = calc_error_saturation(data_all, E_all)
%}
%{
[C_nnls, hbo_C, hbo_S, sum_C, circle_concentration,saturation_error] = calc_error_saturation_k(data_all, E_all);
unmixing_results_all = {C_nnls, hbo_C, hbo_S, sum_C, circle_concentration,saturation_error};
save('unmixing_results_all_v1.mat','unmixing_results_all');
%}