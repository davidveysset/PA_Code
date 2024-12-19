function [full_specunmix_analysis_data] = specunmix_analysis(specunmix_data, mask, wavelength_sets, concentrations, noise_simple, type_names, Nx ,Ny, ...
    save_flag, num_iter, folder_path, plot_flag_c, plot_flag_s, plot_flag_w_avg, plot_flag_diff)
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

if nargin < 9
    save_flag = false;
end
if nargin < 10
    num_iter = 1;
end

if nargin < 11
    folder_path = '/Users/calvinsmith/Bouma_lab/PA_project/PA_Saturation/PA_saturation_data/';
end

if nargin < 12
    plot_flag_s = true;
end
if nargin < 13
    plot_flag_c = true;
end

if nargin < 14
    plot_flag_w_avg = true;
end

if nargin < 15
    plot_flag_diff = true;
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
w_avg_holder = zeros(num_wavelength_sets, num_type,num_noise, num_concentration, num_type);

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


for w = 1:num_wavelength_sets
    for n = 1:num_noise
        for c= 1:num_concentration
            nc_cell_arr = specunmix_data{1,w};
            w_avg_by_type = nc_cell_arr{n,c,4};
            for t = 1:num_type
                w_avg_holder(w,n,c,t) = w_avg_by_type(t);
            end
        end
    end
end


num_c_unique = num_concentration/num_iter;

saturation_mean_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);
saturation_std_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);

concentration_mean_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);
concentration_std_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);

w_avg_mean_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);
w_avg_std_holder = zeros(num_wavelength_sets, num_noise, num_c_unique, num_type);

sat_difference_holder = size(saturation_mean_holder);

for w = 1:num_wavelength_sets
    for n = 1:num_noise
        for c= 1:num_c_unique
            for t = 1:num_type
               saturation_mean_holder(w,n,c,t) = mean(squeeze(saturation_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));
               saturation_std_holder(w,n,c,t) = std(squeeze(saturation_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));
               concentration_mean_holder(w,n,c,t) = mean(squeeze(concentration_avg_holder(w,n, (c-1)*num_iter +1 : c*num_iter, t)));
               concentration_std_holder(w,n,c,t) = std(squeeze(concentration_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));  
                
               w_avg_mean_holder(w,n,c,t) =  mean(squeeze(w_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));
               w_avg_std_holder(w,n,c,t) = std(squeeze(w_avg_holder(w,n, (c-1)*num_iter + 1: c*num_iter, t)));

               sat_difference_holder(w,n,c,t) = saturation_mean_holder(w,n,c,t) - concentrations(c*num_iter,t);
            end
        end
    end
end

full_specunmix_analysis_data = {saturation_mean_holder, saturation_std_holder, concentration_mean_holder, concentration_std_holder, w_avg_mean_holder, w_avg_std_holder, sat_difference_holder};

if save_flag 
    rand_id = num2str(round(1000*rand(1)));
    timestamp = datestr(datetime('now'), 'yyyymmdd');
    timestamp = [timestamp,'_', rand_id];
 
    % add the datetime to str
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
    su_file_name = [folder_path,'full_specunmix_analysis_data', '_', timestamp, '.mat'];
    save(su_file_name,'full_specunmix_analysis_data');


end


%Plot the error_bars

colors = lines(num_wavelength_sets); % Generate distinct colors for each wavelength set

if plot_flag_s
for t = 1:num_type
    for c = 1:num_c_unique
        figure;
        hold on;

        for w = 1:num_wavelength_sets
            sat_mean = squeeze(saturation_mean_holder(w,:,c,t));
            sat_std = squeeze(saturation_std_holder(w,:,c,t));
            h = plot( sat_mean, '-', 'LineWidth', 2, 'Color', colors(w, :));
            errorbar(sat_mean, sat_std, 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', colors(w, :), 'DisplayName', num2str(wavelength_sets(w, :)));
        end

        title(sprintf('%s Saturation Averages for c = %d', type_names{t}, c), 'Interpreter', 'none');
        legend('show', 'Location', 'best', 'FontSize', 14);
        xlabel('Noise Level');
        ylabel('Saturation Mean');
        
        set(gca, 'FontSize',14)
        xticks(1:num_noise);
        xticklabels(arrayfun(@num2str, noise_simple, 'UniformOutput', false));
        grid on;
    end
end
end

if plot_flag_c
for t = 1:num_type
    for c = 1:num_c_unique
        figure;
        hold on;

        for w = 1:num_wavelength_sets
            c_mean = squeeze(concentration_mean_holder(w,:,c,t));
            c_std = squeeze(concentration_std_holder(w,:,c,t));

            % Plot the data with markers and error bars
            plot( c_mean, 'o-', 'LineWidth', 2, 'Color', colors(w, :));
            errorbar( c_mean, c_std, 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', colors(w, :),'DisplayName', num2str(wavelength_sets(w, :)));
        end

        title(sprintf('%s Concentration Averages for c = %d', type_names{t}, c), 'Interpreter', 'none');
        legend('Location', 'best', 'FontSize', 14);
        xlabel('Noise Level');
        ylabel('Concentration Mean');
        set(gca, 'FontSize',14)
        % Set x-ticks to match noise_simple values
        xticks(1:num_noise);
        xticklabels(arrayfun(@num2str, noise_simple, 'UniformOutput', false));

        % Enable grid
        grid on;
    end
end

end

if plot_flag_w_avg
for t = 1:num_type
    for c = 1:num_c_unique
        figure;
        hold on;

        for w = 1:num_wavelength_sets
            % Calculate mean and std for w_avg
            w_mean = squeeze(w_avg_mean_holder(w,:,c,t));
            w_std = squeeze(w_avg_std_holder(w,:,c,t));
            
            % Plot the data with markers and error bars
            plot( w_mean, 'o-', 'LineWidth', 2, 'Color', colors(w, :), 'DisplayName', num2str(wavelength_sets(w, :)));
            errorbar( w_mean, w_std, 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', colors(w, :));
        end

        % Add title, legend, and axis labels
        title(sprintf('%s W_Avg Averages for c = %d', type_names{t}, c), 'Interpreter', 'none');
        legend('Location', 'best', 'FontSize', 14);
        xlabel('Noise Level');
        ylabel('w_{avg} Mean');
        set(gca, 'FontSize',14)
        % Set x-ticks to match noise_simple values
        xticks(1:num_noise);
        xticklabels(arrayfun(@num2str, noise_simple, 'UniformOutput', false));

        % Enable grid
        grid on;
    end
end
end

if plot_flag_diff

for t = 1:num_type
    for c = 1:num_c_unique
        % Create a new figure for each type and concentration
        figure;
        hold on;

        % Prepare data for the bar graph
        bar_data = zeros(num_noise, num_wavelength_sets); % [noise_levels x wavelength_sets]

        % Populate bar_data
        for n = 1:num_noise
            for w = 1:num_wavelength_sets
                bar_data(n, w) = sat_difference_holder(w, n, c, t);
            end
        end

        % Create grouped bar graph
        b = bar(categorical(1:num_noise), bar_data, 'grouped'); % X-axis is categorical for noise levels

        % Set bar colors based on wavelength sets
        for w = 1:num_wavelength_sets
            b(w).FaceColor = 'flat';
            b(w).CData = colors(w, :); % Assign colors for each wavelength set
        end

        % Create legend entries for the wavelength sets
        legend_entries = cell(1, num_wavelength_sets);
        for w = 1:num_wavelength_sets
            legend_entries{w} = sprintf('Wavelength Set: [%s]', num2str(wavelength_sets(w, :)));
        end

        % Add title, axis labels, and legend
        title(sprintf('Difference from Expected Type: %s, C: %d', type_names{t}, c), 'FontSize', 16, 'Interpreter', 'none');
        xlabel('Noise Levels');
        ylabel('Saturation Difference');
        legend(legend_entries, 'Location', 'best', 'FontSize', 12);
        set(gca, 'FontSize',14)
        % Set x-axis tick labels to match noise levels
        xticklabels(arrayfun(@num2str, noise_simple, 'UniformOutput', false));

        % Enable grid
        grid on;
    end
end
end

end