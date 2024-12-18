function [recon_holder_noise, noisy_sensor_data_holder] = build_pressures_w_noise(epsilon, concentrations, wavelengths, noise_levels, types, num_x_pixels, num_y_pixels, dx, dy)
    % Set default values for optional arguments
    if nargin < 1 || isempty(epsilon)
        epsilon = [1361 636; 1075 710];
    end

    if nargin < 2 || isempty(concentrations)
        concentrations = [0.5;0.5];
    end
    if nargin < 3 || isempty(wavelengths)
        wavelengths = [770, 780];
    end
    if nargin < 4 || isempty(noise_levels)
        noise_levels = [1];
    end
    if nargin < 5 || isempty(types)
        types = [1];
    end
    if nargin < 6 || isempty(num_x_pixels)
        num_x_pixels = 200;
    end
    if nargin < 7 || isempty(num_y_pixels)
        num_y_pixels = 200;
    end
    if nargin < 8 || isempty(dx)
        dx = 1e-6;
    end
    if nargin < 9 || isempty(dy)
        dy = 1e-6;
    end

    % Derived parameters
    num_wavelength = length(wavelengths);
    num_type = 2;
    num_concentration = size(concentrations, 2);
    num_noise = length(noise_levels);


    % Calculate the sensor_data array
    test_sd = calculate_sensor_data(1,epsilon(1,1), num_x_pixels, dx, num_y_pixels, dy);
   
    % Build Noise Matrix
    noise = randn(size(test_sd));

    % Pick Noise Matrix
    % need to implement it to pick the max value ???
    noise_mag = zeros(num_wavelength,num_type,num_concentration);

    % Holds the reconstructed images with noise added
    % Structure is (noise level, wavelength, type, concentration, x ,y)
    recon_holder_noise = zeros(num_noise, num_wavelength, num_type, num_concentration, num_x_pixels, num_y_pixels);
    noisy_sensor_data_holder = zeros(num_noise, num_wavelength, num_type, num_concentration, size(test_sd,1), size(test_sd,2));


    % Count how many iterations
    count = 0;
    num_iter = num_noise * num_type * num_wavelength * num_concentration;

    % Create a progress bar
    h = waitbar(0, 'Processing...');

    for n = 1:num_noise
        for w = 1:num_wavelength
            for h_type = 1:num_type
                for c = 1:num_concentration
                    if count == 0
                        timerVal = tic;
                    end

                    % Calculate sensor data
                    
                    sensor_data = calculate_sensor_data(concentrations(h_type, c), epsilon(w, h_type), num_x_pixels, dx, num_y_pixels, dy);
                    if n==1
                        max_val = max(sensor_data,[],'all');
                        noise_mag(w,h_type,c) = max_val;
                    end
                    
                    if max_val < 0
                      noise_mag(w,h_type,c) = 0;
                    end

                    % Add noise to sensor data
                    noise = randn(size(test_sd));
                    noisy_sensor_data = noise_mag(w,h_type,c) * noise_levels(n) * noise + sensor_data;
            
                    %{
                    figure;
                    subplot(2,2,1)
                    imagesc(sensor_data)
                    colormap;
                    colorbar;
                    subplot(2,2,2)
                    imagesc(noisy_sensor_data)
                    colormap;
                    colorbar;
                    sgtitle(sprintf('Noise Level %d',noise_levels(n)));
                    subplot(2,2,3);
                    plot(squeeze(sensor_data(100,:)));
                    subplot(2,2,4)
                    plot(squeeze(noisy_sensor_data(100,:)));
                    %}
                    % store the sensor data
                    
                    noisy_sensor_data_holder(n,w,h_type,c, :,:) = noisy_sensor_data;

                    % Reconstruct the image
                    recon_image_noise = reconstruct_image_k_wave(noisy_sensor_data, num_x_pixels, num_y_pixels, dx, dy);

                    % Store the reconstructed image
                    recon_holder_noise(n, w, h_type, c, :, :) = recon_image_noise;
            

                    figure;
                    imagesc(sensor_data);
                    xlabel('Time Delay (ns)');
                    ylabel('Transducers Postion (um)');
                 

                    figure;
                    imagesc(noisy_sensor_data);
                    xlabel('Time Delay (ns)');
                    ylabel('Transducers Postion (um)');
                 

                    figure;
                    plot((noisy_sensor_data(100,:)))
                    xlabel('Time (ns)');
                    ylabel('Pressue (a.u.)');

                    figure;
                    imagesc(recon_image_noise);
                    xlabel('X-Position (um)');
                    ylabel('Depth from Transducer (um)');

                    figure;
                    plot(recon_image_noise(100,:))
                    xlabel('Depth from Transducer (um)');
                    ylabel('Reconstructed Initial Pressure (a.u.)');
      

                    

                    % Estimate runtime for the first iteration
                    if count == 0
                        runtime = toc(timerVal);
                        total_time = num_iter * runtime;
                        fprintf('Time for one iteration: %.2f seconds\n', runtime);
                        fprintf('Estimated total time: %.2f seconds\n', total_time);
                    end

                    % Update progress
                    count = count + 1;
                    waitbar(count / num_iter, h, sprintf('Processing %d/%d...', count, num_iter));
                    fprintf('count: %d / %d\n', count, num_iter);
                end
            end
        end
    end

    % Close progress bar
    close(h);

    % Save the results
    wavelength_str = strjoin(string(wavelengths), '_');
    save(sprintf('k_recon_concentrations_w_noise_%s.mat', wavelength_str), 'recon_holder_noise');
    save(sprintf('k_sensor_data_w_noise_%s.mat', wavelength_str), 'noisy_sensor_data_holder');
end
