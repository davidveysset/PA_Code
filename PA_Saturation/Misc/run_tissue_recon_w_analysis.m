function [sum_C, hb_C, hbo_C, hb_S, hbo_S, w_avg_hbo_S, w_avg_hb_S] = run_tissue_recon_w_analysis(epsilon, noise, plot_map_og, plot_map_hb_og, ...
    plot_map_hbo_og, types,Nx, Ny, dx,dy)

num_wavelength = length(epsilon);
num_type = length(types);
recon_holder_image = zeros(num_wavelength,Nx, Ny);

count = 1;
noise_map_og = randn(size(plot_map_og));
for w = 1:num_wavelength
    noise_map_og = randn(size(plot_map_og));

    plot_map_hbo = plot_map_hbo_og*epsilon(w,1);
    plot_map_hb = plot_map_hb_og*epsilon(w,2);
    plot_map = plot_map_hb + plot_map_hbo;

    noise_map = noise*max(plot_map,[],'all')*noise_map_og;
    plot_map = plot_map + noise_map;

    sensor_data = calculate_sensor_data_tissue_points(plot_map, Nx,dx,Ny,dy);
    reconstructed_image = reconstruct_k_wave_w_attenuation(sensor_data, Nx, Ny, dx, dy);
    %
    %{ 
    figure(f);
    subplot(2,2,count);
    imagesc(reconstructed_image);
    %}
    recon_holder_image(w,:,:) = reconstructed_image;
    count = count + 1;
end

C_nnls = calc_gen_nnls(recon_holder_image,epsilon, num_wavelength,num_type,Nx, Ny);

C_nnls(C_nnls < 0) = 1e-9;

% Calculate the total concentration
sum_C = squeeze(sum(C_nnls, 1));

% Calculate the hbO and hb Concentration
hbo_C = squeeze(C_nnls(1,:,:));
hb_C = squeeze(C_nnls(2,:,:));

% Get rid of zeros so that in the saturation calculation no dividing by zero
sum_C(sum_C == 0) = 1e-9;

% Calculate the saturation
hbo_S = hbo_C ./ sum_C;
hb_S = hb_C ./ sum_C;

%{
% Plot total concentration
subplot(3, 2, 1);
imagesc(sum_C);
colormap;
colorbar;
title('Total Concentration Map');
xlabel('X-axis');
ylabel('Y-axis');

% Plot HbO concentration
subplot(3, 2, 2);
imagesc(hbo_C);
colormap;
colorbar;
title('HbO Concentration Map');
xlabel('X-axis');
ylabel('Y-axis');

% Plot Hb concentration
subplot(3, 2, 3);
imagesc(hb_C);
colormap;
colorbar;
title('Hb Concentration Map');
xlabel('X-axis');
ylabel('Y-axis');

% Plot HbO saturation
subplot(3, 2, 4);
imagesc(hbo_S);
colormap;
colorbar;
title('HbO Saturation Map');
xlabel('X-axis');
ylabel('Y-axis');

% Plot Hb saturation
subplot(3, 2, 5);
imagesc(hb_S);
colormap;
colorbar;
title('Hb Saturation Map');
xlabel('X-axis');
ylabel('Y-axis');

% Apply mask where plot_map != 0
mask = plot_map ~= 0;
mask(mask ~= 0) = 1;
masked_hbo_S = hbo_S .* mask; % Apply the mask to hbo_S
masked_hb_S = hb_S .* mask;  % Apply the mask to hb_S

% Plotting the masked maps
figure;

% Plot masked Hb saturation
subplot(2, 1, 1);
imagesc(masked_hb_S);
colormap;
colorbar;
title('Masked Hb Saturation Map');
xlabel('X-axis');
ylabel('Y-axis');

% Plot masked HbO saturation
subplot(2, 1, 2);
imagesc(masked_hbo_S);
colormap;
colorbar;
title('Masked HbO Saturation Map');
xlabel('X-axis');
ylabel('Y-axis');

%}
%find_concentration_mask(sum_C, hb_S, hbo_S)

w_avg_hbo_S = sum(hbo_S.*(sum_C.^2))/sum(sum_C.^2);
w_avg_hb_S = sum(hb_S.*(sum_C.^2))/sum(sum_C.^2);

end