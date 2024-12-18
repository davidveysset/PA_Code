% Structure
% two wavelengths, two types hb, hbo
% construct the maps of the 9 points for hb and inverse for hbo
% Input also needs to consider the epsilon values, concentration_map,
% num_x_steps, num_y_steps
epsilon_770_780 = [1361 636; 1075 710]; 
epsilon_750_850 = [1405,518; 691,1050];
epsilon_780_1030 = [1075 710; 1024 206];
epsilon = epsilon_750_850;
num_wavelength = length(epsilon);
types = [1,2];
num_type = length(types);
Nx = 400;
Ny = 200;
dx = 5e-6;
dy = 5e-6;

recon_holder_image = zeros(num_wavelength,Nx, Ny);

f = figure;
hold on;

count = 1;
for w = 1:num_wavelength
    E = epsilon(w,:);
    [plot_map, plot_map_hb, plot_map_hbo]= create_tissue_types_mask(E);
    sensor_data_hb = calculate_sensor_data_tissue_points(plot_map_hb, Nx,dx,Ny,dy);
    sensor_data_hbo = calculate_sensor_data_tissue_points(plot_map_hbo, Nx,dx,Ny,dy);
    reconstructed_image_hb = reconstruct_k_wave_w_attenuation(sensor_data_hb, Nx, Ny, dx, dy);
    reconstructed_image_hbo = reconstruct_k_wave_w_attenuation(sensor_data_hbo, Nx, Ny, dx, dy);
    reconstructed_image  = reconstructed_image_hb + reconstructed_image_hbo;
    figure(f);
    subplot(2,2,count);
    imagesc(reconstructed_image);
    recon_holder_image(w,:,:) = reconstructed_image;
    count = count + 1;
end


%save('recon_tissue_circles.mat','recon_holder_image');

%data = load('recon_dense_circles.mat','recon_holder_image');
%data = data.recon_holder_image;
%data = squeeze(recon_holder_image);
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
figure;

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
