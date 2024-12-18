% Structure
% two wavelengths, two types hb, hbo
% construct the maps of the 9 points for hb and inverse for hbo
% Input also needs to consider the epsilon values, concentration_map,
% num_x_steps, num_y_steps

Nx = 200;
Ny = 200;
dx = 1e-6;
dy = 1e-6;
num_x_steps = 3;
num_y_steps = 3;
hbo_concentration_array = [0.9, 0.6,0.3;
    0.3,0.9,0.6;
    0.6,0.3,0.9];
hb_concentration_array = ones(size(hbo_concentration_array)) - hbo_concentration_array;


epsilon_750_850 = [1405,518;691,1050];
num_wavelength = length(epsilon_750_850);
types = [1,2];
num_type = length(types);

recon_holder_image = zeros(num_wavelength, num_type, Nx, Ny);

f = figure;
hold on;

count = 1;
for w = 1:num_wavelength
    for t = 1:num_type
        epsilon = epsilon_750_850(w,t);
        if t == 1
            concentration_array = hbo_concentration_array;
        elseif t==2
            concentration_array = hb_concentration_array;
        end
        sensor_data = calculate_sensor_data_dense_points(epsilon, concentration_array, num_x_steps, num_y_steps, Nx ,dx , Ny, dy);
        reconstructed_image = reconstruct_k_wave_w_attenuation(sensor_data, Nx, Ny, dx, dy);
        figure(f);
        subplot(2,2,count);
        imagesc(reconstructed_image);
        recon_holder_image(w,t, :,:) = reconstructed_image;
        count = count + 1;
    end
end

recon_holder_image = sum(recon_holder_image,2);

save('recon_dense_circles.mat','recon_holder_image');

%data = load('recon_dense_circles.mat','recon_holder_image');
%data = data.recon_holder_image;
data = squeeze(recon_holder_image);
C_nnls = calc_gen_nnls(data,epsilon_750_850, num_wavelength,num_type,Nx, Ny);


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
subplot(3,2,1)
imagesc(sum_C)
colormap;
colorbar;

subplot(3,2,2)
imagesc(hbo_C)
colormap;
colorbar;

subplot(3,2,3)
imagesc(hb_C)
colormap;
colorbar;


subplot(3,2,4)
imagesc(hbo_S)
colormap;
colorbar;

subplot(3,2,5)
imagesc(hb_S)
colormap;
colorbar;