
num_x_pixels = 200;
num_y_pixels = 200;
dx = 1e-6;
dy = 1e-6;

sensor_data = calculate_sensor_data(1, 1200, num_x_pixels, dx, num_y_pixels, dy);
                 
noise = randn(200,943);

noise_mag = max(sensor_data,[],'all');
noisy_sensor_data_20 = noise_mag * 0.05 * noise + sensor_data;
noisy_sensor_data_2 = noise_mag * 0.5 * noise + sensor_data;

nsd_slice_20 = squeeze(noisy_sensor_data_20(100,:)) -150;
nsd_slice_2 = squeeze(noisy_sensor_data_2(100,:)) +150;

figure;
hold on;
plot(nsd_slice_20, 'DisplayName','SNR 20');
plot(nsd_slice_2, 'DisplayName','SNR 2')
%title('Time Slice Data for SNR of 2 and 20');
xlabel('1e-9 seconds');
ylabel('Signal');
ax = gca; % Get current axes handle
ax.FontSize = 20; % Increase tick label font size

legend('show','FontSize',14);
grid on;
