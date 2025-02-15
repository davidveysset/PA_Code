function snr_data_holder = snr_compute(noisy_sensor_data, wavelengths, Nx, noise_levels, noise_strength)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
num_wavelength = length(wavelengths);
signal_data_holder = zeros(num_wavelength,1);

for w = 1:num_wavelength
   signal = max(squeeze(noisy_sensor_data(1,w,1,round(Nx/2),:)));
   signal_data_holder(w) = signal;
end

signal_min = min(signal_data_holder(:));
snr_data_holder = signal_min./(noise_levels.*noise_strength);

end