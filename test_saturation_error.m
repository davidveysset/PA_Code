data = load('example_recon_data_770_780.mat');
data = data.recon_noise_holder;
[C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error] = calc_error_saturation(data, E, wavelengths, expected_values, plot_hbo, plot_hb, plot_analysis);
data_770_780 = {C_nnls, hbo_C, hb_C, sum_C, hbo_S, hb_S, circle_concentration,saturation_error};