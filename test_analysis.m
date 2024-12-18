results_770_780 = load('unmixing_results_770_780_v1.mat','unmixing_results_770_780');
results_770_780 = results_770_780.unmixing_results_770_780;

results_750_850 = load('unmixing_results_750_850_v1.mat','unmixing_results_750_850');
results_750_850 = results_750_850.unmixing_results_750_850;

results_all = load('unmixing_results_all_v1.mat','unmixing_results_all');

results_all = results_all.unmixing_results_all;

saturation_error_770_780 = results_770_780{6};
saturation_error_750_850 = results_750_850{6};
saturation_error_all = results_all{6};

wavelengths = {'750-850', '770-780', '750-770-780-850'};

data = {results_750_850, results_770_780, results_all};

%multi_wavelength_analysis(data,wavelengths);

%data = load('bad_full_run_data.mat');
%data = data.total_data;

multi_wavelength_analysis(total_data,wavelengths);