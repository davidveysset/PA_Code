data = load('4_wavelength_set_data.mat','total_data');
total_data = data.total_data;
wavelength_names = {'750-850', '770-780', '780-1030', '750-770-780-850'};

multi_wavelength_analysis(total_data,wavelength_names);