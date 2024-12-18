function sensor_data = calculate_sensor_data_tissue_points(plot_map,Nx,dx,Ny,dy)


PML_size = 20;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1540;	% [m/s]
%medium.alpha_coeff = 0.75;       % absorption coefficient [dB/(MHz^y cm)]
%medium.alpha_power = 0.6;        % absorption power y (typically 1–2)


% create initial pressure distribution using makeDisc

p0 = smooth(plot_map ,true);
% assign to the source structure
source.p0 = p0;

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false, 'PlotSim', false};


% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

end