function sensor_data = calculate_sensor_data_v2(concentration, epsilon,Nx,dx,Ny,dy)
PML_size = 20;    

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1540;	% [m/s]

% create initial pressure distribution using makeDisc
disc_magnitude_hbo = epsilon(1) * concentration;         % [Pa]
disc_magnitude_hb = epsilon(2)* (1-concentration);
disc_magnitude = disc_magnitude_hbo + disc_magnitude_hb;

disc_y_pos = Ny/2;            % [grid points]
disc_x_pos = Nx/2;           % [grid points]
disc_radius = 4;            % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

p0 = smooth(disc_1,true);
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