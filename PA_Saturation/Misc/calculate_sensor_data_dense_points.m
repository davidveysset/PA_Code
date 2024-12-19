function sensor_data = calculate_sensor_data_dense_points(epsilon,concentration_array, num_x_steps,num_y_steps,Nx,dx,Ny,dy)

if nargin < 1
    epsilon = [1405];
end
if nargin <2

    concentration_array = [0.9, 0.6,0.3;
        0.3,0.9,0.6;
        0.6,0.3,0.9];
end
if nargin <3
    num_x_steps = 3;
end
if nargin <4
    num_y_steps = 3;
end
if nargin< 5
    Nx = 200;
end
if nargin <6
    dx = 1e-6;
end
if nargin< 7
    Ny = 200;
end
if nargin <8
    dy = 1e-6;
end

PML_size = 20;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1540;	% [m/s]
medium.alpha_coeff = 0.75;       % absorption coefficient [dB/(MHz^y cm)]
medium.alpha_power = 0.6;        % absorption power y (typically 1–2)


% create initial pressure distribution using makeDisc
disc_magnitude = epsilon;       % [Pa]

x_step = round(Nx/(2*num_x_steps));
y_step = round(Ny/(2*num_y_steps));


plot_map = 0;
for n =1:num_x_steps
    for m =1:num_y_steps
        disc_x_pos = (2*n-1)* x_step;
        disc_y_pos = (2*m-1)* y_step;
        disc_radius = 4;
        concentration = concentration_array(n,m);
        disc = concentration*disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

        plot_map = plot_map + disc;
    end
end

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


medium.alpha_sign = [-1 1];
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


end