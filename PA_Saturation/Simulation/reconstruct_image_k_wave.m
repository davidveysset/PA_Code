%> @file reconstruct_image_k_wave.m
%> @brief Reconstructs the original pressure distribution from sensor data using K-Wave.
%> @author Calvin Smith
%> @date 12-20-24
%> @details This function uses K-Wave to reconstruct the original pressure distribution from sensor data 
%> captured by a linear transducer. The reconstruction process may introduce unavoidable artifacts and lobes.
%>
%> @param sensor_data Sensor data from a linear transducer.
%> @param Nx Number of X pixels (horizontal resolution).
%> @param Ny Number of Y pixels (vertical resolution).
%> @param dx Size of each X pixel in meters.
%> @param dy Size of each Y pixel in meters.
%>
%> @return p_xy_rs The reconstructed initial pressure distribution at t = 0 based on the sensor data.

function p_xy_rs = reconstruct_image_k_wave(sensor_data, Nx, Ny, dx, dy)
% DESCRIPTION:
% K-Wave reconstructs from the sensor data the original pressure
% distribution. There are some artifacts and lobes that are unavoidable.
%
% INPUTS:
% sensor data - sensor data from a linear transducer 
% Nx - Number of X Pixels
% Ny - Number of Y Pixels
% dx - size of each x pixel in meters
% dy - size of each y pixel in meters
%
% OUTPUTS:
% p_xy_rs - reconstructed initial pressure distribution at t= 0 based on
% the sensor data

% build kgrid
kgrid = kWaveGrid(Nx, dx, Ny, dy);
PML_size = 20;

% define the properties of the propagation medium
medium.sound_speed = 1540;           % [m/s]


% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

input_args = {'PMLInside', false, 'PMLSize', PML_size,0, 'Smooth', false, 'PlotPML', false};
% create the time array
kgrid.makeTime(medium.sound_speed);


% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% reconstruct the initial pressure
p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed, ...
     'Plot', false, 'PosCond', true, 'Interp', '*linear');

% define a second k-space grid using the dimensions of p_xy
[Nx_recon, Ny_recon] = size(p_xy);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy);

% resample p_xy to be the same size as source.p0
p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));


end