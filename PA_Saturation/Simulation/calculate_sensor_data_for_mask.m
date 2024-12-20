%> @file calculate_sensor_data_for_mask.m
%> @brief Generates sensor data using K-Wave simulation from a pressure mask and linear transducer.
%> @author Calvin Smith
%> @date 12-20-24
%> @details This function uses K-Wave to generate sensor data from a given pressure mask. 
%> The simulation assumes a linear transducer along the x-axis with `Nx` transducers and a pitch of `dx`. 
%> The output is a calculated array where:
%> - The y-axis represents the transducers.
%> - The x-axis represents a time stream of measurements.
%>
%> @param pressure_mask The pressure mask, calculated based on the tissue mask, concentrations, and wavelengths.
%> @param Nx Number of X pixels (horizontal resolution).
%> @param Ny Number of Y pixels (vertical resolution).
%> @param dx Size of each X pixel in meters (transducer pitch).
%> @param dy Size of each Y pixel in meters.
%>
%> @return sensor_data An array of shape (num_transducers, num_time) representing the data stream 
%> from the linear array of transducers. Additional customization for time measurement is planned for future versions.


function sensor_data = calculate_sensor_data_for_mask(pressure_mask,Nx,dx,Ny,dy)
% DESCRIPTION:
% K-Wave generates sensor data from the pressure mask along the x-axis with a linear transducer of
% Nx transducers and dx pitch. Returns a calculated array with transducers
% on the y-axis and the x-axis being a timme stream. 
%
% INPUTS:
% pressure_mask - Calculated based on the tissue mask, concentrations and wavelengths. 
% Nx - Number of X Pixels
% Ny - Number of Y Pixels
% dx - size of each x pixel in meters
% dy - size of each y pixel in meters
%
% OUTPUTS:
% sensor_data - array of (num_transducers, num_time). This represents the data stream from the 
% linear array of transducers. More customization will be added soon to
% increase/decrease time measurement.

% Generates boundary around the grid to prevent wrapping
PML_size = 20;
% Builds kgrid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1540;	% [m/s]
%medium.alpha_coeff = 0.75;       % absorption coefficient [dB/(MHz^y cm)]
%medium.alpha_power = 0.6;        % absorption power y (typically 1–2)


% create initial pressure
p0 = smooth(pressure_mask ,true);
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