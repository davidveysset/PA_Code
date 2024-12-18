function p_xy_rs = reconstruct_image_k_wave_w_attenuation(sensor_data, Nx, Ny, dx, dy)

kgrid = kWaveGrid(Nx, dx, Ny, dy);

PML_size = 20;

% define the properties of the propagation medium
medium.sound_speed = 1540;           % [m/s]


% create the time array
kgrid.makeTime(medium.sound_speed);
% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;


input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false, 'PlotSim', false};

source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

medium.alpha_sign = [1,1];
% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Restore original medium absorption sign for subsequent calculations

% add first order compensation for only recording over a half plane
p0_recon = 2 * p0_recon;

% repeat the FFT reconstruction for comparison
p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed, ...
    'PosCond', true, 'Interp', '*linear');

% define a second k-space grid using the dimensions of p_xy
[Nx_recon, Ny_recon] = size(p_xy);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy);

% resample p_xy to be the same size as source.p0
p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));

%
%{
% p_test = p_xy_rs;
for x = 1:Nx
    for y = 1:Ny
        point = [x,y];
        avg_r = 0;
        for k = 1:Nx
            transducer = [k,0];
            r = norm(point - transducer);
            avg_r = avg_r + r;            
        end
     avg_r = avg_r/Nx;
     p_test(x,y) = p_xy_rs(x,y) * 1/(atan(200/avg_r))^2;
    end
end
%}
%{
figure;
imagesc(p_xy_rs);
colormap;
colorbar;
%}
end