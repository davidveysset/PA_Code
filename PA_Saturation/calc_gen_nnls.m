function C_nnls = calc_gen_nnls(data, E, num_wavelength, num_type, Nx, Ny)
P = zeros(num_wavelength, Nx, Ny);

% C_nnls is concentrations from non-negative least squares as a function of
% (type, x, y)
C_nnls = zeros(num_type, Nx, Ny);

% Change the convergence since it was causing error
opts = optimset('lsqnonneg');
opts.TolX = 1e-9;

% Loop through each noise level, concentration, and optimize Arg_Min{P - EC}
% and Store in C_nnls Matrix

for x = 1:Nx
    for y = 1:Ny
        P_point = squeeze(data(:, x, y));
        C_point = lsqnonneg(E, P_point, opts);
        C_nnls(:, x, y) = C_point;

    end
end
end