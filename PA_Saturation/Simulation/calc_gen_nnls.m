function C_nnls = calc_gen_nnls(data, E, num_wavelength, num_type, Nx, Ny)
%DESCRIPTION:
% Takes in reconstructed data which is pressure maps by wavelength and then
% uses P(w,x,y) = epsilon(w,t)*C(t,x,y) and uses calc_gen_nnls which
% performs a non-negative least squares optimization on Arg_min{P-EC} since
% P-EC = 0, to find C(t,x,y) that fits best. 
% INPUTS:
% data - reconstructed data with shape(num_wavelength, Nx, Ny)
% E - a (num_wavelength, num_types) matrix of the absorption coefficients. must be in order that the types are specified in concentrations
% num_type - number of types
% Nx - Number of X Pixels
% Ny - Number of Y Pixels
%
% OUTPUTS:
% C_nnls = non-negative least squares solution to P = EC where C_nnls has
% shape(num_type, Nx, Ny)

%define pressures
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