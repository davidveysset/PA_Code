%> @file calc_gen_nnls.m
%> @brief Performs spectral unmixing using non-negative least squares optimization.
%> @author Calvin Smith
%> @date 12-20-24
%> @details This function takes reconstructed pressure data and absorption coefficients to solve 
%> the equation \f$ P(w,x,y) = \epsilon(w,t) * C(t,x,y) \f$ using a non-negative least squares (NNLS) 
%> optimization. It minimizes the objective function:
%> \f$ \text{Argmin} \{ P - EC \} \f$
%> where:
%> - \f$ P \f$: Pressure (reconstructed data).
%> - \f$ E \f$: Absorption coefficients by wavelength and type.
%> - \f$ C \f$: Concentration.
%> The function finds the optimal \f$ C(t,x,y) \f$ that fits best to the data.
%>
%> @param data Reconstructed data with shape (num_wavelength, Nx, Ny).
%> @param E A (num_wavelength, num_types) matrix of absorption coefficients. Must match the order in concentrations.
%> @param num_type Number of types (e.g., biomarkers or tissue components).
%> @param Nx Number of X pixels.
%> @param Ny Number of Y pixels.
%>
%> @return C_nnls The NNLS solution to \f$ P = EC \f$, with shape (num_type, Nx, Ny).

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