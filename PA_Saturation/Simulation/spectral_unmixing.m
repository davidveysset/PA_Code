%> @file spectral_unmixing.m
%> @brief Analyzes reconstructed pressure data to calculate concentrations, saturations, and weighted averages.
%> @author Calvin Smith
%> @date 12-20-24
%> @details This function takes reconstructed pressure data (maps by wavelength) and performs spectral unmixing. 
%> Using the equation P(w, x, y) = Epsilon(w, t) * C(t, x, y), it applies `calc_gen_nnls` to perform non-negative 
%> least squares (NNLS) optimization to minimize Argmin{P - EC}, where:
%> - P: Reconstructed pressure data.
%> - Epsilon: Absorption coefficients.
%> - C: Concentration map.
%>
%> Additional Calculations:
%> - **Sum by Type**:
%>   - Computes `sum_C` by summing concentrations by type to create a map of size (Nx, Ny).
%> - **Saturation Maps**:
%>   - Divides each concentration map by `sum_C` to calculate saturation maps.
%> - **Weighted Average**:
%>   - Computes `w_avg_by_type` by weighting each point with `sum_C`, dividing by the scalar sum of `sum_C`, 
%>     and aggregating the result into a scalar value.
%>
%> @param epsilon A (num_wavelength, num_types) matrix of absorption coefficients. Must match the order in `type_names`.
%> @param recon_data_w Reconstructed data with shape (num_wavelength, Nx, Ny).
%> @param type_names A cell array of type names (e.g., {'Hb', 'HbO'}). Must match the order in `epsilon`.
%> @param Nx Number of X pixels (horizontal resolution).
%> @param Ny Number of Y pixels (vertical resolution).
%>
%> @return sum_C A map of size (Nx, Ny) representing summed concentrations by type.
%> @return saturations_by_type An array of size (num_type, Nx, Ny) containing saturation maps for each type.
%> @return concentrations_by_type An array of size (num_type, Nx, Ny) containing concentration maps for each type.
%> @return w_avg_by_type A 1D array of size (num_type) containing the weighted average values for each type.

function [sum_C, saturations_by_type, concentrations_by_type, w_avg_by_type] = spectral_unmixing(epsilon, recon_data_w, ...
    type_names, Nx, Ny)
% DESCRIPTION:
% Takes in reconstructed data which is pressure maps by wavelength and then
% uses P(w,x,y) = epsilon(w,t)*C(t,x,y) and uses calc_gen_nnls which
% performs a non-negative least squares optimization on Arg_min{P-EC} since
% P-EC = 0, to find C(t,x,y) that fits best. 
% Calculates the summing by types for each point to find sum_C and then
% dividing concentration/sum_C
% Also calculates w_avg_by_type by multiplying each point by sum_C and
% then dividing the scalar sum(sum_C). It then sums this total to provide a
% scalar metric of the "weighted average" over the whole map.
% 
% INPUTS:
% epsilon - a (num_wavelength, num_types) matrix of the absorption coefficients. must be in order that the types are specified in concentrations
% recon_data_w - reconstructed data with shape(num_wavelength, Nx, Ny)
% type_names - a "CELL" array of names, ex : {'Hb', 'HbO'}. The order is important must match with epsilon and concentrations
% Nx - Number of X Pixels
% Ny - Number of Y Pixels
%
% OUTPUTS:
% sum_C - summing by type to produce a map(Nx, Ny) of the summed concentrations
% saturations_by_type - array of size(num_type, Nx, Ny) of the saturation maps
% concentrations_by_type - array of size(num_type, Nx, Ny) of the concentration maps
% w_avg_by_type - 1D of size(num_type) which stores the weighted average values

num_wavelength = size(epsilon,1);
num_type = size(epsilon,2);

concentrations_by_type = zeros(num_type,Nx,Ny);
concentrations_by_type = calc_gen_nnls(recon_data_w, epsilon, num_wavelength,num_type,Nx,Ny);
% remove errors so that saturation calculation doesn't have NAN 
concentrations_by_type(concentrations_by_type < 0) = 1e-9;

% Calculate the total concentration
sum_C = squeeze(sum(concentrations_by_type, 1));
sum_C = reshape(sum_C, 1, Nx, Ny); 
sum_C(sum_C == 0) = 1e-9;

% Calculate the saturation
saturations_by_type = concentrations_by_type./ sum_C;

%Calculate the weighted average
w_avg_by_type = squeeze(sum(saturations_by_type.*(sum_C)))/(sum(sum_C,'all'));
end