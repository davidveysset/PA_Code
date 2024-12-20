function [pressure_mask] = convert_mask_to_pressure(mask, epsilon, concentrations, concentration_names, wavelengths, Nx, Ny, plot_flag)
% DESCRIPTION:
% Converts boolean mask, concentrations set, and the absorption coefficient
% matrix into pressure masks that will be converted in to sensor data. 
% Calculate pressure using P(w, x, y) = Epsilon(w, t) * C(type, x, y). 
% P = Pressure, Epsilon = Absorption Coefficients by wavelength and type, C = Concentration
%
% INPUTS:
% mask  -  a boolean mask of size(Nx, Ny) that is the shape of the tissue to simulate 
% epsilon - a (num_wavelength, num_types) matrix of the absorption coefficients. must be in order that the types are specified in concentrations
% wavelengths - 1D array of the wavelengths, ex:[770,780], this must be the same size as epsilon
% concentrations - an array of size (unique_c*num_iter, num_types).concentrations defines the relative saturations of each type in the generated pressure mask
% type_names - a "CELL" array of names, ex : {'Hb', 'HbO'}. The order is important must match with epsilon and concentrations
% Nx - Number of X Pixels
% Ny - Number of Y Pixels
% plotflag -  boolean, if true, plots each pressure mask for all num_concentrations and wavelengths. Only use for testing when num_iter = 1!
%
% OUTPUTS:
% pressure_mask = calculates the pressure using P(w,x,y) = epsilon(w,t) * C(t,x,y)
% This is used for calculating the sensor data



num_concentrations = size(concentrations,1);
num_type = size(concentrations,2);
num_wavelength = size(epsilon,1);


if num_type == size(epsilon, 2)
    pressure_mask = zeros(num_wavelength,num_concentrations,Nx,Ny);
    for c = 1: num_concentrations
        for w = 1:num_wavelength
            for t = 1:num_type
                c_by_type_mask =  epsilon(w,t)*concentrations(c,t)*mask;
                pressure_mask(w,c,:,:) = squeeze(pressure_mask(w, c, :,:)) + c_by_type_mask;
            end
        end
    end
else
    disp('Not enough Concentrations!');
end

if plot_flag
    for w = 1:num_wavelength
        figure;
        max_pressure = max(pressure_mask(w,:,:,:),[],'all');
        count = 1;
        for c = 1:num_concentrations
            for t = 1:num_type
                subplot(num_concentrations, num_type, count)
                imagesc(epsilon(w,t)*concentrations(c,t)*mask);
                colormap;
                colorbar;
                clim([0,max_pressure]);
                title(sprintf('Inital Pressure Distribution for %s at %d nm', concentration_names{t}, wavelengths(w)) );
                count = count + 1;
            end
        end
    end
    for w = 1:num_wavelength
        figure;
        max_pressure = max(pressure_mask(w,:,:,:),[],'all');
        for c = 1:num_concentrations
            subplot(1, num_concentrations, c)
            imagesc(squeeze(pressure_mask(w,c,:,:)));
            colormap;
            colorbar;
            clim([0,max_pressure]);
            title(sprintf('Inital Pressure Distribution for Combined Types at %d', wavelengths(w)) );
           
        end
    end
end
end





