function [pressure_mask] = convert_mask_to_pressure(mask, epsilon, concentrations, concentration_names, wavelengths, Nx, Ny, plot_flag)

%Inputs
%Set of concentrations
%Epsilon where num_concentrations= length of concentrations
%Must ensure that epsilon and the concentrations are in the same order

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





