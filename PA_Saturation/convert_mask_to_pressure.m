function [pressure_mask] = convert_mask_to_pressure(mask, epsilon, concentrations Nx, Ny)

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
        for w = 1:num_wavelengths
            for t = 1:num_type
                pressure_mask(w,c,:,:) = pressure_mask(w, c, :,:) + concentrations(c,t)*mask;
            end
        end
    end
else
    disp('Not enough Concentrations!');
end


end
