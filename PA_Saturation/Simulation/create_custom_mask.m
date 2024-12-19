function [pressure_mask] = convert_mask_to_pressure(mask, epsilon, concentrations Nx, Ny)

%Inputs
%Set of concentrations 
%Epsilon where num_concentrations= length of concentrations
%Must ensure that epsilon and the concentrations are in the same order

num_concentration = length(concentrations);

if num_concentration == length(epsilon)
    pressure_mask = zeros(Nx,Ny);
    for c = 1:num_concentration
        pressure_mask = pressure_mask + concentrations(c)*mask;
    end
else
    disp('Not enough Concentrations!');
end


end
