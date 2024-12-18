function [full_map, full_map_hb, full_map_hbo] = create_tissue_types_mask(epsilon, circle_radius, hbo_C, c_std, fill_pct, full_Nx, full_Ny, Nx, Ny)

% Parameters:
% Filling Percentage
% Circle radius
% HbO Concentration
% Standard Deviation

% Calculate the total area
% While loop that picks points and if not touch another circle places a
% circle until the counted area = filling percentage * total area
% Each circle is given by the hb0 concentration

if nargin < 1
    epsilon = [1405,518; 691,1050];
end
if nargin < 2
    circle_radius = 12;
end
if nargin < 3
    hbo_C = 0.5;
end
if nargin < 4
    c_std = 0.05;
end
if nargin < 5
    fill_pct = 0.1;
end
if nargin < 6
    full_Nx = 400;
end
if nargin < 7
    full_Ny = 200;
end
if nargin < 8
    Nx = 100;
end
if nargin < 9
    Ny = 100;
end

total_area = Nx * Ny;
disc_radius = circle_radius;
loop_flag = true;
plot_map_hb = zeros(size(epsilon,1),Nx,Ny);
plot_map_hbo = zeros(size(epsilon,1),Nx,Ny);
plot_map = zeros(size(epsilon,1),Nx,Ny);

map_area = 0;
circle_area = pi * circle_radius^2;
bool_map = zeros(Nx, Ny);
count = 0;
disp(total_area);

while loop_flag
    count = count + 1;
    x_pos = round((Nx - 2*circle_radius) * rand) + 2*circle_radius;
    y_pos = round((Nx - 2*circle_radius) * rand) + 2*circle_radius;

    good_point = true;

    % Check if the new circle overlaps with existing circles
    for i = -circle_radius:circle_radius
        for j = -circle_radius:circle_radius
            % Ensure indices are within bounds
            x_index = x_pos + i;
            y_index = y_pos + j;
            if x_index > 0 && x_index <= Nx && y_index > 0 && y_index <= Ny
                if bool_map(x_index, y_index) == 1
                    good_point = false;
                    break; 
                end
            end
        end
        if ~good_point
            break; 
        end
    end

    if good_point
        hbo_total = hbo_C + ((2*c_std)*rand - c_std);
        hb_total = 1 - hbo_total;
        hbo_disc = hbo_total * epsilon(1) * makeDisc(Nx, Ny, x_pos, y_pos, disc_radius);
        hb_disc =  hb_total * epsilon(2) * makeDisc(Nx, Ny, x_pos, y_pos, disc_radius);
        total_disc = hbo_disc + hb_disc;
        map_area = map_area + circle_area;
        plot_map = plot_map + total_disc;
        plot_map_hb = plot_map_hb + hb_disc;
        plot_map_hbo = plot_map_hbo + hbo_disc;


        % Mark the circle area in the boolean map
        for i = -circle_radius:circle_radius
            for j = -circle_radius:circle_radius
                x_index = x_pos + i;
                y_index = y_pos + j;
                if x_index > 0 && x_index <= Nx && y_index > 0 && y_index <= Ny
                    bool_map(x_index, y_index) = 1;
                end
            end
        end
    end

    disp(map_area);

    if map_area > total_area * fill_pct || count > 1000
        loop_flag = false;
    end
end

full_map = zeros(full_Nx, full_Ny);
full_map_hb = zeros(full_Nx, full_Ny);
full_map_hbo = zeros(full_Nx, full_Ny);
%{
full_map(full_Nx/2 - Nx/2 +1: full_Nx/2 + Nx/2, full_Ny/2 - Ny/2 +1 : full_Ny/2 + Ny/2) = plot_map;
full_map_hb(full_Nx/2 - Nx/2 +1: full_Nx/2 + Nx/2, full_Ny/2 - Ny/2 +1 : full_Ny/2 + Ny/2) = plot_map_hb;
full_map_hbo(full_Nx/2 - Nx/2 +1: full_Nx/2 + Nx/2, full_Ny/2 - Ny/2 +1 : full_Ny/2 + Ny/2) = plot_map_hbo;
%}

center_x = 100;
center_y = 100;
x_start = round(center_x - Nx / 2) + 1;
x_end = x_start + Nx - 1;
y_start = round(center_y - Ny / 2) + 1;
y_end = y_start + Ny - 1;
% Place the smaller maps into the larger maps at the specified location
full_map(x_start:x_end, y_start:y_end) = plot_map;
full_map_hb(x_start:x_end, y_start:y_end) = plot_map_hb;
full_map_hbo(x_start:x_end, y_start:y_end) = plot_map_hbo;


% Assuming full_map, full_map_hb, full_map_hbo, full_map_hbo_S, and full_map_hb_S are defined

figure;

% Full map of tissue
subplot(1, 3, 1);
imagesc(full_map);
title('Full Tissue Map');
xlabel('X');
ylabel('Y');
colorbar;

% Full map of Hb
subplot(1, 3, 2);
imagesc(full_map_hb);
title('Full Hb Map');
xlabel('X');
ylabel('Y');
colorbar;

% Full map of HbO
subplot(1, 3, 3);
imagesc(full_map_hbo);
title('Full HbO Map');
xlabel('X');
ylabel('Y');
colorbar;

% Adjust the layout for better visualization
sgtitle('Full Tissue and Concentration Maps');

end
