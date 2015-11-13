function plot_solution_contours( filename, equal, contour_lines, num_divisions )
% 
% plot_solution_contours( filename, equal, contour_lines, num_divisions )
%
% This plots the solution from a *.obkcfd output file
% 
% INPUTS
% filename:         'filename' is the name of the *.obkcfd file, without the 
%                   .txt extension, that you want to plot
%
% equal:            0 or 1 (false or true) for plotting with equal axes
%
% contour_lines:    0 or 1 (false or true) for plotting with lines
%                   separating the contours 
%
% num_divisions:    number of contour levels for the contour plot
%
%
%

close all
gamma = 1.4;
Cv = 287/0.4;
colormap_color = 'jet';

M = csvread([filename '.obkcfd'],1,0);
x_node_loc = M(:,2:5);
y_node_loc = M(:,6:9);
rho_long = M(:,15);
rhou = M(:,16);
rhov = M(:,17);
E = M(:,18);

num_cells_x = length(unique(x_node_loc)) - 1;
num_cells_y = (length(M)+4)/num_cells_x;

if (num_cells_x * num_cells_y) ~= (length(M) + 4)
    fprintf('\n\nError in determining number of x and y cells, exiting...\n\n')
    return
end


list_number = 1;
for ix = 1:num_cells_x
    for iy = 1:num_cells_y
        if (ix == 1 || ix == num_cells_x) && (iy == 1 || iy == num_cells_y)
            continue
        elseif (ix == 1 || ix == num_cells_x || iy == 1 || iy == num_cells_y)
            list_number = list_number + 1;
            continue
        end
        
        x_cellcenter(ix,iy) = mean(x_node_loc(list_number,:));
        y_cellcenter(ix,iy) = mean(y_node_loc(list_number,:));

        vel_mag(ix,iy) = sqrt((rhou(list_number)/rho_long(list_number))^2 + (rhov(list_number)/rho_long(list_number))^2);
        u(ix,iy) = rhou(list_number)/rho_long(list_number);
        v(ix,iy) = rhov(list_number)/rho_long(list_number);
        rho(ix,iy) = rho_long(list_number);
        P(ix,iy) = (gamma-1)*(E(list_number) - rho_long(list_number)*(u(ix,iy)^2 + v(ix,iy)^2)/2);
        T(ix,iy) = (E(list_number) - rho_long(list_number)*vel_mag(ix,iy)^2/2)/Cv;

        list_number = list_number + 1;
    end
end

x_cellcenter = x_cellcenter(2:end,2:end);
x_cellcenter = x_cellcenter';
y_cellcenter = y_cellcenter(2:end,2:end);
y_cellcenter = y_cellcenter';
vel_mag = vel_mag(2:end,2:end);
vel_mag = vel_mag';
u = u(2:end,2:end);
u = u';
v = v(2:end,2:end);
v = v';
rho = rho(2:end,2:end);
rho = rho';
P = P(2:end,2:end);
P = P';
T = T(2:end,2:end);
T = T';


figure(1);
[q conthandle] = contourf(x_cellcenter, y_cellcenter, P, num_divisions);
colorbar;
colormap jet;

xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Pressure')
if (equal)
    axis equal
end
if ~contour_lines
    set(conthandle, 'LineStyle','none')
end


figure(2);
[q conthandle] = contourf(x_cellcenter, y_cellcenter, rho, num_divisions);
colorbar;
colormap jet;

set(conthandle, 'LineStyle','none')
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Density')
if (equal)
    axis equal
end
if ~contour_lines
    set(conthandle, 'LineStyle','none')
end

figure(3);
[q conthandle] = contourf(x_cellcenter, y_cellcenter, vel_mag, num_divisions);
colorbar;
colormap jet;

set(conthandle, 'LineStyle','none')
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Velocity Magnitude (m/s)')
if (equal)
    axis equal
end
if ~contour_lines
    set(conthandle, 'LineStyle','none')
end

figure(4);
[q conthandle] = contourf(x_cellcenter, y_cellcenter, u, num_divisions);
colorbar;
colormap jet;

set(conthandle, 'LineStyle','none')
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('X-Velocity (m/s)')
if (equal)
    axis equal
end
if ~contour_lines
    set(conthandle, 'LineStyle','none')
end

figure(5);
[q conthandle] = contourf(x_cellcenter, y_cellcenter, v, num_divisions);
colorbar;
colormap jet;

set(conthandle, 'LineStyle','none')
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Y-Velocity (m/s)')
if (equal)
    axis equal
end
if ~contour_lines
    set(conthandle, 'LineStyle','none')
end

figure(6);
[q conthandle] = contourf(x_cellcenter, y_cellcenter, T, num_divisions);
colorbar;
colormap jet;

set(conthandle, 'LineStyle','none')
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Temperature')
if (equal)
    axis equal
end
if ~contour_lines
    set(conthandle, 'LineStyle','none')
end

end