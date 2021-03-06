function Plot_Solution( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all
gamma = 1.4;

M = csvread(filename,1,0);

%x_node_loc = M(:,1:4);
%y_node_loc = M(:,5:8);
%rho = M(:,9);
%rhou = M(:,10);
%rhov = M(:,11);
%rhoE = M(:,12);

x_node_loc = M(:,2:5);
y_node_loc = M(:,6:9);
rho = M(:,15);
rhou = M(:,16);
rhov = M(:,17);
rhoE = M(:,18);

for k = 1:size(x_node_loc,1)
    x_cellcenter(k) = mean(x_node_loc(k,:));
    y_cellcenter(k) = mean(y_node_loc(k,:));
end

vel_mag = sqrt((rhou./rho).^2+(rhov./rho).^2);
u = rhou./rho;
v = rhov./rho;

P = (gamma-1).*(rhoE - rho.*vel_mag.^2./2);

Cv = 287/0.4;
T = (rhoE - rho.*vel_mag.^2./2)./Cv;


figure(1);
scatter(x_cellcenter, y_cellcenter, 50, P, 'square', 'filled')
colorbar;
colormap(jet);
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Pressure')

figure(2);
scatter(x_cellcenter, y_cellcenter, 50, rho, 'square', 'filled')
colorbar;
colormap(jet);
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Density')

figure(3);
scatter(x_cellcenter, y_cellcenter, 50, vel_mag, 'square', 'filled')
colorbar;
colormap(jet);
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Velocity Magnitude (m/s)')

figure(4);
scatter(x_cellcenter, y_cellcenter, 50, u, 'square', 'filled')
colorbar;
colormap(jet);
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('X-Velocity (m/s)')

figure(5);
scatter(x_cellcenter, y_cellcenter, 50, v, 'square', 'filled')
colorbar;
colormap(jet);
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Y-Velocity (m/s)')

figure(6);
scatter(x_cellcenter, y_cellcenter, 50, T, 'square', 'filled')
colorbar;
colormap(jet);
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Temperature')

% figure(7);
% surf(x_cellcenter, y_cellcenter, P)


end