function plot_grid( filename )
% This is used to plot the grid that is being written for 2D cfd
% calculations

M = csvread([filename '.bkcfd']);
x_node_loc = M(:,1:4);
y_node_loc = M(:,5:8);

cornerlocs_x = x_node_loc;
cornerlocs_y = y_node_loc;

plot(cornerlocs_x, cornerlocs_y,'b.')

end

