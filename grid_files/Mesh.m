clear all
close all
clc

%% Define Boundaries

% What is the horizontal and vertical spacing?

num_cells_x = 150;
num_cells_y = 30;

% What is the rough size of domain
xmin = -0.1;
xmax = 0.6;
ymin = -0.1;
ymax = 0.2;

%%%%%%%%%%%%%% DO NOT EDIT THIS %%%%%%%%%%%%%%%%
x = linspace(xmin,xmax,5000);
y = linspace(ymin,ymax,5000);
syms x_sym y_sym
%%%%%%%%%%%%%% DO NOT EDIT THIS %%%%%%%%%%%%%%%%

% edge_state = 0: interior cell
% edge_state = 1: standard wall ghost cell, aka "wall"
% edge_state = 2: farfield zero gradient ghost cell, aka "pressure outlet"
% edge_state = 3: pressure inlet

% Define left line and boundary type
xL = 0;
boundary.left = 3;

% Define bottom line and boundary type
yB = 0;
yB_sym = 0;
boundary.bottom = 1;

% Define right line and boundary type
xR = 0.5;
boundary.right = 2;

% Define top line and boundary type
yT = -6.1363.*x.^4 + 5.0939.*x.^3 - 0.5857.*x.^2 - 0.1119.*x + 0.1023;
yT_sym = -6.1363.*x_sym.^4 + 5.0939.*x_sym.^3 - 0.5857.*x_sym.^2 - 0.1119.*x_sym + 0.1023;
boundary.top = 1;

% Define global initial conditions
global_initial.rho  = 2.5;
global_initial.rhou = 625;
global_initial.rhov = 0;
global_initial.rhoE = 1906600;

% The initial conditions on the "perturbed" boundary
perturbed_initial.rho  = 2.5;
perturbed_initial.rhou = 625;
perturbed_initial.rhov = 0;
perturbed_initial.rhoE = 1906600;

filename = 'half_CDnozzle_sonic_lowinletvel';

% Enter number of decimal places with which to round your node locations
round_xnodes = 4;
round_ynodes = 4;

% Do you want to plot the grid after it is defined?
plot_grid = 1;

%% Fix input

% Account for case of boundaries specified as scalar values - i.e. straight
% lines

if isscalar(xL)
    xL = xL*ones(size(y));
end
if isscalar(yB)
    yB = yB*ones(size(x));
end
if isscalar(xR)
    xR = xR*ones(size(y));
end
if isscalar(yT)
    yT = yT*ones(size(x));
end



%% Create the mesh of nodes

% Find intersects of each line to define corners of domain
for corner = 1:4
    x_BT = x;
    y_LR = y;
    if ismember(corner,[1 2])
        y_BT = yB;
    else
        y_BT = yT;
    end
    if ismember(corner,[1 4])
        x_LR = xL;
    else 
        x_LR = xR;
    end
    
    min_distance = 1E10;
    for index = 1:length(x)
        distance = sqrt((x_LR(index) - x_BT).^2 + (y_LR(index) - y_BT).^2); 
        if ~isempty(find(distance < min_distance, 1))
            min_distance = min(distance(distance < min_distance));
            x_intersect = x_BT(distance == min_distance);
            y_intersect = y_BT(distance == min_distance);
        else
            break
        end
    end
    
    if corner == 1
        x_BL = x_intersect;
        y_BL = y_intersect;
    elseif corner == 2
        x_BR = x_intersect;
        y_BR = y_intersect;
    elseif corner == 3
        x_TR = x_intersect;
        y_TR = y_intersect;
    elseif corner == 4
        x_TL = x_intersect;
        y_TL = y_intersect;
    end
end

% Equally space on top and bottom

dx_top = (x_TR - x_TL)/(num_cells_x + 1);
dx_bottom = (x_BR - x_BL)/(num_cells_x + 1);

x_top = linspace(x_TL, x_TR, num_cells_x + 1);
x_bottom = linspace(x_BL, x_BR, num_cells_x + 1);

y_top = eval(subs(yT_sym,x_sym,linspace(x_TL,x_TR,num_cells_x+1)));
y_bottom = eval(subs(yB_sym,x_sym,linspace(x_BL,x_BR,num_cells_x+1)));

for ix = 1:num_cells_x+1
    x_nodes(:,ix) = linspace(x_top(ix), x_bottom(ix), num_cells_y + 1);
    y_nodes(:,ix) = linspace(y_top(ix), y_bottom(ix), num_cells_y + 1);
end

x_nodes = [x_nodes(:,1) + (x_nodes(:,1)-x_nodes(:,2)), x_nodes, x_nodes(:,end) + (x_nodes(:,end)-x_nodes(:,end-1))];
y_nodes = [y_nodes(:,1) + (y_nodes(:,1)-y_nodes(:,2)), y_nodes, y_nodes(:,end) + (y_nodes(:,end)-y_nodes(:,end-1))];

x_nodes = [x_nodes(1,:) + (x_nodes(1,:)-x_nodes(2,:)); x_nodes; x_nodes(end,:) + (x_nodes(end,:)-x_nodes(end-1,:))];
y_nodes = [y_nodes(1,:) + (y_nodes(1,:)-y_nodes(2,:)); y_nodes; y_nodes(end,:) + (y_nodes(end,:)-y_nodes(end-1,:))];


% Set the corners to not have a value - there is no ghost cell on corners!
x_nodes(1,1) = nan;
x_nodes(1,end) = nan;
x_nodes(end,1) = nan;
x_nodes(end,end) = nan;
y_nodes(isnan(x_nodes)) = nan;

x_nodes = round(x_nodes.*(10^round_xnodes))./(10^round_xnodes);

y_nodes = round(y_nodes.*(10^round_ynodes))./(10^round_ynodes);



%% Set grid values

cellnum = 1;
for pointnum = 1:numel(y_nodes)
    row = mod(pointnum-1,(num_cells_y + 3)) + 1;
    column = floor(pointnum/(num_cells_y + 3)) + 1;
    if (row == num_cells_y + 3) || (column == num_cells_x + 3)
        continue % Continue for bottom row and right column of nodes
    end
    
    cornerlocs_x = [x_nodes(row+1,column), x_nodes(row+1,column+1), x_nodes(row,column+1), x_nodes(row,column)];
    cornerlocs_y = [y_nodes(row+1,column), y_nodes(row+1,column+1), y_nodes(row,column+1), y_nodes(row,column)];
    if isnan(sum(cornerlocs_x)) %if one of the nodes is a nan, it's that corner ghost cell we don't want
        continue
    else
%------- ASSIGN CELL NUMBER -----------------------------------------------         
        grid.cellnumber(cellnum,1) = cellnum;
        
        % Assign cornerlocs_x and cornerlocs_y
        grid.cornerlocs_x(cellnum,:) = cornerlocs_x;
        grid.cornerlocs_y(cellnum,:) = cornerlocs_y;
             
%------- ASSIGN ADJACENT CELLS --------------------------------------------         
        % Note: at end, need to make all cells < 0 = to -2!
        if cellnum == 19
            disp('')
        end
        % Left
        if (column == 2)
            grid.adjacent_cells(cellnum,1) = cellnum - (num_cells_y + 1);
            if (row == num_cells_y + 2) % Leftmost ghost cell on bottom
                grid.adjacent_cells(cellnum,1) = -2;
            end
        elseif (column == num_cells_x + 2) %Right ghost cells
            grid.adjacent_cells(cellnum,1) = cellnum - (num_cells_y + 1);
        else
            grid.adjacent_cells(cellnum,1) = cellnum - (num_cells_y + 2);
        end
        
        % Bottom
        if (row == num_cells_y + 1) && (column == 1)
            grid.adjacent_cells(cellnum,2) = -2;
        elseif (row == num_cells_y + 2)
            grid.adjacent_cells(cellnum,2) = -2;           
        else
            grid.adjacent_cells(cellnum,2) = cellnum + 1;
        end
        
        % Right CHANGED FROM -2 TO +2
        if (column == num_cells_x + 1) || (column == 1) % Rightmost column that we are working with for nodes or leftmost column
            grid.adjacent_cells(cellnum,3) = cellnum + (num_cells_y + 1);
            if (row == 1) % rightmost ghost cell on top edge
                grid.adjacent_cells(cellnum,3) = -2;
            end
        else
            grid.adjacent_cells(cellnum,3) = cellnum + (num_cells_y + 2);
        end
        
        % Top
        if (row == 1) || (row == 2 && column == num_cells_x + 2)
            grid.adjacent_cells(cellnum,4) = -2;
        else
            grid.adjacent_cells(cellnum,4) = cellnum - 1;
        end
        
        
%------- ASSIGN EDGE STATE ------------------------------------------------         
        if (row == 1)
            grid.edge_state(cellnum,1) = boundary.top;
        elseif (row == num_cells_y + 2)
            grid.edge_state(cellnum,1) = boundary.bottom;
        elseif (column == 1)
            grid.edge_state(cellnum,1) = boundary.left;
        elseif (column == num_cells_x + 2)
            grid.edge_state(cellnum,1) = boundary.right;
        else
            grid.edge_state(cellnum,1) = 0; %Standard interior cell
        end
        
        
%------- ASSIGN INITIAL CONDITIONS-----------------------------------------
        % if the edge state of current cell is 3 (pressure inlet) or in the
        % first or second column (left-most ghost cells or first column of
        % solution
        if (grid.edge_state(cellnum,1) == 3) || ((ismember(column,[1 2])) && (~ismember(row,[1 num_cells_y+2])))
            grid.initial_conditions(cellnum,1:4) = [perturbed_initial.rho perturbed_initial.rhou perturbed_initial.rhov perturbed_initial.rhoE];
        else
            grid.initial_conditions(cellnum,1:4) = [global_initial.rho global_initial.rhou global_initial.rhov global_initial.rhoE];
        end
        
        % Complete the work on this particular node by incrementing the cell number    
        cellnum = cellnum + 1;        
    end
    
end

% grid.edge_state(grid.edge_state == 3) = 2; %set all "pressure inlet" conditions to "zero gradient" conditions

% Assign adjacent cells to -2 where needed
grid.adjacent_cells(grid.adjacent_cells <= 0) = -2;
grid.adjacent_cells(grid.adjacent_cells > max(grid.cellnumber)) = -2;
        

%% Write this all to a .bkcfd file
fprintf('Writing to .bkcfd file...\n')
struct2csv(grid,[filename '.bkcfd'])

if plot_grid
    figure;
    hold on
    plot(xL,y,x,yB,xR,y,x,yT)
    plot(grid.cornerlocs_x,grid.cornerlocs_y,'.','MarkerSize',10)
end

