#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include "helper_functions.h"

double compute_pressure_2D(TDstate state, double gamma) { 
// This computes the pressure for a 2D state
	double vel_mag_squared = pow((state.rhou/state.rho),2) + pow((state.rhov/state.rho),2);
	double pressure = (gamma-1)*(state.rho*state.rhoE - state.rho*vel_mag_squared/2);

	return(pressure);
}

double compute_pressure_1D(ODstate state, double gamma) {
	double pressure = (gamma-1)*(state.rho*state.rhoE - state.rho*pow(state.rhou,2)/2);
	return(pressure);
}

double vectormax(std::vector<double>& locations, bool absolute_value) {
	double max = -pow(10,50);
	
	for (unsigned int it = 0; it < locations.size(); it++) {
		if (absolute_value) {		
			if (it == 0) {
				max = std::abs(locations[it]);
			} 
			else if (std::abs(locations[it]) > max) {
				max = std::abs(locations[it]);
			}
		} else {
			if (it == 0) {
				max = locations[it];
			} 
			else if (locations[it] > max) {
				max = locations[it];
			}
		}
	}
	return(max);
}

double vectormin(std::vector<double>& locations, bool absolute_value) {
	double min = pow(10,50);

	for (unsigned int it = 0; it < locations.size(); it++) {
		if (absolute_value) {		
			if (it == 0) {
				min = std::abs(locations[it]);
			} 
			else if (std::abs(locations[it]) < min) {
				min = std::abs(locations[it]);
			}
		} else {
			if (it == 0) {
				min = locations[it];
			} 
			else if (locations[it] < min) {
				min = locations[it];
			}
		}
	}
	return(min);
}


std::vector<int> find_interior_cells(const std::vector<cell>& grid) {
	std::vector<int> interior_cells;	
	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) {
		if (grid[cellnum].edge_type == 0) {
			interior_cells.push_back(grid[cellnum].cellnumber);
		}
	}
	return(interior_cells);
}


void timestep_calculator(double& dt, double gamma, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, const std::vector<double>& cell_diameter, const std::vector<double>& cell_perimeter) {

	double a_adjacent, a_center, UV_adjacent, UV_center, wave1_speed, wave2_speed;
	dt = pow(10,10);
	int adjacent_cell_index;

	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell
		if (grid[cellnum].edge_type != 0) {
			continue;
		}
		std::vector<double> edge_length;
		double edge_weighted_average_wavespeed = 0;

		TDstate current_state = Up1[cellnum];
		
		a_center = sqrt(gamma*compute_pressure_2D(current_state, gamma)/current_state.rho);

		for (int edge = 0; edge < 3; ++edge) { // for each of the four edge directions - bottom, right, top, left
			adjacent_cell_index = grid[cellnum].adjacent_cells_gridpos[edge];

			a_adjacent = sqrt(gamma*compute_pressure_2D(Up1[adjacent_cell_index], gamma)/Up1[adjacent_cell_index].rho); // a = sqrt(gamma*P/rho)
			
			if ((edge == 0) || (edge == 2)) {
				UV_adjacent = Up1[adjacent_cell_index].rhov/Up1[adjacent_cell_index].rho;
				UV_center = current_state.rhov/current_state.rho;
			} else {
				UV_adjacent = Up1[adjacent_cell_index].rhou/Up1[adjacent_cell_index].rho;
				UV_center = current_state.rhou/current_state.rho;
			}

// There could be some errors introduced here - not sure exactly how to define left and right states here...
			if ((edge == 3) || (edge == 0)) {
				wave1_speed = UV_adjacent - a_adjacent;
				wave2_speed = UV_center + a_center;
			} else {
				wave1_speed = UV_center - a_center;
				wave2_speed = UV_adjacent + a_adjacent;
			}

			switch (edge) {
				case 0:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(std::abs(wave1_speed),std::abs(wave2_speed))*grid[cellnum].edge_lengths.bottom;
					break;
				case 1:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(std::abs(wave1_speed),std::abs(wave2_speed))*grid[cellnum].edge_lengths.right;
					break;
				case 2:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(std::abs(wave1_speed),std::abs(wave2_speed))*grid[cellnum].edge_lengths.top;
					break;
				case 3:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(std::abs(wave1_speed),std::abs(wave2_speed))*grid[cellnum].edge_lengths.left;
					break;
				default:
					std::cout << "Error in computing timesteps... edge number " << edge << '\n';
					break;
				}
			}
		edge_weighted_average_wavespeed = edge_weighted_average_wavespeed/(cell_perimeter[cellnum]);
		dt = std::min(dt, CFL*cell_diameter[cellnum]/edge_weighted_average_wavespeed);

//		std::cout << "old dt = " << dt << ", new dt = " << CFL*cell_diameter[cellnum]/edge_weighted_average_wavespeed << '\n';
	}
}

void outputU(std::vector<TDstate>& U) {
	int i = -1;
	std::cout << '\n';
	for(auto input:U) {
		i++;
		std::cout << "Cell " << i << ": " << input.rho << " " << input.rhou << " " << input.rhov << " " << input.rhoE << '\n';
	}
	std::cout << '\n';
}

void outputODstate(ODstate& state) {
	std::cout << "ODstate = " << state.rho << " " << state.rhou << " " << state.rhoE << '\n';
}

void outputTDstate(TDstate& state) {
	std::cout << "TDstate = " << state.rho << " " << state.rhou << " " << state.rhov << " " << state.rhoE << '\n';
}

// Write grid to grid_file.txt to see how it is being read in and ensure reading is correct
void output_grid(const std::vector<cell>& grid) {
	std::ofstream ofile;
	ofile.open("grid_file.txt");
	for(auto input: grid) {
		ofile << input.cellnumber << '\n';
		ofile << "cornerlocs_x: " << input.cornerlocs_x[0] << ' ';
		ofile << input.cornerlocs_x[1] << ' ';
		ofile << input.cornerlocs_x[2] << ' ';
		ofile << input.cornerlocs_x[3] << '\n';
		ofile << "cornerlocs_y: " << input.cornerlocs_y[0] << ' ';
		ofile << input.cornerlocs_y[1] << ' ';
		ofile << input.cornerlocs_y[2] << ' ';
		ofile << input.cornerlocs_y[3] << '\n';
		ofile << "adjacent_cells: " << input.adjacent_cells[0] << ' ';
		ofile << input.adjacent_cells[1] << ' ';
		ofile << input.adjacent_cells[2] << ' ';
		ofile << input.adjacent_cells[3] << '\n';
		ofile << "edge_type: " << input.edge_type << '\n';
		ofile << "state: " << input.state.rho << ' ';
		ofile << input.state.rhou << ' ';
		ofile << input.state.rhov << ' ';
		ofile << input.state.rhoE << '\n';
	}
}

// Finds at what vector index number (cellposition) the desired cell number is located
int find_cellposition(const std::vector<cell>& grid, int& cellnumber_desired) {
	for(unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		if (grid[cellposition].cellnumber == cellnumber_desired) {
			return(cellposition);
		}
	}
	return(-1);
}


// Computes the outward unit normal vectors for the given "current_cell" and returns them in two vectors of doubles, [0]-[3] corresponding to left, bottom, right, top
// Should only be used in read_grid.cpp, then values referenced during rest of code
void compute_outward_unit_normal(const cell& current_cell, std::vector<double>& x_unit_normals, std::vector<double>& y_unit_normals) {

	// Compute OUTWARD normal to current cell edge
	for (int edge_number = 0; edge_number < 4; edge_number++) {
		int node_index_1, node_index_2; // node_index_1 and _2 are the two points on this edge (0,1,2,3)
		node_index_1 = edge_number-1;
		node_index_2 = edge_number;
		if (node_index_1 < 0) {
			node_index_1 = 3;
		}

		int direction;
		// Determine the direction - 0 for LR and 1 for BT
		if ((edge_number%2) == 0) {
			direction = 0;
		} else {
			direction = 1;
		}

		double dx = current_cell.cornerlocs_x[node_index_2] - current_cell.cornerlocs_x[node_index_1];
		double dy = current_cell.cornerlocs_y[node_index_2] - current_cell.cornerlocs_y[node_index_1];

		std::vector<double> normal;
		if ((dx > 0) && (direction == 0)) { // Make such that x is always postitive
			x_unit_normals.push_back(-dy/sqrt(pow(dx,2) + pow(dy,2))); 	// x unit normal, see notes from 10/16/15, "Incorporate Rusanov.."
			y_unit_normals.push_back(dx/sqrt(pow(dx,2) + pow(dy,2))); // y unit normal
		} else if ((dx <=  0) && (direction == 0)) {
			x_unit_normals.push_back(dy/sqrt(pow(dx,2) + pow(dy,2))); 	
			y_unit_normals.push_back(-dx/sqrt(pow(dx,2) + pow(dy,2)));
		}

		else if ((dy > 0) && (direction == 1)) {
			x_unit_normals.push_back(dy/sqrt(pow(dx,2) + pow(dy,2))); 	
			y_unit_normals.push_back(-dx/sqrt(pow(dx,2) + pow(dy,2)));
		} else if ((dy <= 0) && (direction == 1)) {
			x_unit_normals.push_back(-dy/sqrt(pow(dx,2) + pow(dy,2))); 	
			y_unit_normals.push_back(dx/sqrt(pow(dx,2) + pow(dy,2)));
		}
	}
}

// Computes area of a cell given the cell's cornerlocs information contained in the cell struct
// Area is computed as magnitude of cross product of diagonals
double compute_cell_area(cell& current_cell) {
	double x_vector1, x_vector2, y_vector1, y_vector2, area;
	
	x_vector1 = current_cell.cornerlocs_x.at(2)-current_cell.cornerlocs_x.at(0);
	x_vector2 = current_cell.cornerlocs_x.at(3)-current_cell.cornerlocs_x.at(1);
	y_vector1 = current_cell.cornerlocs_y.at(2)-current_cell.cornerlocs_y.at(0);
	y_vector2 = current_cell.cornerlocs_y.at(3)-current_cell.cornerlocs_y.at(1);

	area = std::abs(0.5*(x_vector1*y_vector2 - x_vector2*y_vector1));

// Alternate way to calculate: compute cell areas using Gauss's Area Formula
//		area = 0.5*std::abs(grid[cellnum].cornerlocs_x[0]*grid[cellnum].cornerlocs_y[1] + grid[cellnum].cornerlocs_x[1]*grid[cellnum].cornerlocs_y[2] + grid[cellnum].cornerlocs_x[2]*grid[cellnum].cornerlocs_y[3] + grid[cellnum].cornerlocs_x[3]*grid[cellnum].cornerlocs_y[0] - grid[cellnum].cornerlocs_x[1]*grid[cellnum].cornerlocs_y[0] - grid[cellnum].cornerlocs_x[2]*grid[cellnum].cornerlocs_y[1] - grid[cellnum].cornerlocs_x[3]*grid[cellnum].cornerlocs_y[2] - grid[cellnum].cornerlocs_x[0]*grid[cellnum].cornerlocs_y[3]);
	return(area);
}


void compute_cell_distances(std::vector<cell>& grid, directional_quantity& cell_distance, unsigned int& cellposition) {
	directional_quantity delta;
	directional_quantity delta_x, delta_y;

	delta_x.left 	= (grid[grid[cellposition].adjacent_cells_gridpos[0]].centroid_x - grid[cellposition].centroid_x);
	delta_y.left 	= (grid[grid[cellposition].adjacent_cells_gridpos[0]].centroid_y - grid[cellposition].centroid_y);
	delta_x.bottom = (grid[grid[cellposition].adjacent_cells_gridpos[1]].centroid_x - grid[cellposition].centroid_x);
	delta_y.bottom = (grid[grid[cellposition].adjacent_cells_gridpos[1]].centroid_y - grid[cellposition].centroid_y);
	delta_x.right 	= (grid[grid[cellposition].adjacent_cells_gridpos[2]].centroid_x - grid[cellposition].centroid_x);
	delta_y.right 	= (grid[grid[cellposition].adjacent_cells_gridpos[2]].centroid_y - grid[cellposition].centroid_y);
	delta_x.top 	= (grid[grid[cellposition].adjacent_cells_gridpos[3]].centroid_x - grid[cellposition].centroid_x);
	delta_y.top 	= (grid[grid[cellposition].adjacent_cells_gridpos[3]].centroid_y - grid[cellposition].centroid_y);

	cell_distance.left 	= sqrt(pow(delta_x.left,2) 	+ pow(delta_y.left,2));
	cell_distance.right 	= sqrt(pow(delta_x.right,2) 	+ pow(delta_y.right,2));
	cell_distance.bottom = sqrt(pow(delta_x.bottom,2) 	+ pow(delta_y.bottom,2));
	cell_distance.top 	= sqrt(pow(delta_x.top,2) 		+ pow(delta_y.top,2));

}

void compute_cell_edge_length(cell& current_cell, directional_quantity& cell_edge_lengths) {
	
	cell_edge_lengths.bottom	= sqrt(pow((current_cell.cornerlocs_x[1]-current_cell.cornerlocs_x[0]),2) + pow((current_cell.cornerlocs_y[1]-current_cell.cornerlocs_y[0]),2));
	cell_edge_lengths.right 	= sqrt(pow((current_cell.cornerlocs_x[2]-current_cell.cornerlocs_x[1]),2) + pow((current_cell.cornerlocs_y[2]-current_cell.cornerlocs_y[1]),2));
	cell_edge_lengths.top 		= sqrt(pow((current_cell.cornerlocs_x[3]-current_cell.cornerlocs_x[2]),2) + pow((current_cell.cornerlocs_y[3]-current_cell.cornerlocs_y[2]),2));
	cell_edge_lengths.left 		= sqrt(pow((current_cell.cornerlocs_x[0]-current_cell.cornerlocs_x[3]),2) + pow((current_cell.cornerlocs_y[0]-current_cell.cornerlocs_y[3]),2));

}
