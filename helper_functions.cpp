#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include "helper_functions.h"

double compute_pressure_2D(TDstate state, double gamma) { 
// This computes the pressure for a 2D state
	double vel_mag_squared = pow((state.rhou/state.rho),2) + pow((state.rhov/state.rho),2);
	double pressure = (gamma-1)*(state.rhoE - state.rho*vel_mag_squared/2);


	return(pressure);
}

double compute_pressure_1D(ODstate state, double gamma) {
	double pressure = (gamma-1)*(state.rhoE - state.rho*pow(state.rhou,2)/2);
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

double timestep_calculator(double gamma, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, const std::vector<double>& cell_diameter, const std::vector<double>& cell_perimeter) {

	double a_adjacent, a_center, U_adjacent, V_adjacent, U_center, V_center, velmag_center, velmag_adjacent, wavespeed_center, wavespeed_adjacent;
	double dt = pow(10,10);
	int adjacent_cell_index;

	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell
		if (grid[cellnum].edge_type != 0) { // Don't take ghost cells into account for dt calculation
			continue;
		}
		std::vector<double> edge_length;
		double edge_weighted_average_wavespeed = 0;

		TDstate current_state = Up1[cellnum];
		
		a_center = sqrt(gamma*compute_pressure_2D(current_state, gamma)/current_state.rho);
		U_center = current_state.rhou/current_state.rho;
		V_center = current_state.rhov/current_state.rho;

		for (int edge = 0; edge < 3; ++edge) { // for each of the four edge directions - left, bottom, right, top 
			adjacent_cell_index = grid[cellnum].adjacent_cells_gridpos[edge];

			a_adjacent = sqrt(gamma*compute_pressure_2D(Up1[adjacent_cell_index], gamma)/Up1[adjacent_cell_index].rho); // a = sqrt(gamma*P/rho)
/*		The Old Way	
			if ((edge == 0) || (edge == 2)) {
				UV_adjacent = Up1[adjacent_cell_index].rhou/Up1[adjacent_cell_index].rho;
				UV_center = current_state.rhou/current_state.rho;
			} else {
				UV_adjacent = Up1[adjacent_cell_index].rhov/Up1[adjacent_cell_index].rho;
				UV_center = current_state.rhov/current_state.rho;
			}

// There could be some errors introduced here - not sure exactly how to define left and right states here...
			if ((edge == 3) || (edge == 0)) {
				wave1_speed = UV_adjacent - a_adjacent;
				wave2_speed = UV_center + a_center;
			} else {
				wave1_speed = UV_center - a_center;
				wave2_speed = UV_adjacent + a_adjacent;
			}
*/

			// The New Way

			// Compute U and V in adjacent cell
			U_adjacent = Up1[adjacent_cell_index].rhou/Up1[adjacent_cell_index].rho;
			V_adjacent = Up1[adjacent_cell_index].rhov/Up1[adjacent_cell_index].rho;

			// Compute normal speed in adjacent cell and center cell - normal to current cell edge
			velmag_adjacent = U_adjacent*grid[cellnum].outward_unit_normals_x[edge] + V_adjacent*grid[cellnum].outward_unit_normals_y[edge];
			velmag_center = U_center*grid[cellnum].outward_unit_normals_x[edge] + V_center*grid[cellnum].outward_unit_normals_y[edge];

			// Determine adjacent max wavespeed and center max wavespeed
			wavespeed_adjacent = std::max(std::abs(velmag_adjacent + a_adjacent),std::abs(velmag_adjacent - a_adjacent));
			wavespeed_center = std::max(std::abs(velmag_center + a_center),std::abs(velmag_center - a_center));
			
			switch (edge) {
				case 0:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(wavespeed_center, wavespeed_adjacent) *grid[cellnum].edge_lengths.left;
					break;
				case 1:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(wavespeed_center, wavespeed_adjacent) *grid[cellnum].edge_lengths.bottom;
					break;
				case 2:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(wavespeed_center, wavespeed_adjacent) *grid[cellnum].edge_lengths.right;
					break;
				case 3:
					edge_weighted_average_wavespeed = edge_weighted_average_wavespeed + std::max(wavespeed_center, wavespeed_adjacent) *grid[cellnum].edge_lengths.top;
					break;
				default:
					std::cout << "Error in computing timesteps... edge number " << edge << '\n';
					break;
				}
			}
		edge_weighted_average_wavespeed = edge_weighted_average_wavespeed/(cell_perimeter[cellnum]);
	//	std::cout << "old dt = " << dt << ", new dt = " << CFL*cell_diameter[cellnum]/edge_weighted_average_wavespeed << '\n';		
		dt = std::min(dt, CFL*cell_diameter[cellnum]/edge_weighted_average_wavespeed);
		//std::cout << "cell number: " << grid[cellnum].cellnumber << " wavespeed: " << edge_weighted_average_wavespeed << ", this C.Dia: " << cell_diameter[cellnum] << ", this dt: " << CFL*cell_diameter[cellnum]/edge_weighted_average_wavespeed <<'\n';

	}
	return(dt);
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
	std::cout << state.rho << " " << state.rhou << " " << state.rhov << " " << state.rhoE << '\n';
}

// Write grid to grid_file.txt to see how it is being read in and ensure reading is correct
void output_grid(const std::vector<cell>& grid) {
	std::ofstream ofile;
	ofile.open("grid_file.txt");
	for(auto input: grid) {
		ofile << "cell number " << input.cellnumber << '\n';
		ofile << "cornerlocs_x: " << input.cornerlocs_x[0] << ' ' << input.cornerlocs_x[1] << ' ' << input.cornerlocs_x[2] << ' ' << input.cornerlocs_x[3] << "  size: " << input.cornerlocs_x.size() << '\n';
		ofile << "cornerlocs_y: " << input.cornerlocs_y[0] << ' ' << input.cornerlocs_y[1] << ' ' << input.cornerlocs_y[2] << ' ' << input.cornerlocs_y[3] << "  size: " << input.cornerlocs_y.size() << '\n';
		ofile << "adjacent_cells: " << input.adjacent_cells[0] << ' ' << input.adjacent_cells[1] << ' ' << input.adjacent_cells[2] << ' ' << input.adjacent_cells[3] << "  size: " << input.adjacent_cells.size() << '\n';
		ofile << "edge_type: " << input.edge_type << '\n';
		ofile << "state: " << input.state.rho << ' ' << input.state.rhou << ' ' << input.state.rhov << ' ' << input.state.rhoE << '\n';
		ofile << "x unit normals: " << input.outward_unit_normals_x[0] << ' ' << input.outward_unit_normals_x[1] << ' ' << input.outward_unit_normals_x[2] << ' ' << input.outward_unit_normals_x[3] << "  size: " << input.outward_unit_normals_x.size() << '\n';
		ofile << "y unit normals: " << input.outward_unit_normals_y[0] << ' ' << input.outward_unit_normals_y[1] << ' ' << input.outward_unit_normals_y[2] << ' ' << input.outward_unit_normals_y[3] << "  size: " << input.outward_unit_normals_y.size() << '\n';
		ofile << "edge lengths: " << input.edge_lengths.left << ' ' << input.edge_lengths.bottom << ' ' << input.edge_lengths.right << ' ' << input.edge_lengths.top << "\n\n";
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
	// std::cout << "Computing normals for cell................. " << current_cell.cellnumber << '\n';
	// Compute OUTWARD normal to current cell edge
	for (int edge_number = 0; edge_number < 4; edge_number++) {
		int node_index_1, node_index_2; // node_index_1 and _2 are the two points on this edge (0,1,2,3)
		node_index_1 = edge_number-1;
		node_index_2 = edge_number;
		if (node_index_1 < 0) {
			node_index_1 = 3;
		}
/*
		int direction;
		// Determine the direction - 0 for LR and 1 for BT
		if ((edge_number%2) == 0) {
			direction = 0;
		} else {
			direction = 1;
		}
*/
		double dx = current_cell.cornerlocs_x[node_index_2] - current_cell.cornerlocs_x[node_index_1];
		double dy = current_cell.cornerlocs_y[node_index_2] - current_cell.cornerlocs_y[node_index_1];
		double cell_edge_length = sqrt(dx*dx + dy*dy);

		x_unit_normals[edge_number] = (dy/cell_edge_length); 	// x unit normal, see notes from 10/16/15, "Incorporate Rusanov.."
		y_unit_normals[edge_number] = (-dx/cell_edge_length); // y unit normal

		if (std::abs(x_unit_normals[edge_number]) < pow(10,-20)) {
			x_unit_normals[edge_number] = 0;
		}
		if (std::abs(y_unit_normals[edge_number]) < pow(10,-20)) {
			y_unit_normals[edge_number] = 0;
		}

	//	double nx = dy/cell_edge_length;
	//	double ny = -dx/cell_edge_length;

/*
		debugging for flux errors 
		std::cout << "cornerlocs_x: " << current_cell.cornerlocs_x[node_index_2] << " - " << current_cell.cornerlocs_x[node_index_1] << '\n';
		std::cout << "cornerlocs_y: " << current_cell.cornerlocs_y[node_index_2] << " - " << current_cell.cornerlocs_y[node_index_1] << '\n';
		std::cout << "dx = " << dx << " dy = " << dy << '\n';		
		
		std::cout << "Check: should equal dx and dy: " << -1.0*y_unit_normals.at(edge_number)*cell_edge_length << ", " << x_unit_normals.at(edge_number)*cell_edge_length << '\n';
	
		std::cout << "edge length: " << cell_edge_length << '\n';
		std::cout << "x: " << x_unit_normals[edge_number] << ", y: " << y_unit_normals[edge_number] << '\n';
*/
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


// Compute distance between cell centroids. Stored/returned in cell_distance, a directional_quantity
void compute_cell_distances(std::vector<cell>& grid, directional_quantity& cell_distance, unsigned int& cellposition) {

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

void compute_edge_normal_flux(TDstate& F, TDstate& G, double x_normal, double y_normal, TDstate& edge_directional_flux) {
	edge_directional_flux.rho  = F.rho *x_normal + G.rho *y_normal;
	edge_directional_flux.rhou = F.rhou*x_normal + G.rhou*y_normal;
	edge_directional_flux.rhov = F.rhov*x_normal + G.rhov*y_normal;
	edge_directional_flux.rhoE = F.rhoE*x_normal + G.rhoE*y_normal;	
}

void weight_flux_with_area(TDstate& flux_to_weight, double edge_length) {
	flux_to_weight.rho  = flux_to_weight.rho *edge_length;
	flux_to_weight.rhou = flux_to_weight.rhou*edge_length;
	flux_to_weight.rhov = flux_to_weight.rhov*edge_length;
	flux_to_weight.rhoE = flux_to_weight.rhoE*edge_length;
}

void output_residual_sum(std::vector<TDstate>& residuals) { //Note that cellnum == 0 is not used, residuals start on cell 1 since there is no cell 0!
	double rho_residual_sum = 0, rhou_residual_sum = 0, rhov_residual_sum = 0, rhoE_residual_sum = 0;
	for (unsigned int cellnum = 0; cellnum < residuals.size(); ++cellnum) {
		rho_residual_sum  = rho_residual_sum  + residuals[cellnum].rho;
		rhou_residual_sum = rhou_residual_sum + residuals[cellnum].rhou;
		rhov_residual_sum = rhov_residual_sum + residuals[cellnum].rhov;
		rhoE_residual_sum = rhoE_residual_sum + residuals[cellnum].rhoE;	
		//std::cout << cellnum << " " << residuals[cellnum].rhov << " " << rhov_residual_sum << '\n';
	}
	std::cout << "  Residuals--rho: " << rho_residual_sum << " rhou: " << rhou_residual_sum << " rhov: " << rhov_residual_sum << " rhoE: " << rhoE_residual_sum << '\n';
}

