#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "file_header.h"
#include "helper_functions.h"

/*
CODE NOTES

- The inlet should be specified as a zero-gradient boundary condition. However, when this happens, the code explodes, with density becoming zero/negative by about timestep 74 (in the case I ran).




*/


// Define helper functions!

// Compute thr first instance of dt, changing values of dt, cell_perimeter and cell_diameter
void initialize_dt(double dt, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, double gamma, std::vector<directional_quantity> & cell_edge_lengths, std::vector<double>& cell_perimeter, std::vector<double>& cell_diameter);

// 
void Reconstruction_Flux(const std::vector<cell>& grid, directional_quantity& delta, TDstate& state_minus, TDstate& state_plus, double& gamma, TDstate& slope_plus, TDstate& slope_minus, int& limiter_number, unsigned int& cellposition, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double& CFL, std::vector<TDstate>& U, bool& debug, std::ofstream& textout);

// Create the correct input and solve the exact riemann problem
void Solve_Riemann_Flux( int& flux_function_choice, unsigned int& cellposition, double& gamma, double& thresh, bool& debug, std::vector<TDstate>& Uph, cell& current_cell, std::vector<TDstate>& F, std::vector<TDstate>& G );

// Update ghost cells to reflect for wall boundaries
void update_ghost_cells(std::vector<cell>& grid, std::vector<TDstate>& U, std::vector<TDstate>& Up1, std::vector<int>& interior_cells, bool debug, std::ofstream& textout);

// Compute the Rusanov flux
TDstate Rusanov(TDstate& state_minus, TDstate& state_plus, std::vector<double>& x_unit_normals, std::vector<double>& y_unit_normals, double& gamma, int edge_number);

// Return bool if cell is an edge cell or not
bool isedge(const std::vector<cell>& grid, int cellposition);

// void compute_flux(std::vector<double> &fluxph, std::vector<double> &fluxmh, std::vector<double> slope, TDstate state, double delta, std::string direction, double gamma);

// This is not math-y, just a formality of grabbing the left/right (state_minus) or bottom/top (state_plus) states from the grid input.
void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<TDstate>& U, const cell& current_cell, int direction);

// Computes the "geometric" slope vector for all conserved quantities between the bottom/left and top/right cells. Used for second order accuracy calculations.
void compute_slope(TDstate &slope_minus, TDstate &slope_plus, int direction, TDstate state_minus, TDstate state_plus, directional_quantity delta, TDstate& U_center);

// This computes the double-value "limiter value" based on the bottom/left and top/right states for a given cell. The limiter value is then used to calculate the flux
void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number);

// Computes the harmonic limiter given r, the ratio of minus to plus slopes
void harmonic_limiter(TDstate &limiter_value, TDstate r);

// Computes the van leers limiter 
void van_leers_limiter(TDstate &limiter_value, TDstate r);

// Computes the Superbee limiter
void superbee_limiter(TDstate &limiter_value, TDstate r);

// Computes the limited flux for two cells (left-right or bottom-top) of the reconstruction step
void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, int& direction, double& gamma, double& CFL, bool& debug, std::ofstream& textout);

// Computes the 2D flux given a TDstate "state" on the edge
void compute_cell_flux(TDstate& flux, const TDstate& state, int direction, double& gamma);

// Reconstructs the solution from t to t+dt/2
void compute_halfway_state(TDstate& Uph, const TDstate& current_U, const std::vector<TDstate>& F, const std::vector<TDstate>& G, double& dt, directional_quantity& delta);

// Computes the min dx or dy (linear) in a grid
double gridmin(const std::vector<cell>& grid, char direction);

// Outputs the grid to the file "grid_file.txt"
//void output_grid(const std::vector<cell>& grid);




int main() {

// Define Variables
	std::string parameter_filename = "solver_parms.txt";

	std::vector<std::string> parameters;
	std::vector<double> Fiph, Fimh, Giph, Gimh;
	std::vector<double> cell_perimeter, cell_diameter;
	directional_quantity delta;
	TDstate state_minus, state_plus, slope_minus, slope_plus, limiter_value;
	std::vector<TDstate> F (2);
	std::vector<TDstate> G (2);

	double thresh = pow(10,-6); // 1E-6 tolerance for exact riemann solver - velocity
	double gamma = 1.4;

	bool save_timesteps = true;	

	std::vector<cell> grid;

	// Read parameters from input file defined by filename string
	parameters = read_parameters(parameter_filename);
	double CFL = std::stod(parameters.at(0));	
	int num_timesteps = std::stod(parameters.at(1));
	std::string input_filename = parameters.at(2);
	std::string output_filename = parameters.at(3);
	
	bool debug;	
	if ((parameters.at(4).compare("True")) == 0) {
		debug = 1;
	} else {
		debug = 0;
	}
	std::cout << "Debug Mode: " << debug << "\n";
	std::ofstream textout;
	if (debug) {
		textout.open("debug_file.txt");
	}
	int flux_function_choice = std::stoi(parameters.at(5));
	int limiter_number = std::stoi(parameters.at(6));
	


// Read initial file defined from Matlab with cell edges and connections
// Read the initial condition file from MATLAB with [rho rho*u rho*v E] defined for each cell
	read_grid(input_filename, grid, gamma);
	
// Generate a vector of ints containing the cell numbers that are interior cells (0, not 1 or 2)
	std::vector<int> interior_cells = find_interior_cells(grid);

// Option to output grid to "grid_file.txt"
//	output_grid(grid);

// Initialize solution U as vector of TD states with initial state at all points on grid
	std::vector<TDstate> U;
	for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		U.push_back(grid[cellposition].state);
	}

// Before starting calculation, make sure ghost cells have been correctly initialized to correct normal condition
	update_ghost_cells(grid, U, U, interior_cells, debug, textout);

// Initialize needed variables
	TDstate Uph;
	std::vector<TDstate> Up1 (grid.size());
	std::vector<TDstate> Uph_vector (grid.size());
	std::vector<directional_quantity> cell_edge_lengths;

// Initialize the dt value
	double dt, t = 0;
	initialize_dt(dt, grid, Up1, CFL, gamma, cell_edge_lengths, cell_perimeter, cell_diameter);

	if (debug) {
		textout << "\nInitial wavespeeds and dt computed, moving on to flux calculation \n \n";
	}

// Start calculation in for loop going from t = 0 to final time
	for (int timestep = 0; timestep < num_timesteps; timestep++) { // for each timestep
		std::cout << "Timestep " << timestep+1 << " of " << num_timesteps << ", dt = " << dt << '\n';
		if (timestep > 0) {
			U = Up1;
		}

	// Go through each cell, calculate fluxes, update to new piecewise linear state, find flux and update to state Uph
		for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) { // For each cell
			if (grid[cellposition].edge_type == 0)  { // Cell is not a ghost cell/edge cell - only do calculations on interior cells

			// Compute F and G for the interpolation step - all else is junk needed for calculation of F and G
   	   	Reconstruction_Flux(grid, delta, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellposition, limiter_value, F, G, CFL, U, debug, textout);

			// Now, use the fluxes calculated above to assign a new TDstate for the current cell, Uph
			//	std::cout << "Flux has been constructed for " << cellposition << ", Time to compute half state \n";
				compute_halfway_state(Uph, U[cellposition], F, G, dt, delta);

				Uph_vector[cellposition] = Uph;
				
				if (debug) {
					textout << "Cell " << grid[cellposition].cellnumber << " Original State: " <<  U[cellposition].rho << " " << U[cellposition].rhou << " " << U[cellposition].rhov << " " << U[cellposition].rhoE << '\n';
					textout << "Reconstruction: \nF[0] = " << F[0].rho << " " << F[0].rhou << " " << F[0].rhov << " " << F[0].rhoE << '\n';
					textout << "F[1] = " << F[1].rho << " " << F[1].rhou << " " << F[1].rhov << " " << F[1].rhoE << '\n';
					textout << "G[0] = " << G[0].rho << " " << G[0].rhou << " " << G[0].rhov << " " << G[0].rhoE << '\n';
					textout << "G[1] = " << G[1].rho << " " << G[1].rhou << " " << G[1].rhov << " " << G[1].rhoE << '\n';
					textout << "Cell " << grid[cellposition].cellnumber << " half-state: " <<  Uph.rho << " " << Uph.rhou << " " << Uph.rhov << " " << Uph.rhoE << '\n';
				}
			} else {
				Uph_vector[cellposition] = U[cellposition];
			}
		}

	// Go through each cell, solve riemann problem (exact or not exact), calculate fluxes and update to state Up1
		for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) { // For each cell
			if (grid[cellposition].edge_type == 0)  { // Cell is not a ghost cell/edge cell - only do calculations on interior cells

			// Solve riemann problem on edge with specified flux function for euler flux
				Solve_Riemann_Flux(flux_function_choice, cellposition, gamma, thresh, debug, Uph_vector, grid[cellposition], F, G);

			// Use all these fluxes to define updated state on each cell 
				compute_halfway_state(Up1[cellposition], Uph_vector[cellposition], F, G, dt, delta); 
		
				if (debug) {
					textout << "Halfway State = " << Uph_vector[cellposition].rho << " " << Uph_vector[cellposition].rhou << " " << Uph_vector[cellposition].rhov << " " << Uph_vector[cellposition].rhoE << '\n';
					textout << "Riemann Fluxes: \nF[0] = " << F[0].rho << " " << F[0].rhou << " " << F[0].rhov << " " << F[0].rhoE << '\n';
					textout << "F[1] = " << F[1].rho << " " << F[1].rhou << " " << F[1].rhov << " " << F[1].rhoE << '\n';
					textout << "G[0] = " << G[0].rho << " " << G[0].rhou << " " << G[0].rhov << " " << G[0].rhoE << '\n';
					textout << "G[1] = " << G[1].rho << " " << G[1].rhou << " " << G[1].rhov << " " << G[1].rhoE << '\n';
					textout << "Cell " << grid[cellposition].cellnumber << " Full-Updated State : " <<  Up1[cellposition].rho << " " << Up1[cellposition].rhou << " " << Up1[cellposition].rhov << " " << Up1[cellposition].rhoE << '\n' << '\n';
				}
			} // end of if (grid[cellposition].edge_type == 0) - loop through all interior cells
		} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) - loop through all timesteps


	// Now that the updated solution has been calculated, we need to redefine the ghost cells to the correct values
		update_ghost_cells(grid, U, Up1, interior_cells, debug, textout);

		//if ((save_timesteps) && ((timestep%10 == 0) || timestep == num_timesteps)) {
		if (save_timesteps) {
			write_to_file(grid,Up1,output_filename);
		}
		
		//Update dt based on new wavespeeds
		timestep_calculator(dt, gamma, grid, Up1, CFL, cell_edge_lengths, cell_diameter, cell_perimeter);
		
		textout << '\n';

		t = t + dt;

	} // end of for(int t = 0; t <= 2*dt; t++) 
	textout.close();
	return 0;
}



// This is the land of the helper functions! ---------------------------------------------------

void initialize_dt(double dt, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, double gamma, std::vector<directional_quantity>& cell_edge_lengths, std::vector<double>& cell_perimeter, std::vector<double>& cell_diameter) {

	double area;

	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell		

		cell_edge_lengths.push_back(directional_quantity());

		cell_edge_lengths[cellnum].bottom 	= sqrt(pow((grid[cellnum].cornerlocs_x[1]-grid[cellnum].cornerlocs_x[0]),2) + pow((grid[cellnum].cornerlocs_y[1]-grid[cellnum].cornerlocs_y[0]),2));
		cell_edge_lengths[cellnum].right 	= sqrt(pow((grid[cellnum].cornerlocs_x[2]-grid[cellnum].cornerlocs_x[1]),2) + pow((grid[cellnum].cornerlocs_y[2]-grid[cellnum].cornerlocs_y[1]),2));
		cell_edge_lengths[cellnum].top 		= sqrt(pow((grid[cellnum].cornerlocs_x[3]-grid[cellnum].cornerlocs_x[2]),2) + pow((grid[cellnum].cornerlocs_y[3]-grid[cellnum].cornerlocs_y[2]),2));
		cell_edge_lengths[cellnum].left 		= sqrt(pow((grid[cellnum].cornerlocs_x[0]-grid[cellnum].cornerlocs_x[3]),2) + pow((grid[cellnum].cornerlocs_y[0]-grid[cellnum].cornerlocs_y[3]),2));

		cell_perimeter.push_back(cell_edge_lengths[cellnum].bottom + cell_edge_lengths[cellnum].right + cell_edge_lengths[cellnum].top + cell_edge_lengths[cellnum].left);

		// Compute cell areas using Gauss's Area Formula
		area = 0.5*std::abs(grid[cellnum].cornerlocs_x[0]*grid[cellnum].cornerlocs_y[1] + grid[cellnum].cornerlocs_x[1]*grid[cellnum].cornerlocs_y[2] + grid[cellnum].cornerlocs_x[2]*grid[cellnum].cornerlocs_y[3] + grid[cellnum].cornerlocs_x[3]*grid[cellnum].cornerlocs_y[0] - grid[cellnum].cornerlocs_x[1]*grid[cellnum].cornerlocs_y[0] - grid[cellnum].cornerlocs_x[2]*grid[cellnum].cornerlocs_y[1] - grid[cellnum].cornerlocs_x[3]*grid[cellnum].cornerlocs_y[2] - grid[cellnum].cornerlocs_x[0]*grid[cellnum].cornerlocs_y[3]);

		cell_diameter.push_back(2*area/cell_perimeter[cellnum]);
	}
	
	// Now that we have cell_diameter and perimeter, compute timestep	
	timestep_calculator(dt, gamma, grid, Up1, CFL, cell_edge_lengths, cell_diameter, cell_perimeter);
}


// Reconstruction_Flux computes the flux for the initial t -> t+1/2 update
// Everything is passed by reference, only items that matter for returning are F and G, which are two-element vectors of TDstates
// This is computed prior to riemann problem for each timestep update
void Reconstruction_Flux(const std::vector<cell>& grid, directional_quantity& delta, TDstate& state_minus, TDstate& state_plus, double& gamma, TDstate& slope_plus, TDstate& slope_minus, int& limiter_number, unsigned int& cellposition, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double& CFL, std::vector<TDstate>& U, bool& debug, std::ofstream& textout) {

	directional_quantity delta_x, delta_y;	
	std::vector<TDstate> temp_F, temp_G;

	delta_x.left 	= (grid[grid[cellposition].adjacent_cells_gridpos[0]].centroid_x - grid[cellposition].centroid_x);
	delta_y.left 	= (grid[grid[cellposition].adjacent_cells_gridpos[0]].centroid_y - grid[cellposition].centroid_y);
	delta_x.bottom = (grid[grid[cellposition].adjacent_cells_gridpos[1]].centroid_x - grid[cellposition].centroid_x);
	delta_y.bottom = (grid[grid[cellposition].adjacent_cells_gridpos[1]].centroid_y - grid[cellposition].centroid_y);
	delta_x.right 	= (grid[grid[cellposition].adjacent_cells_gridpos[2]].centroid_x - grid[cellposition].centroid_x);
	delta_y.right 	= (grid[grid[cellposition].adjacent_cells_gridpos[2]].centroid_y - grid[cellposition].centroid_y);
	delta_x.top 	= (grid[grid[cellposition].adjacent_cells_gridpos[3]].centroid_x - grid[cellposition].centroid_x);
	delta_y.top 	= (grid[grid[cellposition].adjacent_cells_gridpos[3]].centroid_y - grid[cellposition].centroid_y);

	delta.left 		= sqrt(pow(delta_x.left,2) 	+ pow(delta_y.left,2));
	delta.right 	= sqrt(pow(delta_x.right,2) 	+ pow(delta_y.right,2));
	delta.bottom 	= sqrt(pow(delta_x.bottom,2) 	+ pow(delta_y.bottom,2));
	delta.top 		= sqrt(pow(delta_x.top,2) 		+ pow(delta_y.top,2));
	
   // Calculate the limited flux of each conserved quantity for initial t -> t+1/2 update
	// This could be changed to go for each of the 4 directions... but is it needed? May make code simpler to understands
	for (int direction = 0; direction < 2; ++direction) { // For the left/right direction and bottom/top direction, two loop iterations only
		
		// Grabs the bottom/left and the top/right states for current cell, depending on direction
		// state_minus and state_plus are the bottom/left and top/right TDstates, respectively
		// U is the vector of TDstates that is the current U for all states, not the initial U
		assign_edge_state(state_minus, state_plus, U, grid[cellposition], direction);	
		
		// Calculate the flux of each conserved variable, in given direction (0 L-R for first iteration, 1 B-T for second iteration), of that cell.
		// Flux needs to be calculated by first finding the limiter (direction-dependent), then using that limiter to calculate the directional flux, then using the two directional fluxes to update from U(t) to U(t+1/2)
			
		compute_slope(slope_minus, slope_plus, direction, state_minus, state_plus, delta, U[cellposition]);
	
		if ((direction == 0) && (debug)) {
			textout << "Slope LR Minus: " << slope_minus.rho << " " << slope_minus.rhou << " " << slope_minus.rhov << " " << slope_minus.rhoE << '\n';
			textout << "Slope LR Plus: " << slope_plus.rho << " " << slope_plus.rhou << " " << slope_plus.rhov << " " << slope_plus.rhoE << '\n';
		}

//limiter_value is the B term on book pg. 238
		compute_limiter(limiter_value, slope_minus, slope_plus, limiter_number);

		if (debug) {
			textout << "Limiter Num: " << limiter_number << " Limiter Value: " << limiter_value.rho << " " << limiter_value.rhou << " " << limiter_value.rhov << " " << limiter_value.rhoE << '\n';
		}

		if (direction == 0) { // Compute F, x-fluxes
			compute_limited_flux(F, limiter_value, U[cellposition], direction, gamma, CFL, debug, textout);	
		} else { // Compute G, y-fluxes
			compute_limited_flux(G, limiter_value, U[cellposition], direction, gamma, CFL, debug, textout);
		}

/*		F.rho  = F.rho  + temp_F.at(1.rho;
		F.rhou = F.rhou + temp_F.rhou;
		F.rhov = F.rhov + temp_F.rhov;
		F.rhoE = F.rhoE + temp_F.rhoE;

		G.rho  = G.rho  + temp_G.rho;
		G.rhou = G.rhou + temp_G.rhou;
		G.rhov = G.rhov + temp_G.rhov;
		G.rhoE = G.rhoE + temp_G.rhoE;*/

	} // end of for (int direction = 0; direction < 2; ++direction)
}


void Solve_Riemann_Flux( int& flux_function_choice, unsigned int& cellposition, double& gamma, double& thresh, bool& debug, std::vector<TDstate>& Uph, cell& current_cell, std::vector<TDstate>& F, std::vector<TDstate>& G ) {

	std::vector<double> x_unit_normals = current_cell.unit_normals_x;
	std::vector<double> y_unit_normals = current_cell.unit_normals_y;	

	TDstate state_center_2D = Uph[cellposition];
	TDstate tempflux, state_leftcell_2D, state_rightcell_2D, state_bottomcell_2D, state_topcell_2D;
	ODstate state_leftcell_1D, center_lr, center_bt, state_rightcell_1D, state_bottomcell_1D, state_topcell_1D;
				
	// Assign what states are located on the left and right/bottom and top of the current cell. Easier to grab them here and store in a variable changed each iteration than later
	assign_edge_state(state_leftcell_2D, state_rightcell_2D, Uph, current_cell, 0); //direction 0 == left-right
	assign_edge_state(state_bottomcell_2D, state_topcell_2D, Uph, current_cell, 1); //direction 1 == bottom-top
	
	if (flux_function_choice == 1) { //Roe Flux

	} else if (flux_function_choice == 2) { // Rusanov Flux - this takes cell normals into account! Yey!	
		F[0] = Rusanov(state_leftcell_2D,   state_center_2D,    x_unit_normals, y_unit_normals, gamma, 0); 	// Edge 0 - left
		F[1] = Rusanov(state_center_2D,     state_rightcell_2D, x_unit_normals, y_unit_normals, gamma, 2); 	// Edge 2 - right 
		G[0] = Rusanov(state_bottomcell_2D, state_center_2D,    x_unit_normals, y_unit_normals, gamma, 1);		// Edge 1 - bottom
		G[1] = Rusanov(state_center_2D,     state_topcell_2D,   x_unit_normals, y_unit_normals, gamma, 3); 	// Edge 3 - top
	} else if (flux_function_choice == 3) { // Exact Riemann Problem

// All this mess is because the riemann problem is 1D and we store 2D states.
// So we have to make new, temporary 1D states and fixes (all the doubles) to make the Riemann Problem work.
		center_lr.rho = state_center_2D.rho;
		center_lr.rhou = state_center_2D.rhou;
		double center_lr_parallelvel = state_center_2D.rhov/state_center_2D.rho; 
		center_lr.rhoE = state_center_2D.rhoE;

		center_bt.rho = state_center_2D.rho;
		double center_bt_parallelvel = state_center_2D.rhou/state_center_2D.rho; 
		center_bt.rhou = state_center_2D.rhov;	
		center_bt.rhoE = state_center_2D.rhoE;
					
		state_leftcell_1D.rho = state_leftcell_2D.rho;
		state_leftcell_1D.rhou = state_leftcell_2D.rhou;
		double leftcell_parallelvel = state_leftcell_2D.rhov/state_leftcell_2D.rho;
		state_leftcell_1D.rhoE = state_leftcell_2D.rhoE;				

		state_rightcell_1D.rho = state_rightcell_2D.rho; 
		state_rightcell_1D.rhou = state_rightcell_2D.rhou;
		double rightcell_parallelvel = state_rightcell_2D.rhov/state_rightcell_2D.rho; 
		state_rightcell_1D.rhoE = state_rightcell_2D.rhoE;	

		state_bottomcell_1D.rho = state_bottomcell_2D.rho;	
		double bottomcell_parallelvel = state_bottomcell_2D.rhou/state_bottomcell_2D.rho;
		state_bottomcell_1D.rhou = state_bottomcell_2D.rhov;
		state_bottomcell_1D.rhoE = state_bottomcell_2D.rhoE;

		state_topcell_1D.rho = state_topcell_2D.rho;
		double topcell_parallelvel = state_topcell_2D.rhou/state_topcell_2D.rho;
		state_topcell_1D.rhou = state_topcell_2D.rhov;
		state_topcell_1D.rhoE = state_topcell_2D.rhoE;


// Need to read in overall velocity magnitude to Exact Riemann Solver so we can get an accurate pressure reading!
// Reconstruct, from the state_[  ] ODstates, the TDstates corresponding assuming that the non-OD direction velocity is constant
// Then compute the flux from these new TDstates (4x), and update solution from t+1/2 to t+1 !!!
		ODstate state_leftedge_1D = Exact_Riemann_Solver(state_leftcell_1D, center_lr, leftcell_parallelvel, center_lr_parallelvel, thresh, gamma, debug);
		ODstate state_rightedge_1D = Exact_Riemann_Solver(center_lr, state_rightcell_1D, center_lr_parallelvel, rightcell_parallelvel, thresh, gamma, debug);
		ODstate state_bottomedge_1D = Exact_Riemann_Solver(state_bottomcell_1D, center_bt, bottomcell_parallelvel, center_bt_parallelvel, thresh, gamma, debug);
		ODstate state_topedge_1D = Exact_Riemann_Solver(center_bt, state_topcell_1D, center_bt_parallelvel, topcell_parallelvel, thresh, gamma, debug);

// initialize new TD states that will be computed from the riemann 1D states
		TDstate state_leftedge_2D, state_rightedge_2D, state_bottomedge_2D, state_topedge_2D;

		state_leftedge_2D.rho = state_leftedge_1D.rho;
		state_leftedge_2D.rhou = state_leftedge_1D.rhou;
		state_leftedge_2D.rhov = state_center_2D.rhov;
		state_leftedge_2D.rhoE = state_leftedge_1D.rhoE;
					
		state_rightedge_2D.rho = state_rightedge_1D.rho;
		state_rightedge_2D.rhou = state_rightedge_1D.rhou;
		state_rightedge_2D.rhov = state_center_2D.rhov;
		state_rightedge_2D.rhoE = state_rightedge_1D.rhoE;

		state_bottomedge_2D.rho = state_bottomedge_1D.rho;
		state_bottomedge_2D.rhou = state_center_2D.rhou; 
		state_bottomedge_2D.rhov = state_bottomedge_1D.rhou;
		state_bottomedge_2D.rhoE = state_bottomedge_1D.rhoE;

		state_topedge_2D.rho = state_topedge_1D.rho;
		state_topedge_2D.rhou = state_center_2D.rhou; 
		state_topedge_2D.rhov = state_topedge_1D.rhou;
		state_topedge_2D.rhoE = state_topedge_1D.rhoE;
			
		compute_cell_flux(tempflux, state_leftedge_2D, 0, gamma);
		F[0] = tempflux;
		compute_cell_flux(tempflux, state_rightedge_2D, 0, gamma);
		F[1] = tempflux;
		compute_cell_flux(tempflux, state_bottomedge_2D, 1, gamma);
		G[0] = tempflux;
		compute_cell_flux(tempflux, state_topedge_2D, 1, gamma);
		G[1] = tempflux;
	}
}



// Computes the TDstate "flux" from the Rusanov Flux Function
TDstate Rusanov(TDstate& state_minus, TDstate& state_plus, std::vector<double>& x_unit_normals, std::vector<double>& y_unit_normals, double& gamma, int edge_number) {
		double H_minus, H_plus, roe_average_H, roe_average_vel, roe_soundspeed;
		int direction;
		TDstate flux_minus, flux_plus;

// Determine the direction - 0 for LR and 1 for BT
		if ((edge_number%2) == 0) {
			direction = 0;
		} else {
			direction = 1;
		}
	
		double vel_normal_minus = (state_minus.rhou/state_minus.rho)*x_unit_normals[edge_number] + (state_minus.rhov/state_minus.rho)*y_unit_normals[edge_number];
		double vel_normal_plus = (state_plus.rhou/state_plus.rho)*x_unit_normals[edge_number] + (state_plus.rhov/state_plus.rho)*y_unit_normals[edge_number];

		// Compute enthalpy
		H_minus = state_minus.rhoE + compute_pressure_2D(state_minus, gamma)/state_minus.rho;
		H_plus = state_plus.rhoE + compute_pressure_2D(state_plus, gamma)/state_plus.rho;
		
		// Compute roe-averaged states
		roe_average_vel = (sqrt(state_minus.rho)*vel_normal_minus + sqrt(state_plus.rho)*vel_normal_plus)/(sqrt(state_minus.rho) + sqrt(state_plus.rho)); 
		roe_average_H = (sqrt(state_minus.rho)*H_minus + sqrt(state_plus.rho)*H_plus)/(sqrt(state_minus.rho) + sqrt(state_plus.rho));
		roe_soundspeed = sqrt((gamma-1)*(roe_average_H - 0.5*pow(roe_average_vel,2)));
		
		// Compute all eigenvalues (wavespeeds) for L cell, R cell and roe state
		std::vector<double> eigenvalues;	
		eigenvalues.push_back(roe_average_vel);
		eigenvalues.push_back(roe_average_vel + roe_soundspeed);
		eigenvalues.push_back(roe_average_vel - roe_soundspeed);

		double soundspeed, vel_minus, vel_plus;
		if (direction == 1) {
			vel_minus = state_minus.rhou/state_minus.rho;
			vel_plus  = state_plus.rhou/state_plus.rho;
		} else {
			vel_minus = state_minus.rhov/state_minus.rho;
			vel_plus  = state_plus.rhov/state_plus.rho;
		}
		soundspeed = sqrt(gamma*compute_pressure_2D(state_minus, gamma)/state_minus.rho);
		eigenvalues.push_back(vel_minus);
		eigenvalues.push_back(vel_minus + soundspeed);
		eigenvalues.push_back(vel_minus - soundspeed);

		soundspeed = sqrt(gamma*compute_pressure_2D(state_plus, gamma)/state_minus.rho);
		eigenvalues.push_back(vel_plus);
		eigenvalues.push_back(vel_plus + soundspeed);
		eigenvalues.push_back(vel_plus - soundspeed);

		// Entropy Fix - see "handout" from Prof. Fidkowski in Aero 623 Notes folder, pg. 4
		double epsilon = roe_soundspeed/100;
		for (unsigned int i = 0; i < eigenvalues.size(); i++) {
			if (std::abs(eigenvalues.at(i)) < epsilon) {
				//eigenvalues.at(i) = (pow(epsilon,2)+pow(eigenvalues.at(i),2))/(2*epsilon);
				eigenvalues.at(i) = eigenvalues.at(i) + 0.5*(pow(eigenvalues.at(i),2)/pow(epsilon,2) + epsilon);
			}
		}
		
		bool should_I_compute_max_abs_value_of_vector = true;
		double max_wavespeed = vectormax(eigenvalues, should_I_compute_max_abs_value_of_vector);
		
		// Compute flux values in the state_minus and state_plus TDstates
		compute_cell_flux(flux_minus, state_minus, direction, gamma);
		compute_cell_flux(flux_plus,  state_plus,  direction, gamma);

		// Compute Flux (F or G)
		TDstate flux;
		flux.rho  = 0.5*(flux_minus.rho + flux_plus.rho)   - 0.5*max_wavespeed*(state_plus.rho  - state_minus.rho);
		flux.rhou = 0.5*(flux_minus.rhou + flux_plus.rhou) - 0.5*max_wavespeed*(state_plus.rhou - state_minus.rhou);
		flux.rhov = 0.5*(flux_minus.rhov + flux_plus.rhov) - 0.5*max_wavespeed*(state_plus.rhov - state_minus.rhov);
		flux.rhoE = 0.5*(flux_minus.rhoE + flux_plus.rhoE) - 0.5*max_wavespeed*(state_plus.rhoE - state_minus.rhoE);

	return(flux);
}

// Returns bool 1 if grid[cellposition], (NOT CELL INDEX NUMBER), is a -2 edge (i.e. farfield/zero gradient BC ghost cell)
bool isedge(const std::vector<cell>& grid, int cellposition) {
	int special_value = -2;
	if ((grid[cellposition].adjacent_cells[0] == special_value) || (grid[cellposition].adjacent_cells[1] == special_value) || (grid[cellposition].adjacent_cells[2] == special_value)  || (grid[cellposition].adjacent_cells[3] == special_value)) {
		return(1);
	} else {
		return(0);
	}
}


// Assigns two TDstates "state_minus" and "state_plus" to the left and right OR top and bottom states based on the connectivity defined in the grid. 
//This is called WAY TOO MANY TIMES and there has to be a better way to do this! Lots of repeating tasks
// direction == 0 if left-right, direction == 1 if bottom-top
void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<TDstate>& U, const cell& current_cell, int direction) {

	int adjacent_bl_cellindex, adjacent_tr_cellindex;

	adjacent_bl_cellindex = current_cell.adjacent_cells_gridpos[direction]; // bottom (1) or left (0) cell location in grid vector
	adjacent_tr_cellindex = current_cell.adjacent_cells_gridpos[direction+2]; // Top (3) or right (2) cell location in grid vector

	state_minus = U[adjacent_bl_cellindex];
	state_plus = U[adjacent_tr_cellindex];
}


//  Computes the slope between two cells' TDstate and saves it to TDstate "slope_minus" and "slope_plus." Direction is input as well for either left-right (direction == 0) or bottom-top (direction == 1)
void compute_slope(TDstate &slope_minus, TDstate &slope_plus, int direction, TDstate state_minus, TDstate state_plus, directional_quantity delta, TDstate& U_center) {
	if (direction == 0) { // left-right
		slope_minus.rho = ((U_center.rho - state_minus.rho)/delta.left); // was 2*delta_x, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta.left);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta.left);
		slope_minus.rhoE = ((U_center.rhoE - state_minus.rhoE)/delta.left);

		slope_plus.rho = ((state_plus.rho - U_center.rho)/delta.right);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta.right);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta.right);
		slope_plus.rhoE = ((state_plus.rhoE - U_center.rhoE)/delta.right);
	} else { // bottom-top
		slope_minus.rho = ((U_center.rho - state_minus.rho)/delta.bottom); // was 2*delta_y, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta.bottom);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta.bottom);
		slope_minus.rhoE = ((U_center.rhoE - state_minus.rhoE)/delta.bottom);

		slope_plus.rho = ((state_plus.rho - U_center.rho)/delta.top);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta.top);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta.top);
		slope_plus.rhoE = ((state_plus.rhoE - U_center.rhoE)/delta.top);
	}
}


// The wrapper for all limiters - takes in user-defined limiter number and calls the appropriate limiter function
void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number) {
	
	double dd = 1E-20;
	TDstate r;
	r.rho = slope_minus.rho/(slope_plus.rho + dd);
	r.rhou = slope_minus.rhou/(slope_plus.rhou + dd);
	r.rhov = slope_minus.rhov/(slope_plus.rhov + dd);
	r.rhoE = slope_minus.rhoE/(slope_plus.rhoE + dd);

	if (limiter_number == 1) { // Harmonic
		harmonic_limiter(limiter_value, r);
	} else if (limiter_number == 2) { // First order upwind
		limiter_value.rho = 0;
		limiter_value.rhou = 0;
		limiter_value.rhov = 0;
		limiter_value.rhoE = 0;
	} else if (limiter_number == 3) { // Lax-Wendroff
		limiter_value = slope_plus;
	} else if (limiter_number == 4) { // Van Leers
		van_leers_limiter(limiter_value, r);
	} else if (limiter_number == 5) { // Superbee
		superbee_limiter(limiter_value, r);
	}
}


// Harmonic Limiter - Computes TDstate "limiter_value" for the harmonic limiter
void harmonic_limiter(TDstate &limiter_value, TDstate r) {
/*	double delta = pow(10,-12);
	limiter_value.rho = (std::abs(slope_minus.rho)*slope_plus.rho + slope_minus.rho*std::abs(slope_plus.rho))/(std::abs(slope_minus.rho) + std::abs(slope_plus.rho) + delta);
	limiter_value.rhou = (std::abs(slope_minus.rhou)*slope_plus.rhou + slope_minus.rhou*std::abs(slope_plus.rhou))/(std::abs(slope_minus.rhou) + std::abs(slope_plus.rhou) + delta);
	limiter_value.rhov = (std::abs(slope_minus.rhov)*slope_plus.rhov + slope_minus.rhov*std::abs(slope_plus.rhov))/(std::abs(slope_minus.rhov) + std::abs(slope_plus.rhov) + delta);
	limiter_value.rhoE = (std::abs(slope_minus.rhoE)*slope_plus.rhoE + slope_minus.rhoE*std::abs(slope_plus.rhoE))/(std::abs(slope_minus.rhoE) + std::abs(slope_plus.rhoE) + delta);
*/

	limiter_value.rho = (2*r.rho)/(1+r.rho);
	limiter_value.rhou = (2*r.rhou)/(1+r.rhou);
	limiter_value.rhov = (2*r.rhov)/(1+r.rhov);
	limiter_value.rhoE = (2*r.rhoE)/(1+r.rhoE);
}

// This computes the van leers limiter and returns TDstate limiter_value
void van_leers_limiter(TDstate &limiter_value, TDstate r) {

	limiter_value.rho  = (r.rho + std::abs(r.rho))/(1+std::abs(r.rho));
	limiter_value.rhou = (r.rhou + std::abs(r.rhou))/(1+std::abs(r.rhou));
	limiter_value.rhov = (r.rhov + std::abs(r.rhov))/(1+std::abs(r.rhov));
	limiter_value.rhoE = (r.rhoE + std::abs(r.rhoE))/(1+std::abs(r.rhoE));
}

void superbee_limiter(TDstate &limiter_value, TDstate r) {
	limiter_value.rho  = std::max(0.0, std::max(std::min(2*r.rho, 1.0), std::min(r.rho, 2.0)));
	limiter_value.rhou = std::max(0.0, std::max(std::min(2*r.rhou, 1.0), std::min(r.rhou, 2.0)));
	limiter_value.rhov = std::max(0.0, std::max(std::min(2*r.rhov, 1.0), std::min(r.rhov, 2.0)));
	limiter_value.rhoE = std::max(0.0, std::max(std::min(2*r.rhoE, 1.0), std::min(r.rhoE, 2.0)));
}

// Used to compute cell fluxes using the specified limiter for the F and G calculation to update from U(t) to U(t+1/2)
// Computes both the left and right fluxes for F and the bottom and top fluxes for G, returns each as [0] and [1] of fluxi and fluxj
// flux is for the (i-1/2 and i+1/2) faces [or (j-1/2 and j+1/2) faces]
void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, int& direction, double& gamma, double& CFL, bool& debug, std::ofstream& textout) {

	TDstate minus_limited_state, plus_limited_state, tempflux;

	assert(flux.size() == 2);

	minus_limited_state.rho  = center_state.rho  - 0.5*limiter_value.rho*(1-CFL);
	minus_limited_state.rhou = center_state.rhou - 0.5*limiter_value.rhou*(1-CFL);
	minus_limited_state.rhov = center_state.rhov - 0.5*limiter_value.rhov*(1-CFL);
	minus_limited_state.rhoE = center_state.rhoE - 0.5*limiter_value.rhoE*(1-CFL);

	plus_limited_state.rho  = center_state.rho  + 0.5*limiter_value.rho*(1-CFL);
	plus_limited_state.rhou = center_state.rhou + 0.5*limiter_value.rhou*(1-CFL);
	plus_limited_state.rhov = center_state.rhov + 0.5*limiter_value.rhov*(1-CFL);
	plus_limited_state.rhoE = center_state.rhoE + 0.5*limiter_value.rhoE*(1-CFL);

	if (debug && (direction == 0)) {
		textout << "Left Edge Reconstructed State: " << minus_limited_state.rho << " " << minus_limited_state.rhou << " " << minus_limited_state.rhov << " " << minus_limited_state.rhoE << "\n";
		textout << "Right Edge Reconstructed State: " << plus_limited_state.rho << " " << plus_limited_state.rhou << " " << plus_limited_state.rhov << " " << plus_limited_state.rhoE << "\n";
	}
	
	compute_cell_flux(tempflux, minus_limited_state, direction, gamma);
	flux[0] = tempflux;
	compute_cell_flux(tempflux, plus_limited_state, direction, gamma);
	flux[1] = tempflux;

}


// Computes the interface flux between two cells for a given TDstate "state" that is on the border between them. Used both for riemann fluxes and for interpolation fluxes
void compute_cell_flux(TDstate& flux, const TDstate& state, int direction, double& gamma) {
	assert(state.rho > 0);
	assert(state.rhoE > 0);

	double U = state.rhou/state.rho;
	double V = state.rhov/state.rho;
	double pressure = compute_pressure_2D(state, gamma);

/*	fluxi.rho  = (state.rhou)													*current_x_unit_normal*sign;
	fluxi.rhou = (state.rho*pow(U,2) + pressure)							*current_x_unit_normal*sign;
	fluxi.rhov = (state.rho*U*V)												*current_x_unit_normal*sign;
	fluxi.rhoE = (state.rho*U*(state.rhoE + pressure/state.rho))	*current_x_unit_normal*sign;

	fluxj.rho  = (state.rhov)													*current_y_unit_normal*sign;
	fluxj.rhou = (state.rho*U*V)												*current_y_unit_normal*sign;
	fluxj.rhov = (state.rho*pow(V,2) + pressure)							*current_y_unit_normal*sign;
	fluxj.rhoE = (state.rho*V*(state.rhoE + pressure/state.rho))	*current_y_unit_normal*sign;
*/

	if (direction == 0) { //direction == 0, left/right
		flux.rho  = state.rhou;
		flux.rhou = state.rho*pow(U,2) + pressure;
		flux.rhov = state.rho*U*V;
		flux.rhoE = state.rho*U*(state.rhoE + pressure/state.rho);

	} else { //direction == 1, bottom/top
		flux.rho  = state.rhov;
		flux.rhou = state.rho*V*U;
		flux.rhov = state.rho*pow(V,2) + pressure;
		flux.rhoE = state.rho*V*(state.rhoE + pressure/state.rho);
	}
}



// Computes the intermediate U(t+1/2) state, using the fluxes computed from the interpolation step.
// This is prior to the computation of the riemann flux
// THIS COULD BE A SOURCE OF ERROR - ORIGINALLY WAS DT/(2DX) CHANGED TO DT/(DELTA.RIGHT + DELTA.LEFT)!!!
void compute_halfway_state(TDstate& Uph, const TDstate& current_U, const std::vector<TDstate>& F, const std::vector<TDstate>& G, double& dt, directional_quantity& delta) {// cell& current_cell) {
	
//	double cell_area = current_cell.area;

	Uph.rho 	= (current_U.rho 	- dt/(delta.left + delta.right)*(F[1].rho-F[0].rho)   - dt/(delta.bottom + delta.top)*(G[1].rho-G[0].rho));
	Uph.rhou = (current_U.rhou - dt/(delta.left + delta.right)*(F[1].rhou-F[0].rhou) - dt/(delta.bottom + delta.top)*(G[1].rhou-G[0].rhou));
	Uph.rhov = (current_U.rhov - dt/(delta.left + delta.right)*(F[1].rhov-F[0].rhov) - dt/(delta.bottom + delta.top)*(G[1].rhov-G[0].rhov));
	Uph.rhoE = (current_U.rhoE - dt/(delta.left + delta.right)*(F[1].rhoE-F[0].rhoE) - dt/(delta.bottom + delta.top)*(G[1].rhoE-G[0].rhoE));
}


// Computes and returns the minimum spacing in the grid. For dt calculation using CFL number.
double gridmin(const std::vector<cell>& grid, char direction) {
	double dxy, min_dxy;
	if (direction == 'x') {
		min_dxy = std::abs(grid[0].cornerlocs_x[0] - grid[0].cornerlocs_x[1]);
		for (unsigned int it = 0; it < grid.size(); ++it) {
			dxy = std::abs(grid[it].cornerlocs_x[0] - grid[it].cornerlocs_x[1]);
			if (dxy < min_dxy) {
				min_dxy = dxy;
			}
			dxy = std::abs(grid[it].cornerlocs_x[3] - grid[it].cornerlocs_x[2]);
			if (dxy < min_dxy) {
				min_dxy = dxy;
			}
		}
	} else { // direction == "y"
		min_dxy = std::abs(grid[0].cornerlocs_y[3] - grid[0].cornerlocs_y[0]);
		for (unsigned int it = 0; it < grid.size(); ++it) {
			dxy = std::abs(grid[it].cornerlocs_y[3] - grid[it].cornerlocs_y[0]);
			if (dxy < min_dxy) {
				min_dxy = dxy;
			}
			dxy = std::abs(grid[it].cornerlocs_y[2] - grid[it].cornerlocs_y[1]);
			if (dxy < min_dxy) {
				min_dxy = dxy;
			}
		}
	}
	return(min_dxy);
}

void update_ghost_cells(std::vector<cell>& grid, std::vector<TDstate>& U, std::vector<TDstate>& Up1, std::vector<int>& interior_cells, bool debug, std::ofstream& textout) {
	for(unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		if (grid[cellposition].edge_type != 0) { // It is a ghost cell, either type 1, 2 or 3 

			int num_adjacent_interior_cells = 0;		
			std::vector<int> adjacent_edge; //vector of ints; 0 is left, 1 is bottom, 2 is right, 3 is top. Cardinal direction of adjacent cells
			std::vector<int> adjacent_interior_cellnumber;
			for(int adjacentcell_index = 0; adjacentcell_index < 4; ++adjacentcell_index) {
				auto interior_cell_loc = std::find(std::begin(interior_cells), std::end(interior_cells), grid[cellposition].adjacent_cells[adjacentcell_index]);
				if (interior_cell_loc != std::end(interior_cells)) { // if we actually find a match in the interior cells for the adjacent cell
					num_adjacent_interior_cells++;
					adjacent_edge.push_back(adjacentcell_index); //0 and is for l, 1 and is for bottom, [2] is r, [3] is top
					adjacent_interior_cellnumber.push_back(grid[cellposition].adjacent_cells[adjacentcell_index]);
				}
			}
			if (debug) {
				textout << "Cell Number: " << grid[cellposition].cellnumber << '\n';
				textout << "Num adjacent interior cells: " << num_adjacent_interior_cells << '\n';
			}
				// If this assert is tripped, your grid somehow has at least one cell that borders three or more interior cells, and we can't make boundary conditions for this!
			assert(num_adjacent_interior_cells < 3);
		
			double Uvel,Vvel;
			std::vector<int> directions = adjacent_edge;
				
			// Directions is vector of int of 0 or 1
			if (directions[0] >= 2) {
				directions[0] = (directions[0] - 2);
			}

			if (num_adjacent_interior_cells == 1) { // If it is a normal edge along a wall or far field, not a corner
				if (grid[cellposition].edge_type == 3) { // Pressure inlet - set conditions exactly equal to previous state's onditions
					Up1[cellposition] = U[cellposition];
					continue; //We are done with this cell!					
				}					
			
				// Assign the ghost cell state to be equal to its' adjacent cell's state
				Up1[cellposition] = Up1[find_cellposition(grid, grid[cellposition].adjacent_cells[adjacent_edge[0]])];

				Uvel = Up1[cellposition].rhou/Up1[cellposition].rho;
				Vvel = Up1[cellposition].rhov/Up1[cellposition].rho;			

//				std::vector<double> cartesian_x_unit_normal, cartesian_y_unit_normal;
//				compute_outward_unit_normal(grid[cellposition], cartesian_x_unit_normal, cartesian_y_unit_normal);
		
				if (grid[cellposition].edge_type == 1) { // wall BC, not farfield - set momentum equal and opposite to interior cell by setting a velocity reflection using V = V - 2*(V.n)n where every variable here is a vector. See 10/28/2015 - Making Walls w/ Actual Normal Vectors
					
					Up1[cellposition].rhou = (Uvel - 2*(Uvel*grid[cellposition].unit_normals_x[adjacent_edge[0]] + Vvel*grid[cellposition].unit_normals_y[adjacent_edge[0]])*grid[cellposition].unit_normals_x[adjacent_edge[0]])*Up1[cellposition].rho;
					Up1[cellposition].rhov = (Vvel - 2*(Uvel*grid[cellposition].unit_normals_x[adjacent_edge[0]] + Vvel*grid[cellposition].unit_normals_y[adjacent_edge[0]])*grid[cellposition].unit_normals_y[adjacent_edge[0]])*Up1[cellposition].rho;
				}
			} else { // This is an interior corner - two or more adjacent interior cells

				// Find where in the grid indexing the adjacent interior cells are, save to cellposition_adjacent1 and _adjacent2
				int cellposition_adjacent1 = find_cellposition(grid, grid[cellposition].adjacent_cells[directions[0]]);
				int cellposition_adjacent2 = find_cellposition(grid, grid[cellposition].adjacent_cells[directions[1]]);
	
				Up1[cellposition].rho = (Up1[cellposition_adjacent1].rho + Up1[cellposition_adjacent2].rho)/2;
				Up1[cellposition].rhoE = (Up1[cellposition_adjacent1].rhoE + Up1[cellposition_adjacent2].rhoE)/2;
					
				std::vector<int> directions01 = directions;

				if (directions01[0] >= 2) {
					directions01[0] = (directions01[0] - 2);
				}
				if (directions01[1] >= 2) {
					directions01[1] = (directions01[1] - 2);
				}

				// Here is where we actually update the ghost cell state
				//FIX THIS FOR NON-XY ONLY MESHES
				if (directions01[0] == 0) { //first one is the left-right boundary cell
					Uvel = -(Up1[cellposition_adjacent1].rhou/Up1[cellposition_adjacent1].rho);
					Vvel = -(Up1[cellposition_adjacent2].rhov/Up1[cellposition_adjacent2].rho);
				} else {
					Vvel = -(Up1[cellposition_adjacent1].rhov/Up1[cellposition_adjacent1].rho);
					Uvel = -(Up1[cellposition_adjacent2].rhou/Up1[cellposition_adjacent2].rho);
				}

				Up1[cellposition].rhou = Up1[cellposition].rho*Uvel;
				Up1[cellposition].rhov = Up1[cellposition].rho*Vvel;
			}

			if (debug) {
				textout << "Updated state: " << Up1[cellposition].rho << " " << Up1[cellposition].rhou << " " << Up1[cellposition].rhov << " " << Up1[cellposition].rhoE << '\n' << '\n';				
			}
		} // end of if (isedge(grid,cell)) 
	} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) (second one)
}
