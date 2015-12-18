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

// Compute the first instance of dt, changing values of dt, cell_perimeter and cell_diameter
void initialize_dt(double& dt, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, double gamma, std::vector<double>& cell_perimeter, std::vector<double>& cell_diameter);

// 
void Reconstruction_Flux(const std::vector<cell>& grid, double& gamma, int& limiter_number, unsigned int& cellposition, std::vector<TDstate>& cell_x_flux, std::vector<TDstate>& cell_y_flux, double& CFL, std::vector<TDstate>& U, bool& debug, std::ofstream& textout);

// Create the correct input and solve the exact riemann problem
void Solve_Riemann_Flux( int& flux_function_choice, unsigned int& cellposition, double& gamma, double& thresh, bool& debug, std::vector<TDstate>& Uph, cell& current_cell, std::vector<TDstate>& F, std::vector<TDstate>& G );

// Update ghost cells to reflect for wall boundaries
void update_ghost_cells(std::vector<cell>& grid, std::vector<TDstate>& old_state_vector, std::vector<TDstate>& new_state_vector, std::vector<int>& interior_cells, bool debug, std::ofstream& textout, bool& freestream_case);

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
void harmonic_limiter(TDstate& limiter_value, TDstate& slope_minus, TDstate& slope_plus);

// Computes the van leers limiter 
void van_leers_limiter(TDstate &limiter_value, TDstate r);

// Computes the Superbee limiter
void superbee_limiter(TDstate &limiter_value, TDstate r);

// Computes the limited flux for two cells (left-right or bottom-top) of the reconstruction step
void compute_limited_flux(TDstate& cell_directional_flux, const TDstate& limiter_value, const TDstate& center_state, int direction, double& gamma, double& CFL, bool& debug, std::ofstream& textout, double x_unit_normal, double y_unit_normal);

// Computes the 2D flux given a TDstate "state" on the edge
void compute_cell_edge_flux(TDstate& F, TDstate& G, const TDstate& state, double& gamma);

// Reconstructs the solution from t to t+dt/2
void compute_halfway_state(TDstate& Uph, TDstate& residual, const TDstate& current_U, const std::vector<TDstate>& cell_x_flux, const std::vector<TDstate>& cell_y_flux, double& dt, cell& current_cell);

// Computes the min dx or dy (linear) in a grid
//double gridmin(const std::vector<cell>& grid, char direction);



int main() {

// Define Variables
	std::string parameter_filename = "solver_parms.txt";

	std::vector<std::string> parameters;
	std::vector<double> Fiph, Fimh, Giph, Gimh;
	std::vector<TDstate> cell_x_flux (2);
	std::vector<TDstate> cell_y_flux (2);

	double thresh = pow(10,-6); // 1E-6 tolerance for exact riemann solver - velocity
	double gamma = 1.4;

	bool save_timesteps = true;	


	// Read parameters from input file defined by filename string
	parameters = read_parameters(parameter_filename);
	double CFL = std::stod(parameters.at(0));	
	int num_timesteps = std::stod(parameters.at(1));
	std::string input_filename = parameters.at(2);
	std::string output_filename = parameters.at(3);
	
	bool debug, freestream_case;	
	if ((parameters.at(4).compare("True")) == 0) {
		freestream_case = 1;
	} else {
		freestream_case = 0;
	}
std::cout << "Freestream test? " << freestream_case << '\n';
	if ((parameters.at(5).compare("True")) == 0) {
		debug = 1;
	} else {
		debug = 0;
	}
	std::cout << "Debug Mode: " << debug << "\n";
	std::ofstream textout;
	if (debug) {
		textout.open("debug_file.txt");
	}
	int flux_function_choice = std::stoi(parameters.at(6));
	int limiter_number = std::stoi(parameters.at(7));
	
	


// Read initial file defined from Matlab with cell edges and connections
// Read the initial condition file from MATLAB with [rho rho*u rho*v E] defined for each cell
	std::vector<cell> grid;	

	read_grid(input_filename, grid, gamma);
	
// Generate a vector of ints containing the cell numbers that are interior cells (0, not 1 or 2)
	std::vector<int> interior_cells = find_interior_cells(grid);

// Option to output grid to "grid_file.txt"
	output_grid(grid);

// Initialize solution U as vector of TD states with initial state at all points on grid
	std::vector<TDstate> U;
	for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		U.push_back(grid[cellposition].state);
	}

//	textout << "Updating ghost cells...\n";

// Before starting calculation, make sure ghost cells have been correctly initialized to correct normal condition
// This is not done during grid generation - ghost and interior cells will have the same state. We need to recompute the ghost cells to make sure walls are correctly accounted for.
//	update_ghost_cells(grid, U, U, interior_cells, debug, textout);
//	textout << "Ghost cells updated!\n";	

// Initialize needed variables
	TDstate Uph;
	std::vector<TDstate> Up1 (grid.size());
	std::vector<TDstate> Uph_vector (grid.size());
	std::vector<TDstate> residuals;
	std::vector<double> cell_perimeter (grid.size()), cell_diameter (grid.size());
	for (unsigned int i = 0; i <= grid.size(); ++i) {
		residuals.push_back(Uph); //just an empty TDstate 
		residuals[i].rho = 0;
		residuals[i].rhou = 0;
		residuals[i].rhov = 0;
		residuals[i].rhoE = 0;
	}

// Initialize the dt value
	double dt, t = 0;
	initialize_dt(dt, grid, U, CFL, gamma, cell_perimeter, cell_diameter);

	if (debug) {
		textout << "\nInitial wavespeeds and dt computed, moving on to flux calculation \n \n";
	}
	std::cout << "\nInitial wavespeeds and dt computed, moving on to flux calculation \n \n";

// Start calculation in for loop going from t = 0 to final time
	for (int timestep = 0; timestep < num_timesteps; timestep++) { // for each timestep
		std::cout << "Timestep " << timestep+1 << " of " << num_timesteps << ", dt = " << dt << '\n';;
		if (timestep > 0) {
			U = Up1;
		}

	// Go through each cell, calculate fluxes, update to new piecewise linear state, find flux and update to state Uph
		for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) { // For each cell
			if (grid[cellposition].edge_type == 0)  { // Cell is not a ghost cell/edge cell - only do calculations on interior cells

			// Compute cell_x_flux and cell_y_flux for the interpolation step - all else is junk needed for calculation of cell_x_flux and cell_y_flux
   	   	Reconstruction_Flux(grid, gamma, limiter_number, cellposition, cell_x_flux, cell_y_flux, CFL, U, debug, textout);

			// Now, use the fluxes calculated above to assign a new TDstate for the current cell, Uph
				compute_halfway_state(Uph, residuals[grid[cellposition].cellnumber], U[cellposition], cell_x_flux, cell_y_flux, dt, grid[cellposition]);

				Uph_vector[cellposition] = Uph;
				
				if (debug) {
					textout << "Cell " << grid[cellposition].cellnumber << " Original State: " <<  U[cellposition].rho << " " << U[cellposition].rhou << " " << U[cellposition].rhov << " " << U[cellposition].rhoE << '\n';
					textout << "Reconstruction: \ncell_x_flux[0] = " << cell_x_flux[0].rho << " " << cell_x_flux[0].rhou << " " << cell_x_flux[0].rhov << " " << cell_x_flux[0].rhoE << '\n';
					textout << "cell_x_flux[1] = " << cell_x_flux[1].rho << " " << cell_x_flux[1].rhou << " " << cell_x_flux[1].rhov << " " << cell_x_flux[1].rhoE << '\n';
					textout << "cell_y_flux[0] = " << cell_y_flux[0].rho << " " << cell_y_flux[0].rhou << " " << cell_y_flux[0].rhov << " " << cell_y_flux[0].rhoE << '\n';
					textout << "cell_y_flux[1] = " << cell_y_flux[1].rho << " " << cell_y_flux[1].rhou << " " << cell_y_flux[1].rhov << " " << cell_y_flux[1].rhoE << '\n';
					textout << "Cell " << grid[cellposition].cellnumber << " half-state: " <<  Uph.rho << " " << Uph.rhou << " " << Uph.rhov << " " << Uph.rhoE << "\n\n";
				}
			} else {
				Uph_vector[cellposition] = U[cellposition];
			}
		}

		output_residual_sum(residuals);
		if (limiter_number != 1) { // If using a higher order metho
			update_ghost_cells(grid, U, Uph_vector, interior_cells, debug, textout, freestream_case);
		}

	// Go through each cell, solve riemann problem (exact or not exact), calculate fluxes and update to state Up1
		for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) { // For each cell
			if (grid[cellposition].edge_type == 0)  { // Cell is not a ghost cell/edge cell - only do calculations on interior cells
			// Solve riemann problem on edge with specified flux function for euler flux
				Solve_Riemann_Flux(flux_function_choice, cellposition, gamma, thresh, debug, Uph_vector, grid[cellposition], cell_x_flux, cell_y_flux);

			// Use all these fluxes to define updated state on each cell 
				compute_halfway_state(Up1[cellposition], residuals[grid[cellposition].cellnumber], Uph_vector[cellposition], cell_x_flux, cell_y_flux, dt, grid[cellposition]); 
		
				if (debug) {
					textout << "Current Cell: " << grid[cellposition].cellnumber << '\n';
					textout << "Halfway State = " << Uph_vector[cellposition].rho << " " << Uph_vector[cellposition].rhou << " " << Uph_vector[cellposition].rhov << " " << Uph_vector[cellposition].rhoE << '\n';
					textout << "Riemann Fluxes: \ncell_x_flux[0] = " << cell_x_flux[0].rho << " " << cell_x_flux[0].rhou << " " << cell_x_flux[0].rhov << " " << cell_x_flux[0].rhoE << '\n';
					textout << "cell_x_flux[1] = " << cell_x_flux[1].rho << " " << cell_x_flux[1].rhou << " " << cell_x_flux[1].rhov << " " << cell_x_flux[1].rhoE << '\n';
					textout << "cell_y_flux[0] = " << cell_y_flux[0].rho << " " << cell_y_flux[0].rhou << " " << cell_y_flux[0].rhov << " " << cell_y_flux[0].rhoE << '\n';
					textout << "cell_y_flux[1] = " << cell_y_flux[1].rho << " " << cell_y_flux[1].rhou << " " << cell_y_flux[1].rhov << " " << cell_y_flux[1].rhoE << '\n';
					textout << "Cell " << grid[cellposition].cellnumber << " Full-Updated State : " <<  Up1[cellposition].rho << " " << Up1[cellposition].rhou << " " << Up1[cellposition].rhov << " " << Up1[cellposition].rhoE << "\n\n";
				}
				assert(Up1[cellposition].rho >= 0);
				assert(Up1[cellposition].rhoE >= 0);
			} // end of if (grid[cellposition].edge_type == 0) - loop through all interior cells
		} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) - loop through all timesteps


	// Now that the updated solution has been calculated, we need to redefine the ghost cells to the correct values
		update_ghost_cells(grid, Uph_vector, Up1, interior_cells, debug, textout, freestream_case);

		//if ((save_timesteps) && ((timestep%10 == 0) || timestep == num_timesteps)) {
		if (save_timesteps) {
			write_to_file(grid,Up1,output_filename);
		}
		
		//Update dt based on new wavespeeds
		dt = timestep_calculator(gamma, grid, Up1, CFL, cell_diameter, cell_perimeter);
		
		textout << '\n';
		t = t + dt;
		output_residual_sum(residuals);
	} // end of for(int t = 0; t <= 2*dt; t++) 
	textout.close();
	return 0;
}



// This is the land of the helper functions! ---------------------------------------------------

void initialize_dt(double& dt, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, double gamma, std::vector<double>& cell_perimeter, std::vector<double>& cell_diameter) {

	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell		

/*		cell_edge_lengths[cellnum].bottom 	= sqrt(pow((grid[cellnum].cornerlocs_x[1]-grid[cellnum].cornerlocs_x[0]),2) + pow((grid[cellnum].cornerlocs_y[1]-grid[cellnum].cornerlocs_y[0]),2));
		cell_edge_lengths[cellnum].right 	= sqrt(pow((grid[cellnum].cornerlocs_x[2]-grid[cellnum].cornerlocs_x[1]),2) + pow((grid[cellnum].cornerlocs_y[2]-grid[cellnum].cornerlocs_y[1]),2));
		cell_edge_lengths[cellnum].top 		= sqrt(pow((grid[cellnum].cornerlocs_x[3]-grid[cellnum].cornerlocs_x[2]),2) + pow((grid[cellnum].cornerlocs_y[3]-grid[cellnum].cornerlocs_y[2]),2));
		cell_edge_lengths[cellnum].left 		= sqrt(pow((grid[cellnum].cornerlocs_x[0]-grid[cellnum].cornerlocs_x[3]),2) + pow((grid[cellnum].cornerlocs_y[0]-grid[cellnum].cornerlocs_y[3]),2));
*/

		cell_perimeter[cellnum] = (grid[cellnum].edge_lengths.left + grid[cellnum].edge_lengths.bottom + grid[cellnum].edge_lengths.right + grid[cellnum].edge_lengths.top);
		cell_diameter[cellnum] = (4*grid[cellnum].area/cell_perimeter[cellnum]);
		//std::cout << "cell " << grid[cellnum].cellnumber << " area: " << grid[cellnum].area << ", perimeter: " << cell_perimeter[cellnum] << '\n';
	}
	
	// Now that we have cell_diameter and perimeter, compute timestep	
	dt = timestep_calculator(gamma, grid, Up1, CFL, cell_diameter, cell_perimeter);

}


// Reconstruction_Flux computes the flux for the initial t -> t+1/2 update
// Everything is passed by reference, only items that matter for returning are F and G, which are two-element vectors of TDstates
// This is computed prior to riemann problem for each timestep update from the state limited reconstruction
void Reconstruction_Flux(const std::vector<cell>& grid, double& gamma, int& limiter_number, unsigned int& cellposition, std::vector<TDstate>& cell_x_flux, std::vector<TDstate>& cell_y_flux, double& CFL, std::vector<TDstate>& U, bool& debug, std::ofstream& textout) {
	TDstate state_minus, state_plus, slope_minus, slope_plus, limiter_value_LR, limiter_value_BT;

   // Calculate the limited flux of each conserved quantity for initial t -> t+1/2 update

	for (int direction = 0; direction < 2; ++direction) { // For each of the two directions - LR and BT - compute the slope between neighboring states

		// Grabs the bottom/left and the top/right states for current cell, depending on direction
		// state_minus and state_plus are the bottom/left and top/right TDstates, respectively
		// U is the vector of TDstates that is the current U for all states, not the initial U
		assign_edge_state(state_minus, state_plus, U, grid[cellposition], direction);	
		
		// Calculate the flux of each conserved variable, in given direction (0 L-R for first iteration, 1 B-T for second iteration), of that cell.
		// Flux needs to be calculated by first finding the limiter (direction-dependent), then using that limiter to calculate the directional flux, then using the two directional fluxes to update from U(t) to U(t+1/2)
			
		compute_slope(slope_minus, slope_plus, direction, state_minus, state_plus, grid[cellposition].cell_distance, U[cellposition]);
	
//limiter_value is the B term on book pg. 238
		if (direction == 0) {
			compute_limiter(limiter_value_LR, slope_minus, slope_plus, limiter_number);
		} else { 
			compute_limiter(limiter_value_BT, slope_minus, slope_plus, limiter_number);
		}

		if ((direction == 0) && (debug)) {
			textout << "Slope LR Minus: " << slope_minus.rho << " " << slope_minus.rhou << " " << slope_minus.rhov << " " << slope_minus.rhoE << '\n';
			textout << "Slope LR Plus: " << slope_plus.rho << " " << slope_plus.rhou << " " << slope_plus.rhov << " " << slope_plus.rhoE << '\n';
			textout << "Limiter Num: " << limiter_number << " Limiter Value (LR): " << limiter_value_LR.rho << " " << limiter_value_LR.rhou << " " << limiter_value_LR.rhov << " " << limiter_value_LR.rhoE << '\n';	
			textout << "Weighting flux with areas (L,R): " << grid[cellposition].edge_lengths.left << " " << grid[cellposition].edge_lengths.right << '\n';
			//textout << "x unit normals (left edge, right edge): " << x_unit_normals_onedirection[0] << " " << x_unit_normals_onedirection[1] << '\n';
			//textout << "y unit normals (left edge, right edge): " << y_unit_normals_onedirection[0] << " " << y_unit_normals_onedirection[1] << '\n';	
		}
		if ((direction == 1) && (debug)) {
			textout << "Slope BT Minus: " << slope_minus.rho << " " << slope_minus.rhou << " " << slope_minus.rhov << " " << slope_minus.rhoE << '\n';
			textout << "Slope BT Plus: " << slope_plus.rho << " " << slope_plus.rhou << " " << slope_plus.rhov << " " << slope_plus.rhoE << '\n';
			textout << "Limiter Num: " << limiter_number << " Limiter Value (BT): " << limiter_value_BT.rho << " " << limiter_value_BT.rhou << " " << limiter_value_BT.rhov << " " << limiter_value_BT.rhoE << '\n';		
			textout << "Weighting flux with areas (B,T): " << grid[cellposition].edge_lengths.bottom << " " << grid[cellposition].edge_lengths.top << '\n';
			//textout << "x unit normals (bottom edge, top edge): " << x_unit_normals_onedirection[0] << " " << x_unit_normals_onedirection[1] << '\n';
			//textout << "y unit normals (bottom edge, top edge): " << y_unit_normals_onedirection[0] << " " << y_unit_normals_onedirection[1] << '\n';
		}
	}

	// For each edge, compute the flux
	// For the left and right edges:		
	compute_limited_flux(cell_x_flux[0], limiter_value_LR, U[cellposition], 0, gamma, CFL, debug, textout, grid[cellposition].outward_unit_normals_x[0], grid[cellposition].outward_unit_normals_y[0]);
	compute_limited_flux(cell_x_flux[1], limiter_value_LR, U[cellposition], 2, gamma, CFL, debug, textout, grid[cellposition].outward_unit_normals_x[2], grid[cellposition].outward_unit_normals_y[2]);
	// For the bottom and top edges:
	compute_limited_flux(cell_y_flux[0], limiter_value_BT, U[cellposition], 1, gamma, CFL, debug, textout, grid[cellposition].outward_unit_normals_x[1], grid[cellposition].outward_unit_normals_y[1]);
	compute_limited_flux(cell_y_flux[1], limiter_value_BT, U[cellposition], 3, gamma, CFL, debug, textout, grid[cellposition].outward_unit_normals_x[3], grid[cellposition].outward_unit_normals_y[3]);

	//Now, we need to weight the cell_x_flux and cell_y_flux based on cell edge areas!
	weight_flux_with_area(cell_x_flux[0], grid[cellposition].edge_lengths.left);
	weight_flux_with_area(cell_x_flux[1], grid[cellposition].edge_lengths.right);
	weight_flux_with_area(cell_y_flux[0], grid[cellposition].edge_lengths.bottom);
	weight_flux_with_area(cell_y_flux[1], grid[cellposition].edge_lengths.top);
}


void Solve_Riemann_Flux( int& flux_function_choice, unsigned int& cellposition, double& gamma, double& thresh, bool& debug, std::vector<TDstate>& Uph, cell& current_cell, std::vector<TDstate>& cell_x_flux, std::vector<TDstate>& cell_y_flux ) {
	// debugging for flux errors	
	//std::cout << "Cell Number: " << current_cell.cellnumber << '\n';
	std::vector<double> x_unit_normals = current_cell.outward_unit_normals_x;
	std::vector<double> y_unit_normals = current_cell.outward_unit_normals_y;

	TDstate state_center_2D = Uph[cellposition];
	TDstate state_leftcell_2D, state_rightcell_2D, state_bottomcell_2D, state_topcell_2D;
	ODstate state_leftcell_1D, center_lr, center_bt, state_rightcell_1D, state_bottomcell_1D, state_topcell_1D;
				
	// Assign what states are located on the left and right/bottom and top of the current cell. Easier to grab them here and store in a variable changed each iteration than later
	assign_edge_state(state_bottomcell_2D, state_topcell_2D, Uph, current_cell, 1); //direction 1 == bottom-top

	assign_edge_state(state_leftcell_2D, state_rightcell_2D, Uph, current_cell, 0); //direction 0 == left-right
	if (flux_function_choice == 1) { //Roe Flux

	} else if (flux_function_choice == 2) { // Rusanov Flux - this takes cell normals into account! Yey!	
		cell_x_flux[0] = Rusanov(state_center_2D, state_leftcell_2D,	x_unit_normals, y_unit_normals, gamma, 0); 	// Edge 0 - left
		cell_x_flux[1] = Rusanov(state_center_2D, state_rightcell_2D, 	x_unit_normals, y_unit_normals, gamma, 2); 	// Edge 2 - right 
		cell_y_flux[0] = Rusanov(state_center_2D, state_bottomcell_2D, x_unit_normals, y_unit_normals, gamma, 1);	// Edge 1 - bottom
		cell_y_flux[1] = Rusanov(state_center_2D, state_topcell_2D,   	x_unit_normals, y_unit_normals, gamma, 3); 	// Edge 3 - top
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

		TDstate F,G;
		
		compute_cell_edge_flux(F, G, state_leftedge_2D, gamma); //Computes the F and G for the state on the left edge
		compute_edge_normal_flux(F, G, x_unit_normals[0], y_unit_normals[0], cell_x_flux[0]); // computes the cell's directional left edge flux, cell_x_flux[0] 

		compute_cell_edge_flux(F, G, state_rightedge_2D, gamma);
		compute_edge_normal_flux(F, G, x_unit_normals[2], y_unit_normals[2], cell_x_flux[1]);

		compute_cell_edge_flux(F, G, state_bottomedge_2D, gamma);
		compute_edge_normal_flux(F, G, x_unit_normals[1], y_unit_normals[1], cell_y_flux[0]);

		compute_cell_edge_flux(F, G, state_topedge_2D, gamma);
		compute_edge_normal_flux(F, G, x_unit_normals[3], y_unit_normals[3], cell_y_flux[1]);
	}

// Use area weighting on each of the fluxes to account for the area of the face it is coming into or out of
	weight_flux_with_area(cell_x_flux[0], current_cell.edge_lengths.left);
	weight_flux_with_area(cell_x_flux[1], current_cell.edge_lengths.right);
	weight_flux_with_area(cell_y_flux[0], current_cell.edge_lengths.bottom);
	weight_flux_with_area(cell_y_flux[1], current_cell.edge_lengths.top);

	//std::cout << current_cell.edge_lengths.left << " " << current_cell.edge_lengths.right << " " << current_cell.edge_lengths.bottom << " " << current_cell.edge_lengths.top << '\n';
	//std::cout << "\nUnit Normals X: " << x_unit_normals[0] << " " << x_unit_normals[1] << " " << x_unit_normals[2] << " " << x_unit_normals[3] << '\n';
	//std::cout << "Unit Normals Y: " << y_unit_normals[0] << " " << y_unit_normals[1] << " " << y_unit_normals[2] << " " << y_unit_normals[3] << '\n';

}


// Computes the TDstate "flux" from the Rusanov Flux Function
// Left is the center cell ALWAYS, right is the adjacent cell ALWAYS - per outward unit normal convention
TDstate Rusanov(TDstate& state_left, TDstate& state_right, std::vector<double>& x_unit_normals, std::vector<double>& y_unit_normals, double& gamma, int edge_number) {
		double H_left, H_right, roe_average_H, roe_average_vel_x, roe_average_vel_y, roe_average_vel, roe_soundspeed;
	//	int direction;
		TDstate F_left, G_left, F_right, G_right;

/*// Determine the direction - 0 for LR and 1 for BT
		if ((edge_number%2) == 0) {
			direction = 0;
		} else {
			direction = 1;
		}*/
	
		double vel_normal_left = (state_left.rhou/state_left.rho)*x_unit_normals[edge_number] + (state_left.rhov/state_left.rho)*y_unit_normals[edge_number];
		double vel_normal_right = (state_right.rhou/state_right.rho)*x_unit_normals[edge_number] + (state_right.rhov/state_right.rho)*y_unit_normals[edge_number];

		// Compute enthalpy
		H_left = state_left.rhoE + compute_pressure_2D(state_left, gamma)/state_left.rho;
		H_right = state_right.rhoE + compute_pressure_2D(state_right, gamma)/state_right.rho;
		
		// Compute roe-averaged states
		roe_average_vel_x = (sqrt(state_left.rho)*(state_left.rhou/state_left.rho) + sqrt(state_right.rho)*(state_right.rhou/state_right.rho))/(sqrt(state_left.rho) + sqrt(state_right.rho));
		roe_average_vel_y = (sqrt(state_left.rho)*(state_left.rhov/state_left.rho) + sqrt(state_right.rho)*(state_right.rhov/state_right.rho))/(sqrt(state_left.rho) + sqrt(state_right.rho));
		roe_average_H = (sqrt(state_left.rho)*H_left + sqrt(state_right.rho)*H_right)/(sqrt(state_left.rho) + sqrt(state_right.rho));
		roe_average_vel = sqrt(roe_average_vel_x*roe_average_vel_x + roe_average_vel_y*roe_average_vel_y);
		roe_soundspeed = sqrt((gamma-1)*(roe_average_H - 0.5*pow(roe_average_vel,2)));
		
		// Compute all eigenvalues (wavespeeds) for L cell, R cell and roe state
		std::vector<double> eigenvalues;	
		eigenvalues.push_back(roe_average_vel);
		eigenvalues.push_back(roe_average_vel + roe_soundspeed);
		eigenvalues.push_back(roe_average_vel - roe_soundspeed);

		double soundspeed;//, vel_left, vel_right;
	/*	if (direction == 1) {
			vel_left = (state_left.rhou/state_left.rho;
			vel_right  = state_right.rhou/state_right.rho;
		} else {
			vel_left = state_left.rhov/state_left.rho;
			vel_right  = state_right.rhov/state_right.rho;
		}*/
		soundspeed = sqrt(gamma*compute_pressure_2D(state_left, gamma)/state_left.rho);
	//	eigenvalues.push_back(vel_left);
	//	eigenvalues.push_back(vel_left + soundspeed);
	//	eigenvalues.push_back(vel_left - soundspeed);
		eigenvalues.push_back(vel_normal_left);
		eigenvalues.push_back(vel_normal_left + soundspeed);
		eigenvalues.push_back(vel_normal_left - soundspeed);

		soundspeed = sqrt(gamma*compute_pressure_2D(state_right, gamma)/state_left.rho);
//		eigenvalues.push_back(vel_right);
//		eigenvalues.push_back(vel_right + soundspeed);
//		eigenvalues.push_back(vel_right - soundspeed);
		eigenvalues.push_back(vel_normal_right);
		eigenvalues.push_back(vel_normal_right + soundspeed);
		eigenvalues.push_back(vel_normal_right - soundspeed);


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
		
		// Compute flux values in the state_left and state_right TDstates
		compute_cell_edge_flux(F_left, G_left, state_left, gamma);
		compute_cell_edge_flux(F_right,  G_right,  state_right,  gamma);

		// Compute Directional Flux on Cell Edge taking into account edge normals - see Notebook 2 pg. 41
		// This takes the F and G fluxes that were already computed and projects them onto the cell UNIT normal of the current cell.
		TDstate directional_flux_left, directional_flux_right, directional_flux_celledge;

		compute_edge_normal_flux(F_left, G_left, x_unit_normals[edge_number], y_unit_normals[edge_number], directional_flux_left);
		compute_edge_normal_flux(F_right, G_right, x_unit_normals[edge_number], y_unit_normals[edge_number], directional_flux_right);
		
		// debugging for flux errors
		//std::cout << "Edge Number: " << edge_number << '\n';
		//std::cout << "Edge Normals: x = " << x_unit_normals[edge_number] << '\n';
		//std::cout << "Edge Normals: y = " << y_unit_normals[edge_number] << '\n';

		// Compute the overall flux on this edge, NOT area-weighted here! DOES take normals, not areas, into account.
		directional_flux_celledge.rho  = 0.5*(directional_flux_left.rho  + directional_flux_right.rho)  - 0.5*max_wavespeed*(state_right.rho  - state_left.rho);
		directional_flux_celledge.rhou = 0.5*(directional_flux_left.rhou + directional_flux_right.rhou) - 0.5*max_wavespeed*(state_right.rhou - state_left.rhou);
		directional_flux_celledge.rhov = 0.5*(directional_flux_left.rhov + directional_flux_right.rhov) - 0.5*max_wavespeed*(state_right.rhov - state_left.rhov);
		directional_flux_celledge.rhoE = 0.5*(directional_flux_left.rhoE + directional_flux_right.rhoE) - 0.5*max_wavespeed*(state_right.rhoE - state_left.rhoE);

	return(directional_flux_celledge);
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
		slope_minus.rho  = ((U_center.rho  - state_minus.rho)/delta.left); // was 2*delta_x, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta.left);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta.left);
		slope_minus.rhoE = ((U_center.rhoE - state_minus.rhoE)/delta.left);

		slope_plus.rho  = ((state_plus.rho  - U_center.rho)/delta.right);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta.right);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta.right);
		slope_plus.rhoE = ((state_plus.rhoE - U_center.rhoE)/delta.right);
	} else { // bottom-top
		slope_minus.rho  = ((U_center.rho  - state_minus.rho)/delta.bottom);
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta.bottom);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta.bottom);
		slope_minus.rhoE = ((U_center.rhoE - state_minus.rhoE)/delta.bottom);

		slope_plus.rho  = ((state_plus.rho  - U_center.rho)/delta.top);
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

	if (limiter_number == 2) { // Harmonic
		harmonic_limiter(limiter_value, slope_minus, slope_plus);
	} else if (limiter_number == 1) { // First order upwind
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
//void harmonic_limiter(TDstate &limiter_value, TDstate r) {
void harmonic_limiter(TDstate& limiter_value, TDstate& slope_minus, TDstate& slope_plus) {

	double delta = pow(10,-12);
	limiter_value.rho = (std::abs(slope_minus.rho)*slope_plus.rho + slope_minus.rho*std::abs(slope_plus.rho))/(std::abs(slope_minus.rho) + std::abs(slope_plus.rho) + delta);
	limiter_value.rhou = (std::abs(slope_minus.rhou)*slope_plus.rhou + slope_minus.rhou*std::abs(slope_plus.rhou))/(std::abs(slope_minus.rhou) + std::abs(slope_plus.rhou) + delta);
	limiter_value.rhov = (std::abs(slope_minus.rhov)*slope_plus.rhov + slope_minus.rhov*std::abs(slope_plus.rhov))/(std::abs(slope_minus.rhov) + std::abs(slope_plus.rhov) + delta);
	limiter_value.rhoE = (std::abs(slope_minus.rhoE)*slope_plus.rhoE + slope_minus.rhoE*std::abs(slope_plus.rhoE))/(std::abs(slope_minus.rhoE) + std::abs(slope_plus.rhoE) + delta);

/*
	limiter_value.rho = (2*r.rho)/(1+r.rho);
	limiter_value.rhou = (2*r.rhou)/(1+r.rhou);
	limiter_value.rhov = (2*r.rhov)/(1+r.rhov);
	limiter_value.rhoE = (2*r.rhoE)/(1+r.rhoE);*/
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

// Used to compute cell fluxes using the specified limiter for the cell_x_flux and cell_y_flux calculation to update from U(t) to U(t+1/2)
// Note that this "limited flux" is computed ONLY for the limited reconstructions, not with riemann solver fluxes.
// This is called once for each edge. 
// right_state is named such as to remind us that the center cell is the "left" state and the neighboring cell is the "right" state for outward unit normals
// int direction is either [0 1 2 3] for [left bottom right top] direction
void compute_limited_flux(TDstate& cell_directional_flux, const TDstate& limiter_value, const TDstate& center_state, int direction, double& gamma, double& CFL, bool& debug, std::ofstream& textout, double x_unit_normal, double y_unit_normal) {

	TDstate limited_state, F, G;

	if ((direction == 0) || (direction == 2)) {
		limited_state.rho  = center_state.rho  - 0.5*limiter_value.rho*(1-CFL);
		limited_state.rhou = center_state.rhou - 0.5*limiter_value.rhou*(1-CFL);
		limited_state.rhov = center_state.rhov - 0.5*limiter_value.rhov*(1-CFL);
		limited_state.rhoE = center_state.rhoE - 0.5*limiter_value.rhoE*(1-CFL);
	} else if ((direction == 1) || (direction == 3)) {
		limited_state.rho  = center_state.rho  + 0.5*limiter_value.rho*(1-CFL);
		limited_state.rhou = center_state.rhou + 0.5*limiter_value.rhou*(1-CFL);
		limited_state.rhov = center_state.rhov + 0.5*limiter_value.rhov*(1-CFL);
		limited_state.rhoE = center_state.rhoE + 0.5*limiter_value.rhoE*(1-CFL);
	} else {
		assert(0 && "direction is not [0 1 2 3], exiting");
	}

	if (debug) {
		textout << "Edge " << direction << " Reconstructed State: " << limited_state.rho << " " << limited_state.rhou << " " << limited_state.rhov << " " << limited_state.rhoE << "\n";
	}
	
	// F and G are temporary variables
	// Compute the cell edge flux from the limited state
	compute_cell_edge_flux(F, G, limited_state, gamma);
	compute_edge_normal_flux(F, G, x_unit_normal, y_unit_normal, cell_directional_flux);
}


// Computes the interface flux between two cells for a given TDstate "state" that is on the border between them. Also used to compute the flux values for a given "state" that is a cell state, not an edge state. Used both for riemann fluxes and for interpolation fluxes
void compute_cell_edge_flux(TDstate& F, TDstate& G, const TDstate& state, double& gamma) {
	//assert(state.rho > 0);
	//assert(state.rhoE > 0);

	double U = state.rhou/state.rho;
	double V = state.rhov/state.rho;
	double pressure = compute_pressure_2D(state, gamma);

	F.rho  = state.rhou;
	F.rhou = state.rho*pow(U,2) + pressure;
	F.rhov = state.rho*U*V;
	F.rhoE = state.rho*U*(state.rhoE + pressure/state.rho);

	G.rho  = state.rhov;
	G.rhou = state.rho*V*U;
	G.rhov = state.rho*pow(V,2) + pressure;
	G.rhoE = state.rho*V*(state.rhoE + pressure/state.rho);
}



// Computes the intermediate U(t+1/2) state, using the fluxes computed from the interpolation step.
void compute_halfway_state(TDstate& Uph, TDstate& residual, const TDstate& current_U, const std::vector<TDstate>& cell_x_flux, const std::vector<TDstate>& cell_y_flux, double& dt, cell& current_cell) {
	
	double cell_area = current_cell.area;
	residual.rho  = cell_x_flux[1].rho  + cell_x_flux[0].rho  + cell_y_flux[1].rho  + cell_y_flux[0].rho;
	residual.rhou = cell_x_flux[1].rhou + cell_x_flux[0].rhou + cell_y_flux[1].rhou + cell_y_flux[0].rhou;
	residual.rhov = cell_x_flux[1].rhov + cell_x_flux[0].rhov + cell_y_flux[1].rhov + cell_y_flux[0].rhov;
	residual.rhoE = cell_x_flux[1].rhoE + cell_x_flux[0].rhoE + cell_y_flux[1].rhoE + cell_y_flux[0].rhoE;
	
	/*if ((current_cell.cellnumber) == 12) {
		std::cout << "rho-fluxes: " ;
		std::cout << cell_x_flux[1].rho << " " << cell_x_flux[0].rho << " " << cell_y_flux[1].rho << " " << cell_y_flux[0].rho <<'\n';
		std::cout << "These should be equal: " << cell_x_flux[1].rho  + cell_x_flux[0].rho  + cell_y_flux[1].rho  + cell_y_flux[0].rho << " " << residual.rho << '\n';
		std::cout << "residual: ";
		outputTDstate(residual);
		std::cout << "Half-Updated State for Cell 12: ";
		outputTDstate(Uph);
	}*/

// The actual flux update...
	Uph.rho 	= (current_U.rho 	- dt/(2*cell_area)*residual.rho);
	Uph.rhou	= (current_U.rhou - dt/(2*cell_area)*residual.rhou);
	Uph.rhov	= (current_U.rhov - dt/(2*cell_area)*residual.rhov);
	Uph.rhoE	= (current_U.rhoE - dt/(2*cell_area)*residual.rhoE);

 //Temporary flux update for debugging
/*	Uph.rho 	= (current_U.rho 	- dt/(cell_area)*residual.rho);
	Uph.rhou	= (current_U.rhou - dt/(cell_area)*residual.rhou);
	Uph.rhov	= (current_U.rhov - dt/(cell_area)*residual.rhov);
	Uph.rhoE	= (current_U.rhoE - dt/(cell_area)*residual.rhoE);

/*	Old flux update
	Uph.rho 	= (current_U.rho 	- dt/(2*cell_area)*(cell_x_flux[1].rho -cell_x_flux[0].rho)  - dt/(2*cell_area)*(cell_y_flux[1].rho -cell_y_flux[0].rho));
	Uph.rhou = (current_U.rhou - dt/(2*cell_area)*(cell_x_flux[1].rhou-cell_x_flux[0].rhou) - dt/(2*cell_area)*(cell_y_flux[1].rhou-cell_y_flux[0].rhou));
	Uph.rhov = (current_U.rhov - dt/(2*cell_area)*(cell_x_flux[1].rhov-cell_x_flux[0].rhov) - dt/(2*cell_area)*(cell_y_flux[1].rhov-cell_y_flux[0].rhov));
	Uph.rhoE = (current_U.rhoE - dt/(2*cell_area)*(cell_x_flux[1].rhoE-cell_x_flux[0].rhoE) - dt/(2*cell_area)*(cell_y_flux[1].rhoE-cell_y_flux[0].rhoE));
*/
}

/*
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
*/

// "new_state_vector" has old ghost cells and you want it updated to new ghost cells
void update_ghost_cells(std::vector<cell>& grid, std::vector<TDstate>& old_state_vector, std::vector<TDstate>& new_state_vector, std::vector<int>& interior_cells, bool debug, std::ofstream& textout, bool& freestream_case) {

	for(unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		if (grid[cellposition].edge_type != 0) { // It is a ghost cell, either type 1, 2 or 3 
		//	std::cout << "Updating ghost cell " << grid[cellposition].cellnumber << '\n';
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

			if (num_adjacent_interior_cells == 1) { // If it is a "nominal" edge along a wall or far field, not a corner
				if (grid[cellposition].edge_type == 3) { // Pressure inlet or freestream test - set conditions exactly equal to previous state's conditions
					new_state_vector[cellposition] = old_state_vector[cellposition];
					if (debug) {
						textout << "Updated state: " << new_state_vector[cellposition].rho << " " << new_state_vector[cellposition].rhou << " " << new_state_vector[cellposition].rhov << " " << new_state_vector[cellposition].rhoE << '\n' << '\n';				
					}
					continue; //We are done with this cell!					
				}					
			
				// Assign the ghost cell state to be equal to its' adjacent cell's state
				//Up1[cellposition] = Up1[grid[cellposition].adjacent_cells_gridpos[adjacent_edge[0]]];				
				new_state_vector[cellposition] = new_state_vector[find_cellposition(grid, grid[cellposition].adjacent_cells[adjacent_edge[0]])];

				if (freestream_case) continue;

				Uvel = new_state_vector[cellposition].rhou/new_state_vector[cellposition].rho;
				Vvel = new_state_vector[cellposition].rhov/new_state_vector[cellposition].rho;			

//				std::vector<double> cartesian_x_unit_normal, cartesian_y_unit_normal;
//				compute_outward_unit_normal(grid[cellposition], cartesian_x_unit_normal, cartesian_y_unit_normal);
		
				if (grid[cellposition].edge_type == 1) { // wall BC, not farfield - set momentum equal and opposite to interior cell by setting a velocity reflection using V = V - 2*(V.n)n where every variable here is a vector. See 10/28/2015 - Making Walls w/ Actual Normal Vectors
					new_state_vector[cellposition].rhou = (Uvel - 2*(Uvel*grid[cellposition].outward_unit_normals_x[adjacent_edge[0]] + Vvel*grid[cellposition].outward_unit_normals_y[adjacent_edge[0]])*grid[cellposition].outward_unit_normals_x[adjacent_edge[0]])*new_state_vector[cellposition].rho;
					new_state_vector[cellposition].rhov = (Vvel - 2*(Uvel*grid[cellposition].outward_unit_normals_x[adjacent_edge[0]] + Vvel*grid[cellposition].outward_unit_normals_y[adjacent_edge[0]])*grid[cellposition].outward_unit_normals_y[adjacent_edge[0]])*new_state_vector[cellposition].rho;
				}
			} else { // This is an interior corner - two or more adjacent interior cells

				// Find where in the grid indexing the adjacent interior cells are, save to cellposition_adjacent1 and _adjacent2
				int cellposition_adjacent1 = find_cellposition(grid, grid[cellposition].adjacent_cells[directions[0]]);
				int cellposition_adjacent2 = find_cellposition(grid, grid[cellposition].adjacent_cells[directions[1]]);
	
				new_state_vector[cellposition].rho = (new_state_vector[cellposition_adjacent1].rho + new_state_vector[cellposition_adjacent2].rho)/2;
				new_state_vector[cellposition].rhoE = (new_state_vector[cellposition_adjacent1].rhoE + new_state_vector[cellposition_adjacent2].rhoE)/2;
					
				std::vector<int> directions01 = directions;

				if (directions01[0] >= 2) {
					directions01[0] = (directions01[0] - 2);
				}
				if (directions01[1] >= 2) {
					directions01[1] = (directions01[1] - 2);
				}

				// Here is where we actually update the ghost cell state
				//FIX THIS FOR NON-XY ONLY MESHES - although interior corner not currently being used on non-xy meshes :)
				if (directions01[0] == 0) { //first one is the left-right boundary cell
					Uvel = -(new_state_vector[cellposition_adjacent1].rhou/new_state_vector[cellposition_adjacent1].rho);
					Vvel = -(new_state_vector[cellposition_adjacent2].rhov/new_state_vector[cellposition_adjacent2].rho);
				} else {
					Vvel = -(new_state_vector[cellposition_adjacent1].rhov/new_state_vector[cellposition_adjacent1].rho);
					Uvel = -(new_state_vector[cellposition_adjacent2].rhou/new_state_vector[cellposition_adjacent2].rho);
				}

				new_state_vector[cellposition].rhou = new_state_vector[cellposition].rho*Uvel;
				new_state_vector[cellposition].rhov = new_state_vector[cellposition].rho*Vvel;
			}

			if (debug) {
				textout << "Updated state: " << new_state_vector[cellposition].rho << " " << new_state_vector[cellposition].rhou << " " << new_state_vector[cellposition].rhov << " " << new_state_vector[cellposition].rhoE << '\n' << '\n';				
			}
		} // end of if (isedge(grid,cell)) 
	} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) (second one)
}
