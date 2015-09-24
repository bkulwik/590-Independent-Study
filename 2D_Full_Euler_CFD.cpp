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

void LimitedFluxSolver(const std::vector<cell>& grid, double& delta_x, double& delta_y, TDstate& state_minus, TDstate& state_plus, double gamma, TDstate& slope_plus, TDstate& slope_minus, int limiter_number, int cellposition, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double CFL, double& max_wavespeed, std::vector<TDstate>& U, bool debug, std::ofstream& textout);

// Define helper functions!
bool isedge(const std::vector<cell>& grid, int cellposition);

// void compute_flux(std::vector<double> &fluxph, std::vector<double> &fluxmh, std::vector<double> slope, TDstate state, double delta, std::string direction, double gamma);

// This is not math-y, just a formality of grabbing the left/right (state_minus) or bottom/top (state_plus) states from the grid input.
void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<TDstate>& U, const std::vector<cell>& grid, int cellposition, int direction);

// Computes the "geometric" slope vector for all conserved quantities between the bottom/left and top/right cells. Used for second order accuracy calculations.
void compute_slope(TDstate &slope_minus, TDstate &slope_plus, int direction, TDstate state_minus, TDstate state_plus, double delta_x, double delta_y, TDstate& U_center);

// This computes the double-value "limiter value" based on the bottom/left and top/right states for a given cell. The limiter value is then used to calculate the flux
void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number);

// Computes the harmonic limiter given r, the ratio of minus to plus slopes
void harmonic_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, TDstate r);

// Computes the van leers limiter 
void van_leers_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, TDstate r);

void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, const TDstate& slope_plus, const TDstate& slope_minus, int direction, double gamma, double CFL, bool debug, std::ofstream& textout);

// Computes the value of flux in the center cell, not along the edges but in the cell
void compute_cell_flux(TDstate &flux, const TDstate& state, int direction, double gamma);

// Reconstructs the solution from t to t+dt/2
void compute_halfway_state(TDstate& Uph, const TDstate& current_U, const std::vector<TDstate>& F, const std::vector<TDstate>& G, double dt, double dx, double dy);

// Computes the min dx or dy in a grid
double gridmin(const std::vector<cell>& grid, char direction);

// Find the position in the grid vector of cellnumber_desired
double find_cellposition(const std::vector<cell>& grid, int& cellnumber_desired);

// Outputs the grid to the file "grid_file.txt"
void output_grid(const std::vector<cell>& grid);




int main() {

// Define Variables
	std::string parameter_filename = "solver_parms.txt";

	std::vector<std::string> parameters;
	std::vector<double> Fiph, Fimh, Giph, Gimh;
	double dt, delta_x, delta_y, min_delta_x, min_delta_y, min_delta;
	TDstate state_minus, state_plus, slope_minus, slope_plus, limiter_value;
	TDstate state_leftcell_2D, state_rightcell_2D, state_bottomcell, state_topcell;
	std::vector<TDstate> F (2);
	std::vector<TDstate> G (2);

	double thresh = pow(10,-6); // 1E-8 tolerance for riemann solver - velocity
	double gamma = 1.4;

	bool save_timesteps = true;	

// 1 = harmonic mean for slopes, 2 = first order upwind (no limiter), 3 = lax-wendroff ...
// 2 is the only one that currently works - need to fix limiters; by timestep 2 we have E < 0 assert error (1) or rho < 0 assert error (3)
//	int limiter_number = 2; 

	std::vector<cell> grid;

	std::ofstream textout;
	textout.open("temp.txt");

	// Read parameters from input file defined by filename string
	parameters = read_parameters(parameter_filename);
	double CFL = std::stod(parameters.at(0));	
	int num_timesteps = std::stod(parameters.at(1));
	std::string input_filename = parameters.at(2);
	std::string output_filename = parameters.at(3);
	
	int debug;	
	if ((parameters.at(4).compare("True")) == 0) {
		debug = 1;
	} else {
		debug = 0;
	}
	std::cout << "Debug Mode: " << debug << "\n";
	int limiter_number = std::stoi(parameters.at(5));
	


// Read initial file defined from Matlab with cell edges and connections
// Read the initial condition file from MATLAB with [rho rho*u rho*v E] defined for each cell
	read_grid(input_filename, grid, gamma);
	
// Generate a vector of ints containing the cell numbers that are interior cells (0, not 1 or 2)
	std::vector<int> interior_cells = find_interior_cells(grid);

// Compute the min grid spacing for use in timestep calculations
	min_delta_x = gridmin(grid, 'x');
	min_delta_y = gridmin(grid, 'y');
	min_delta = std::min(min_delta_x,min_delta_y);

// Option to output grid to "grid_file.txt"
	output_grid(grid);

// Initialize solution U as vector of TD states with initial state at all points on grid
	std::vector<TDstate> U;
	for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		U.push_back(grid[cellposition].state);
	}


// Based on initial conditions, compute fastest wavespeed between any two cells in order to compute dt
// Using LimitedFLuxSolver is overkill here, but it also computes the max wavespeed in addition to all other calculations it is doing. Therefore it is used here for simplicity and to not necessitate rewriting code.

	double max_wavespeed = 0; 	
	for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) { // For each cell
		if (grid[cellposition].edge_type == 0) { // Cell is not a ghost cell/edge cell
         LimitedFluxSolver(grid, delta_x, delta_y, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellposition, limiter_value, F, G, CFL, max_wavespeed, U, debug, textout);
		}
	}

	dt = CFL*min_delta/max_wavespeed;
	
	if (debug) {
		textout << "\nInitial wavespeeds and dt computed, moving on to flux calculation \n \n";
	}

	TDstate Uph;
	std::vector<TDstate> Up1 (grid.size());




// Start calculation in for loop going from t = 0 to final time
	for (int timestep = 0; timestep < num_timesteps; timestep++) { // for each timestep
		std::cout << "Timestep " << timestep << " of " << num_timesteps << '\n';
		if (timestep > 0) {
			U = Up1;
		}

// Go through each cell, calculate fluxes, update to new piecewise linear state
		for (unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) { // For each cell
			double max_wavespeed = 0;
			if (grid[cellposition].edge_type == 0)  { // Cell is not a ghost cell/edge cell - only do calculations on interior cells

				//Compute F and G for the interpolation step - all else is junk needed for calculation of F and G
   	   	LimitedFluxSolver(grid, delta_x, delta_y, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellposition, limiter_value, F, G, CFL, max_wavespeed, U, debug, textout);

			// Now, use the fluxes calculated above to assign a new TDstate for the current cell, Uph
				compute_halfway_state(Uph, U[cellposition], F, G, dt, delta_x, delta_y);

				if (debug) {				
/*					std::cout << "Cell number: " << grid[cellposition].cellnumber << '\n';
					std::cout << "bottom/top Limiter - " << limiter_value.rho << " " << limiter_value.rhou << " " << limiter_value.rhov << " "  << limiter_value.E << '\n';
					std::cout << "Current State: " << U[cellposition].rho << " " << U[cellposition].rhou << " " << U[cellposition].rhov << " " << U[cellposition].E << '\n';
					std::cout << "F, cell number " << grid[cellposition].cellnumber << ": " << '\n';
					std::cout << F[0].rho << " " << F[0].rhou << " " << F[0].rhov << " " << F[0].E << '\n';
					std::cout << F[1].rho << " " << F[1].rhou << " " << F[1].rhov << " " << F[1].E << '\n';
	
					std::cout << "G, cell number: " << grid[cellposition].cellnumber << '\n';
					std::cout << G[0].rho << " " << G[0].rhou << " " << G[0].rhov << " " << G[0].E << '\n';
					std::cout << G[1].rho << " " << G[1].rhou << " " << G[1].rhov << " " << G[1].E << '\n';
	
					std::cout << "Cell " << grid[cellposition].cellnumber << " half-state: " <<  Uph.rho << " " << Uph.rhou << " " << Uph.rhov << " " << Uph.E << '\n';
*/
					textout << "Cell " << grid[cellposition].cellnumber << " Original State: " <<  U[cellposition].rho << " " << U[cellposition].rhou << " " << U[cellposition].rhov << " " << U[cellposition].E << '\n';
					textout << "F[0] = " << F[0].rho << " " << F[0].rhou << " " << F[0].rhov << " " << F[0].E << '\n';
					textout << "F[1] = " << F[1].rho << " " << F[1].rhou << " " << F[1].rhov << " " << F[1].E << '\n';
					textout << "G[0] = " << G[0].rho << " " << G[0].rhou << " " << G[0].rhov << " " << G[0].E << '\n';
					textout << "G[1] = " << G[0].rho << " " << G[1].rhou << " " << G[1].rhov << " " << G[1].E << '\n';
					textout << "Cell " << grid[cellposition].cellnumber << " half-state: " <<  Uph.rho << " " << Uph.rhou << " " << Uph.rhov << " " << Uph.E << '\n';
				}

// Solve riemann problem on edge for euler flux
				TDstate state_center_2D = U[cellposition];
				TDstate tempflux;
				ODstate state_leftcell_1D, center_lr, center_bt, state_rightcell_1D, state_bottomcell_1D, state_topcell_1D;
				
				// Assign what states are located on the left and right/bottom and top of the current cell. Easier to grab them here and store in a variable changed each iteration than later
				assign_edge_state(state_leftcell_2D, state_rightcell_2D, U, grid, cellposition, 0);
				assign_edge_state(state_bottomcell, state_topcell, U, grid, cellposition, 1); 

// All this mess is because the riemann problem is 1D and we store 2D states.
// So we have to make new, temporary 1D states and fixes (all the doubles) to make the Riemann Problem work.
				center_lr.rho = U[cellposition].rho;
				center_lr.rhou = U[cellposition].rhou;
				double center_lr_parallelvel = U[cellposition].rhov/U[cellposition].rho; 
				center_lr.E = U[cellposition].E;

				center_bt.rho = U[cellposition].rho;
				double center_bt_parallelvel = U[cellposition].rhou/U[cellposition].rho; 
				center_bt.rhou = U[cellposition].rhov;	
				center_bt.E = U[cellposition].E;
				
				state_leftcell_1D.rho = state_leftcell_2D.rho;
				state_leftcell_1D.rhou = state_leftcell_2D.rhou;
				double leftcell_parallelvel = state_leftcell_2D.rhov/state_leftcell_2D.rho;
				state_leftcell_1D.E = state_leftcell_2D.E;

				state_rightcell_1D.rho = state_rightcell_2D.rho; 
				state_rightcell_1D.rhou = state_rightcell_2D.rhou;
				double rightcell_parallelvel = state_rightcell_2D.rhov/state_rightcell_2D.rho; 
				state_rightcell_1D.E = state_rightcell_2D.E;

				state_bottomcell_1D.rho = state_bottomcell.rho;
				double bottomcell_parallelvel = state_bottomcell.rhou/state_bottomcell.rho;
				state_bottomcell_1D.rhou = state_bottomcell.rhov;
				state_bottomcell_1D.E = state_bottomcell.E;

				state_topcell_1D.rho = state_topcell.rho;
				double topcell_parallelvel = state_topcell.rhou/state_topcell.rho;
				state_topcell_1D.rhou = state_topcell.rhov;
				state_topcell_1D.E = state_topcell.E;


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
				state_leftedge_2D.E = state_leftedge_1D.E;
				
				state_rightedge_2D.rho = state_rightedge_1D.rho;
				state_rightedge_2D.rhou = state_rightedge_1D.rhou;
				state_rightedge_2D.rhov = state_center_2D.rhov;
				state_rightedge_2D.E = state_rightedge_1D.E;

				state_bottomedge_2D.rho = state_bottomedge_1D.rho;
				state_bottomedge_2D.rhou = state_center_2D.rhou; 
				state_bottomedge_2D.rhov = state_bottomedge_1D.rhou;
				state_bottomedge_2D.E = state_bottomedge_1D.E;

				state_topedge_2D.rho = state_topedge_1D.rho;
				state_topedge_2D.rhou = state_center_2D.rhou; 
				state_topedge_2D.rhov = state_topedge_1D.rhou;
				state_topedge_2D.E = state_topedge_1D.E;
				
				compute_cell_flux(tempflux, state_leftedge_2D, 0, gamma);
				F[0] = tempflux;
				compute_cell_flux(tempflux, state_rightedge_2D, 0, gamma);
				F[1] = tempflux;
				compute_cell_flux(tempflux, state_bottomedge_2D, 1, gamma);
				G[0] = tempflux;
				compute_cell_flux(tempflux, state_topedge_2D, 1, gamma);
				G[1] = tempflux;
//Use all these fluxes to define updated state on each cell 
				compute_halfway_state(Up1[cellposition], Uph, F, G, dt, delta_x, delta_y); 
				
				if (debug) {
					textout << "Cell " << grid[cellposition].cellnumber << " Full-Updated State : " <<  Up1[cellposition].rho << " " << Up1[cellposition].rhou << " " << Up1[cellposition].rhov << " " << Up1[cellposition].E << '\n' << '\n';
				}
			} // end of if (grid[cellposition].edge_type == 0) - loop through all interior cells
		} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) - loop through all timesteps

		//Update dt based on current wavespeeds
		dt = CFL*min_delta/max_wavespeed;

		std::cout << "dt = " << dt << '\n';



		// Now that the updated solution has been calculated, we need to redefine the ghost cells to the correct values
		for(unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
			if (grid[cellposition].edge_type != 0) { // It is a ghost cell, either type 1 or 2 

				int num_adjacent_interior_cells = 0;		
				std::vector<int> adjacent_directions; //vector of ints; 0 is left, 1 is bottom, 2 is right, 3 is top. Cardinal direction of adjacent cells
				std::vector<int> adjacent_interior_cellnumber;
				for(int adjacentcell_index = 0; adjacentcell_index < 4; ++adjacentcell_index) {
					auto interior_cell_loc = std::find(std::begin(interior_cells), std::end(interior_cells), grid[cellposition].adjacent_cells[adjacentcell_index]);
					if (interior_cell_loc != std::end(interior_cells)) { // if we actually find a match in the interior cells for the adjacent cell
						num_adjacent_interior_cells++;
						adjacent_directions.push_back(adjacentcell_index); //0 is l, 1 is bottom, 2 is r, 3 is top
						adjacent_interior_cellnumber.push_back(grid[cellposition].adjacent_cells[adjacentcell_index]);
					}
				}
				if (debug) {
					textout << "Cell Number: " << grid[cellposition].cellnumber << '\n';
					textout << "Num adjacent interior cells: " << num_adjacent_interior_cells << '\n';
				}

				// If this assert is tripped, your grid somehow has at least one cell that borders three or more interior cells, and we can't make boundary conditions for this!
				assert(num_adjacent_interior_cells < 3);
			
				if (num_adjacent_interior_cells == 1) { // If it is a normal edge along a wall or far field, not a corner

				// Assign the ghost cell state to be equal to its' adjacent cell's state
					Up1[cellposition] = Up1[find_cellposition(grid, grid[cellposition].adjacent_cells[adjacent_directions[0]])];

					if (adjacent_directions[0] >= 2) {
						adjacent_directions[0] = (adjacent_directions[0] - 2);
					}

					assert((adjacent_directions[0] == 0) || (adjacent_directions[0] == 1)); //should have been set just above

					if (grid[cellposition].edge_type == 1) { // wall BC, not farfield - set momentum equal and opposite to interior cell
						if ((adjacent_directions[0] == 0)){ // left-right
							if (Up1[cellposition].rhou != 0) {
								Up1[cellposition].rhou = -(Up1[cellposition].rhou);
							}
						} else if (adjacent_directions[0] == 1) { //up-down
							if (Up1[cellposition].rhov != 0) {
								Up1[cellposition].rhov = -(Up1[cellposition].rhov);
							}
						}			
					}
				} else { // This is an interior corner - two or more adjacent interior cells

					// Find where in the grid indexing the adjacent interior cells are, save to cellposition_adjacent1 and _adjacent2
					int cellposition_adjacent1 = find_cellposition(grid, grid[cellposition].adjacent_cells[adjacent_directions[0]]);
					int cellposition_adjacent2 = find_cellposition(grid, grid[cellposition].adjacent_cells[adjacent_directions[1]]);
	
					Up1[cellposition].rho = (Up1[cellposition_adjacent1].rho + Up1[cellposition_adjacent2].rho)/2;
					Up1[cellposition].E = (Up1[cellposition_adjacent1].E + Up1[cellposition_adjacent2].E)/2;
					
					std::vector<int> adjacent_directions01 = adjacent_directions;
					double U,V;

					if (adjacent_directions01[0] >= 2) {
						adjacent_directions01[0] = (adjacent_directions01[0] - 2);
					}
					if (adjacent_directions01[1] >= 2) {
						adjacent_directions01[1] = (adjacent_directions01[1] - 2);
					}

					// Here is where we actually update the ghost cell state
					if (adjacent_directions01[0] == 0) { //first one is the left-right boundary cell
						U = -(Up1[cellposition_adjacent1].rhou/Up1[cellposition_adjacent1].rho);
						V = -(Up1[cellposition_adjacent2].rhov/Up1[cellposition_adjacent2].rho);
					} else {
						V = -(Up1[cellposition_adjacent1].rhov/Up1[cellposition_adjacent1].rho);
						U = -(Up1[cellposition_adjacent2].rhou/Up1[cellposition_adjacent2].rho);
					}

					Up1[cellposition].rhou = Up1[cellposition].rho*U;
					Up1[cellposition].rhov = Up1[cellposition].rho*V;
				}

				if (debug) {
					textout << "Updated state: " << Up1[cellposition].rho << " " << Up1[cellposition].rhou << " " << Up1[cellposition].rhov << " " << Up1[cellposition].E << '\n' << '\n';
				}
			} // end of if (isedge(grid,cell)) 
		} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) (second one)

		if (save_timesteps) {
			write_to_file(grid,Up1,output_filename);
		}
	} // end of for(int t = 0; t <= 2*dt; t++) 

	textout.close();

return 0;

}



// This is the land of the helper functions! ---------------------------------------------------


// LimitedFluxSolver computes the flux for the initial t -> t+1/2 update
// Everything is passed by reference, only items that matter for returning are F and G, which are two-element vectors of TDstates
// This is computed prior to riemann problem for each timestep update

void LimitedFluxSolver(const std::vector<cell>& grid, double& delta_x, double& delta_y, TDstate& state_minus, TDstate& state_plus, double gamma, TDstate& slope_plus, TDstate& slope_minus, int limiter_number, int cellposition, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double CFL, double& max_wavespeed, std::vector<TDstate>& U, bool debug, std::ofstream& textout) {

	delta_x = (vectormax(grid[cellposition].cornerlocs_x) - vectormin(grid[cellposition].cornerlocs_x));
	delta_y = (vectormax(grid[cellposition].cornerlocs_y) - vectormin(grid[cellposition].cornerlocs_y));

   // Calculate the limited flux of each conserved quantity for initial t -> t+1/2 update
	for (int direction = 0; direction < 2; ++direction) { // For the left/right direction and bottom/top direction, two loop iterations only
		
		// Assign the bottom/left and the top/right states for current cell
		// state_minus and state_plus are the bottom/left and top/right TDstates, respectively
		// U is the vector of TDstates that is the current U for all states, not the initial U

		assign_edge_state(state_minus, state_plus, U, grid, cellposition, direction);	
		
		// Calculate the flux of each conserved variable, in given direction (0 L-R for first iteration, 1 B-T for second iteration), of that cell.
		// Flux needs to be calculated by first finding the limiter (direction-dependent), then using that limiter to calculate the directional flux, then using the two directional fluxes to update from U(t) to U(t+1/2)
			
		compute_slope(slope_minus, slope_plus, direction, state_minus, state_plus, delta_x, delta_y, U[cellposition]);
	
		if ((direction == 0) && (debug)) {
			textout << "Slope LR Minus: " << slope_minus.rho << " " << slope_minus.rhou << " " << slope_minus.rhov << " " << slope_minus.E << '\n';
			textout << "Slope LR Plus: " << slope_plus.rho << " " << slope_plus.rhou << " " << slope_plus.rhov << " " << slope_plus.E << '\n';
		}

//limiter_value is the B term on book pg. 238
		compute_limiter(limiter_value, slope_minus, slope_plus, limiter_number);

		if (debug) {
			textout << "Limiter Num: " << limiter_number << " Limiter Value: " << limiter_value.rho << " " << limiter_value.rhou << " " << limiter_value.rhov << " " << limiter_value.E << '\n';
		}

		if (direction == 0) { // Compute F, x-fluxes
			compute_limited_flux(F, limiter_value, U[cellposition], slope_plus, slope_minus, direction, gamma, CFL, debug, textout);
//			std::cout << "rho limiter x: " << limiter_value.rho << '\n';
		} else { // Compute G, y-fluxes
//			std::cout << "rho limiter y: " << limiter_value.rho << '\n';
			compute_limited_flux(G, limiter_value, U[cellposition], slope_plus, slope_minus, direction, gamma, CFL, debug, textout);
		}
	} // end of for (int direction = 0; direction < 2; ++direction)

	max_wavespeed_calculator_riemann(max_wavespeed, U[cellposition], gamma);
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


// Finds at what vector index number (cellposition) the desired cell number is located
double find_cellposition(const std::vector<cell>& grid, int& cellnumber_desired) {
	for(unsigned int cellposition = 0; cellposition < grid.size(); ++cellposition) {
		if (grid[cellposition].cellnumber == cellnumber_desired) {
			return(cellposition);
		}
	}
}


// Assigns two TDstates "state_minus" and "state_plus" to the left and right OR top and bottom states based on the connectivity defined in the grid. 
//This is called WAY TOO MANY TIMES and there has to be a better way to do this! Lots of repeating tasks
void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<TDstate>& U, const std::vector<cell>& grid, int cellposition, int direction) {

	int adjacent_bl, adjacent_tr;

	adjacent_bl = (grid[cellposition].adjacent_cells[direction]); //bottom or left cell index
	adjacent_tr = (grid[cellposition].adjacent_cells[direction+2]); //top or right cell index

	// Now find what the adjacent cell's cellposition is and grab its current state
	//cover simple case without a loop - if the bottom/left adjacent cell is the cell directly to the left (i.e. for left/right), grab its state; otherwise, go through a loop to find it (i.e. for bottom/top)

	if (grid[cellposition-1].cellnumber == adjacent_bl) { 
		state_minus = U[cellposition-1];
	} else { //go though a loop
		state_minus = U[find_cellposition(grid,adjacent_bl)];
// Can something be done here with this???
//auto interior_cell_loc = std::find(std::begin(interior_cells), std::end(interior_cells), grid[cellposition].adjacent_cells[adjacentcell_index]);

	}

//	std::cout << "In AES: adjacent_bl = " << adjacent_bl << '\n';	

	if (grid[cellposition+1].cellnumber == adjacent_tr) {
		state_plus = U[cellposition+1];
	} else { //go though a loop
		state_plus = U[find_cellposition(grid,adjacent_tr)];	
	}
}


//  Computes the slope between two cells' TDstate and saves it to TDstate "slope_minus" and "slope_plus." Direction is input as well for either left-right (direction == 0) or bottom-top (direction == 1)
void compute_slope(TDstate &slope_minus, TDstate &slope_plus, int direction, TDstate state_minus, TDstate state_plus, double delta_x, double delta_y, TDstate& U_center) {
	if (direction == 0) { // left-right
		slope_minus.rho = ((U_center.rho - state_minus.rho)/delta_x); // was 2*delta_x, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta_x);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta_x);
		slope_minus.E = ((U_center.E - state_minus.E)/delta_x);

		slope_plus.rho = ((state_plus.rho - U_center.rho)/delta_x);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta_x);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta_x);
		slope_plus.E = ((state_plus.E - U_center.E)/delta_x);
	} else { // bottom-top
		slope_minus.rho = ((U_center.rho - state_minus.rho)/delta_y); // was 2*delta_y, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta_y);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta_y);
		slope_minus.E = ((U_center.E - state_minus.E)/delta_y);

		slope_plus.rho = ((state_plus.rho - U_center.rho)/delta_y);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta_y);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta_y);
		slope_plus.E = ((state_plus.E - U_center.E)/delta_y);
	}
}


// The wrapper for all limiters - takes in user-defined limiter number and calls the appropriate limiter function
void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number) {
	
	TDstate r;
	r.rho = slope_minus.rho/slope_plus.rho;
	r.rhou = slope_minus.rhou/slope_plus.rhou;
	r.rhov = slope_minus.rhov/slope_plus.rhov;
	r.E = slope_minus.E/slope_plus.E;

	if (limiter_number == 1) { // Harmonic
		harmonic_limiter(limiter_value, slope_minus, slope_plus, r);
	} else if (limiter_number == 2) { // First order upwind
		limiter_value.rho = 0;
		limiter_value.rhou = 0;
		limiter_value.rhov = 0;
		limiter_value.E = 0;
	} else if (limiter_number == 3) { // Lax-Wendroff
		limiter_value = slope_plus;
	} else if (limiter_number == 4) { // Van Leers
		van_leers_limiter(limiter_value, slope_minus, slope_plus, r);
	}
}


// Harmonic Limiter - Computes TDstate "limiter_value" for the harmonic limiter
void harmonic_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, TDstate r) {
/*	double delta = pow(10,-12);
	limiter_value.rho = (std::abs(slope_minus.rho)*slope_plus.rho + slope_minus.rho*std::abs(slope_plus.rho))/(std::abs(slope_minus.rho) + std::abs(slope_plus.rho) + delta);
	limiter_value.rhou = (std::abs(slope_minus.rhou)*slope_plus.rhou + slope_minus.rhou*std::abs(slope_plus.rhou))/(std::abs(slope_minus.rhou) + std::abs(slope_plus.rhou) + delta);
	limiter_value.rhov = (std::abs(slope_minus.rhov)*slope_plus.rhov + slope_minus.rhov*std::abs(slope_plus.rhov))/(std::abs(slope_minus.rhov) + std::abs(slope_plus.rhov) + delta);
	limiter_value.E = (std::abs(slope_minus.E)*slope_plus.E + slope_minus.E*std::abs(slope_plus.E))/(std::abs(slope_minus.E) + std::abs(slope_plus.E) + delta);
*/

	limiter_value.rho = (2*r.rho)/(1+r.rho);
	if (isnan(limiter_value.rho)) {limiter_value.rho = 0;}
	limiter_value.rhou = (2*r.rhou)/(1+r.rhou);
	if (isnan(limiter_value.rhou)) {limiter_value.rhou = 0;}
	limiter_value.rhov = (2*r.rhov)/(1+r.rhov);
	if (isnan(limiter_value.rhov)) {limiter_value.rhov = 0;}
	limiter_value.E = (2*r.E)/(1+r.E);
	if (isnan(limiter_value.E)) {limiter_value.E = 0;}
}

// This computes the van leers limiter and returns TDstate limiter_value
void van_leers_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, TDstate r) {

	limiter_value.rho = (r.rho + std::abs(r.rho))/(1+std::abs(r.rho));
	limiter_value.rhou = (r.rhou + std::abs(r.rhou))/(1+std::abs(r.rhou));
	limiter_value.rhov = (r.rhov + std::abs(r.rhov))/(1+std::abs(r.rhov));
	limiter_value.E = (r.E + std::abs(r.E))/(1+std::abs(r.E));

	if (isnan(limiter_value.rho)) {limiter_value.rho = 0;}
	if (isnan(limiter_value.rhou)) {limiter_value.rhou = 0;}
	if (isnan(limiter_value.rhov)) {limiter_value.rhov = 0;}
	if (isnan(limiter_value.E)) {limiter_value.E = 0;}
}

// Used to compute cell fluxes using the specified limiter for the F and G calculation to update from U(t) to U(t+1/2)
void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, const TDstate& slope_plus, const TDstate& slope_minus, int direction, double gamma, double CFL, bool debug, std::ofstream& textout) {
// flux is for the (i-1/2 and i+1/2) faces [or (j-1/2 and j+1/2) faces]
	TDstate minus_limited_state, plus_limited_state, tempflux;

	assert(flux.size() == 2);

//	std::cout << "In CLF: center_state = " << center_state.rho << " " << center_state.rhou << " " << center_state.rhov << " " << center_state.E << '\n';
/*
	double a = CFL*0.625/0.00113534;//CFL*dx/dt;

	tempflux.rho = a*(center_state.rho + 1/2*limiter_value.rho*(1-CFL));
	tempflux.rhou = a*(center_state.rhou + 1/2*limiter_value.rhou*(1-CFL));
	tempflux.rhov = a*(center_state.rhov + 1/2*limiter_value.rhov*(1-CFL));
	tempflux.E = a*(center_state.E + 1/2*limiter_value.E*(1-CFL));
	flux[0] = tempflux;

	tempflux.rho = a*(center_state.rho - 1/2*limiter_value.rho*(1-CFL));
	tempflux.rhou = a*(center_state.rhou - 1/2*limiter_value.rhou*(1-CFL));
	tempflux.rhov = a*(center_state.rhov - 1/2*limiter_value.rhov*(1-CFL));
	tempflux.E = a*(center_state.E - 1/2*limiter_value.E*(1-CFL));
	flux[1] = tempflux;
*/


// This is the code that was used for correctly completing the first order calculatons for the 623 project.

	minus_limited_state.rho = center_state.rho - 0.5*slope_minus.rho*limiter_value.rho*(1-CFL);
	minus_limited_state.rhou = center_state.rhou - 0.5*slope_minus.rhou*limiter_value.rhou*(1-CFL);
	minus_limited_state.rhov = center_state.rhov - 0.5*slope_minus.rhov*limiter_value.rhov*(1-CFL);
	minus_limited_state.E = center_state.E - 0.5*slope_minus.E*limiter_value.E*(1-CFL);

	plus_limited_state.rho = center_state.rho + 0.5*slope_plus.rho*limiter_value.rho*(1-CFL);
	plus_limited_state.rhou = center_state.rhou + 0.5*slope_plus.rhou*limiter_value.rhou*(1-CFL);
	plus_limited_state.rhov = center_state.rhov + 0.5*slope_plus.rhov*limiter_value.rhov*(1-CFL);
	plus_limited_state.E = center_state.E + 0.5*slope_plus.E*limiter_value.E*(1-CFL);

	if (debug && (direction == 0)) {
		textout << "Left Edge Reconstructed State: " << minus_limited_state.rho << " " << minus_limited_state.rhou << " " << minus_limited_state.rhov << " " << minus_limited_state.E << "\n";
		textout << "Right Edge Reconstructed State: " << plus_limited_state.rho << " " << plus_limited_state.rhou << " " << plus_limited_state.rhov << " " << plus_limited_state.E << "\n";
	}

	compute_cell_flux(tempflux, minus_limited_state, direction, gamma);
	flux[0] = (tempflux);
	compute_cell_flux(tempflux, plus_limited_state, direction, gamma);
	flux[1] = (tempflux);

}



// Computes the interface flux between two cells for a given TDstate "state" that is on the border between them. Used both for riemann fluxes and for interpolation fluxes
void compute_cell_flux(TDstate &flux, const TDstate& state, int direction, double gamma) {

//	std::cout << "In CCF: state = " << state.rho << " " << state.rhou << " " << state.rhov << " " << state.E << '\n';
	assert(state.rho > 0);
	assert(state.E > 0);

	double U = state.rhou/state.rho;
	double V = state.rhov/state.rho;
	double pressure = compute_pressure_2D(state, gamma);
	
	if (direction == 0) { // left/right direction
		flux.rho = state.rhou;
		flux.rhou = state.rho*pow(U,2) + pressure;
		if ((U == 0) || (V == 0)) {
			flux.rhov = 0;
		} else {
			flux.rhov = state.rho*U*V;
		}
		flux.E = U*(state.E + pressure);
	} else { //direction == 1, bottom/top
		flux.rho = state.rhov;
		if ((U == 0) || (V == 0)) {
			flux.rhou = 0;
		} else {
			flux.rhou = state.rho*U*V;
		}
//		std::cout << "In CCF: rho,V, V^2, P = " << state.rho << " " << V << " " << pow(V,2) << " " << pressure << '\n';
		flux.rhov = state.rho*pow(V,2) + pressure;
		flux.E = V*(state.E + pressure);
	}
}



// Computes the intermediate U(t+1/2) state, using the fluxes computed from the interpolation step.
// This is prior to the computation of the riemann flux
void compute_halfway_state(TDstate& Uph, const TDstate& current_U, const std::vector<TDstate>& F, const std::vector<TDstate>& G, double dt, double dx, double dy) {
	Uph.rho = (current_U.rho - dt/(2*dx)*(F[1].rho-F[0].rho) - dt/(2*dy)*(G[1].rho-G[0].rho));
	Uph.rhou = (current_U.rhou - dt/(2*dx)*(F[1].rhou-F[0].rhou) - dt/(2*dy)*(G[1].rhou-G[0].rhou));
	Uph.rhov = (current_U.rhov - dt/(2*dx)*(F[1].rhov-F[0].rhov) - dt/(2*dy)*(G[1].rhov-G[0].rhov));
	Uph.E = (current_U.E - dt/(2*dx)*(F[1].E-F[0].E) - dt/(2*dy)*(G[1].E-G[0].E));
}


// Computes and returns the minimum spacing in the grid. For dt calculation using CFL number.
double gridmin(const std::vector<cell>& grid, char direction) {
	double delta, min_delta;
	if (direction == 'x') {
		min_delta = std::abs(grid[0].cornerlocs_x[0] - grid[0].cornerlocs_x[1]);
		for (unsigned int it = 0; it < grid.size(); ++it) {
			delta = std::abs(grid[it].cornerlocs_x[0] - grid[it].cornerlocs_x[1]);
			if (delta < min_delta) {
				min_delta = delta;
			}
			delta = std::abs(grid[it].cornerlocs_x[3] - grid[it].cornerlocs_x[2]);
			if (delta < min_delta) {
				min_delta = delta;
			}
		}
	} else { // direction == "y"
		min_delta = std::abs(grid[0].cornerlocs_y[3] - grid[0].cornerlocs_y[0]);
		for (unsigned int it = 0; it < grid.size(); ++it) {
			delta = std::abs(grid[it].cornerlocs_y[3] - grid[it].cornerlocs_y[0]);
			if (delta < min_delta) {
				min_delta = delta;
			}
			delta = std::abs(grid[it].cornerlocs_y[2] - grid[it].cornerlocs_y[1]);
			if (delta < min_delta) {
				min_delta = delta;
			}
		}
	}
	return(min_delta);
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
		ofile << input.state.E << '\n';
	}
}

