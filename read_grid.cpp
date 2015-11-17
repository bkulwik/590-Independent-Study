#include <iostream>
#include <fstream>
#include <sstream>
#include "file_header.h"
#include "helper_functions.h"

/*
 This reads in the grid from the input file specified in input_filename
 It also computes other necessary items from the grid to minimize calculation repetition later in the code.

 returns grid, a vector of cells with the following entries:
 grid[cellnum].cellnumber, .cornerlocs_x, .centroid_x, .unit_normals_x([0],[1],[2],[3]), .cornerlocs_y, .centroid_y, .unit_normals_y([0],[1],[2],[3]), .edge_type, .state (.rho,.rhou,.rhov,.E)


 Input file should be in this format:

 "cellnumber","cornerlocs_x","cornerlocs_y","adjacent_cells","initial_conditions",
 0,0,0.125,0.125,0,0,0,1.3333,1.3333,0,0,2,9,1,0,0,250000,
*/

void read_grid(std::string input_filename, std::vector<cell> &grid, double gamma)  {
	std::string line;
	std::cout << "reading grid... \n";
	int sep_loc_1;
	int sep_loc_old, sep_loc_new;
	int linenumber = 1;
	int cell_num = 0;
	std::vector<double> unit_normals_x_temp, unit_normals_y_temp;

	std::ifstream input_file(input_filename);
	if (input_file.is_open()) {
		while (getline(input_file,line)) {
			if (linenumber == 1) {
				linenumber++;
				continue;
			}
			else if (!line.empty()) {
				grid.push_back(cell());

				// Read cell number
				sep_loc_1 = line.find(',');
				std::string cell_number_str = line.substr(0, sep_loc_1);
				grid[cell_num].cellnumber = (std::stoi(cell_number_str));
				

				// Read corner node x locations and compute x centroid
				sep_loc_old = sep_loc_1;
				double x_point_sum = 0;
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string node_x_temp = line.substr(sep_loc_old+1, sep_loc_new);
					grid[cell_num].cornerlocs_x.push_back(std::stod(node_x_temp));
					x_point_sum = x_point_sum + (std::stod(node_x_temp));
					sep_loc_old = sep_loc_new;
				}
				grid[cell_num].centroid_x = (x_point_sum/4);


				// Read corner node y locations and compute y centroid	
				double y_point_sum = 0;	
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string node_y_temp = line.substr(sep_loc_old+1, sep_loc_new);	
					grid[cell_num].cornerlocs_y.push_back(std::stod(node_y_temp));
					y_point_sum = y_point_sum + (std::stod(node_y_temp));
					sep_loc_old = sep_loc_new;
				}
				grid[cell_num].centroid_y = (y_point_sum/4);

				// Compute cell area
				grid[cell_num].area = compute_cell_area(grid[cell_num]);

				// Compute unit_normals_x and unit_normals_y now that we have cell position
				compute_outward_unit_normal(grid[cell_num], unit_normals_x_temp, unit_normals_y_temp);
				grid[cell_num].unit_normals_x = unit_normals_x_temp;
				grid[cell_num].unit_normals_y = unit_normals_y_temp;


				// Read adjacent cells				
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string adjacentcell_temp = line.substr(sep_loc_old+1, sep_loc_new);	
					grid[cell_num].adjacent_cells.push_back(std::stod(adjacentcell_temp));
					sep_loc_old = sep_loc_new;
				}

				// Compute cell edge lengths and assign
				compute_cell_edge_length(grid[cell_num], grid[cell_num].edge_lengths);

				// Read edge_cell boolean
				sep_loc_new = line.find(',',sep_loc_old+1);
				std::string edge_type_temp = line.substr(sep_loc_old+1, sep_loc_new);		
				grid[cell_num].edge_type = (std::stoi(edge_type_temp));
				sep_loc_old = sep_loc_new;

				// Read initial conditions				
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string initialcondition_temp = line.substr(sep_loc_old+1, sep_loc_new);		
					if (iX == 1) {
						grid[cell_num].state.rho = (std::stod(initialcondition_temp));
					} else if (iX == 2) {
						grid[cell_num].state.rhou = (std::stod(initialcondition_temp));
					} else if (iX == 3) {
						grid[cell_num].state.rhov = (std::stod(initialcondition_temp));
					} else if (iX == 4) {
						grid[cell_num].state.rhoE = (std::stod(initialcondition_temp));
					}
					sep_loc_old = sep_loc_new;
				}
			//	grid[cell_num].state.pressure = compute_pressure_2D(grid[cell_num].state, gamma);
			}
			else {
			// We have reached the end of the input file - time to compute adjacent cell indices!		
				for(unsigned int cell_num = 0; cell_num <= grid.size()-1; cell_num++) {

					grid[cell_num].adjacent_cells_gridpos.push_back(find_cellposition(grid, grid[cell_num].adjacent_cells[0])); // left
					grid[cell_num].adjacent_cells_gridpos.push_back(find_cellposition(grid, grid[cell_num].adjacent_cells[1])); // bottom
					grid[cell_num].adjacent_cells_gridpos.push_back(find_cellposition(grid, grid[cell_num].adjacent_cells[2])); // right
					grid[cell_num].adjacent_cells_gridpos.push_back(find_cellposition(grid, grid[cell_num].adjacent_cells[3])); // top
					
					if (grid[cell_num].edge_type == 0) {
						compute_cell_distances(grid, grid[cell_num].cell_distance, cell_num);
					}
				}
				break;
			}
			linenumber++;
			cell_num++;
		}
	} else {
		std::cout << "Input file not open... file does not exist." << '\n';
	}
}
