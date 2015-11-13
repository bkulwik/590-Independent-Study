#include "file_header.h"
#include "helper_functions.h"
#include <iostream>
#include <fstream>


void output_4doubles(const std::vector<double>& cornerlocs, std::ofstream& ofile);

void output_4ints(const std::vector<int>& adjacent_cells, std::ofstream& ofile);

void output_state(const TDstate& U, std::ofstream& ofile);

void write_to_file(const std::vector<cell>& grid, const std::vector<TDstate>& U, std::string& filename) {
	std::ofstream ofile;
	ofile.open(filename);
	
	ofile << "'cellnumber','cornerlocs_x','cornerlocs_y','adjacent_cells','edge_state','initial_conditions'" << '\n';

	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) {
		ofile << grid[cellnum].cellnumber << ",";
		output_4doubles(grid[cellnum].cornerlocs_x, ofile);
		output_4doubles(grid[cellnum].cornerlocs_y, ofile);
		output_4ints(grid[cellnum].adjacent_cells, ofile);
		ofile << grid[cellnum].edge_type << ",";
		output_state(U[cellnum], ofile);
		ofile << '\n';
	}
	ofile << '\n';
	ofile.close();
}

void output_4doubles(const std::vector<double>& to_output, std::ofstream& ofile) {
	for(unsigned int index = 0; index < 4; ++index) {
		ofile << to_output[index] << ",";
	}
}

void output_4ints(const std::vector<int>& to_output, std::ofstream& ofile) {
	for(unsigned int index = 0; index < 4; ++index) {
		ofile << to_output[index] << ",";
	}
}

void output_state(const TDstate& U, std::ofstream& ofile) {
	ofile << U.rho << "," << U.rhou << "," << U.rhov << "," << U.rhoE << ",";
}

