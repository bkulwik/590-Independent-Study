// This is used to read in all solver parameters in the solver_parms.txt file
// It may need to be updated with other parameters as we go, and as other parameters need to be taken into account in this code.


#ifndef FILE_HEADER_H
#define FILE_HEADER_H
#include <vector>
#include <string>


struct ODstate { // ODstate for 1-dimensional state
	double rho;
	double rhou;
	double rhoE;
	double pressure;
};

struct TDstate { // TDstate for 2-dimensional state
	double rho;
	double rhou;
	double rhov;
	double rhoE;
};

struct directional_quantity {
	double left;
	double bottom;
	double right;
	double top;
};
 
struct cell {
	int cellnumber;
	std::vector<double> cornerlocs_x; //bottom left, bottom right, top right, top left
	double centroid_x;
	std::vector<double> unit_normals_x; //left, bottom, right, top - outward unit normals
	std::vector<double> cornerlocs_y; //bottom left, bottom right, top right, top left
	double centroid_y;
	std::vector<double> unit_normals_y; //left, bottom, right, top - outward unit normals
	std::vector<int> adjacent_cells; //left, bottom, right, top
	std::vector<int> adjacent_cells_gridpos; //left, bottom, right, top
	int edge_type;
	TDstate state;
	double area;
	directional_quantity cell_distance; //distance between left, bottom, right and top cell's centroids and current cell's centroids
	directional_quantity edge_lengths;
};


//This is used to read the parameters from the parameter file
std::vector<std::string> read_parameters(std::string filename);

//This is used to read the input node locations and cell connectivity from the given input_filename
void read_grid(std::string input_filename, std::vector<cell> &grid, double gamma);

// This is the 1D exact riemann solver
ODstate Exact_Riemann_Solver(ODstate left, ODstate right, double left_parallelvel, double right_parallelvel, double thresh, double gamma, bool debug);

//This is used to write solver outputs to a .txt file to be read into and plotted by MATLAB
void write_to_file(const std::vector<cell>& grid, const std::vector<TDstate>& U, std::string& filename);

#endif
//==============================================================
