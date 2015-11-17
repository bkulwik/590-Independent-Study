#include "file_header.h"

// Computes 2D pressure from a TDstate
double compute_pressure_2D(TDstate state, double gamma);

// Computes 1D pressure from a ODstate
double compute_pressure_1D(ODstate state, double gamma);

// Computes the max value in a vector of doubles
double vectormax(std::vector<double>& locations, bool absolute_value);

// Computes the min value in a vector of doubles
double vectormin(std::vector<double>& locations, bool absolute_value);

// Computes the min delta_(x or y) value in grid
//double gridmin(std::vector<cell> grid, char direction);

// Computes and returns the cell numbers of all non-edge cells in the grid
std::vector<int> find_interior_cells(const std::vector<cell>& grid);

// Compute dt based on current wavespeeds and the cell hydraulic diameter and CFL number
void timestep_calculator(double& dt, double gamma, const std::vector<cell>& grid, const std::vector<TDstate>& Up1, double CFL, const std::vector<double>& cell_diameter, const std::vector<double>& cell_perimeter);

// Outputs the global state U to the terminal
void outputU(std::vector<TDstate>& U);

// Outputs a ODstate to the terminal
void outputODstate(ODstate& state);

// Outputs a TDstate to the terminal
void outputTDstate(TDstate& state);

// Write grid to grid_file.txt to see how it is being read in and ensure reading is correct
void output_grid(const std::vector<cell>& grid);

// Finds at what vector index number (cellposition) the desired cell number is located
int find_cellposition(const std::vector<cell>& grid, int& cellnumber_desired);

// Computes the outward unit normal vectors for the given "current_cell" and returns them in two vectors of doubles, [0]-[3] corresponding to left, bottom, right, top
void compute_outward_unit_normal(const cell& current_cell, std::vector<double>& x_unit_normals, std::vector<double>& y_unit_normals);

// Computes area of a cell given the cell's cornerlocs information contained in the cell struct
double compute_cell_area(cell& current_cell);

// Compute and return cell edge lengths in a directional_quantity
void compute_cell_distances(std::vector<cell>& grid, directional_quantity& cell_distance, unsigned int& cellposition);

// Compute and return (in cell_edge_lengths) the edge lengths of this cell. Used in read_grid.cpp
void compute_cell_edge_length(cell& current_cell, directional_quantity& cell_edge_lengths);
//==============================================================
