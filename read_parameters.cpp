
//This is used to read all parameters from the input file
#include "file_header.h"
#include <iostream>
#include <fstream>

std::vector<std::string> read_parameters(std::string filename) {
	std::vector<std::string> inputs;
	std::string line;
	int linenum = 0;

// Read input file with defined solver parameters
// First - CFL, second - tmax, third - input file name (*.bkcfd), fourth - output file name (*.txt), fifth - debug mode boolean

	std::ifstream parameter_file(filename);
	if (parameter_file.is_open()) {
		while (getline (parameter_file,line)) {
			if (line[0] != '%') {
				inputs.push_back(line);
				linenum++;
			}
		}
	parameter_file.close();
	}
	else {
		std::cout << "Cannot Open File" << '\n';
	}
return inputs;
}

