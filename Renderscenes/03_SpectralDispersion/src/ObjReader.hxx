#ifndef OBJ_READER_HXX_
#define OBJ_READER_HXX_

#include "math.hxx"
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

int readObjFile(const char* file_path, std::vector<Vec3f>& vertices, std::vector<int>& indices) {
	
	std::ifstream obj_file(file_path);
	std::string line;
	std::string line_type;
	if (obj_file.fail()) {
		return -4;
	}

	float x, y, z;
	int ix, iy, iz;

	//read obj vertices from file
	while (std::getline(obj_file, line)) {
		std::istringstream iss(line);
		iss >> line_type;

		if (line_type.compare("v") == 0) {
			if (!(iss >> x >> y >> z)) {
				//error with file
				return -1;
			}
			//add vertices to vector
			vertices.push_back(Vec3f(x, y, z));
		}
		else if (line_type.compare("f") == 0) {
			char slash;
			int tmp;
			if (!(iss >> ix >> slash >> tmp >> slash >> tmp >> iy >> slash >> tmp >> slash >> tmp >> iz)) {
				//error with file
				return -3;
			}
			//add vertices to vector
			indices.push_back(ix-1);
			indices.push_back(iy-1);
			indices.push_back(iz-1);
		}	
	}

	return 0;
}

#endif
