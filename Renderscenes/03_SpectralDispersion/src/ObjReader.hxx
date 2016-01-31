#ifndef OBJ_READER_HXX_
#define OBJ_READER_HXX_

#include "math.hxx"
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

int readObjFile(const char* file_path, std::vector<Vec3f>& vertices, std::vector<int>& indices) {
	
	std::ifstream obj_file(file_path);
	std::string line;
	std::string line_type;

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
			if (!(iss >> ix >> slash >> iy >> slash >> iz)) {
				//error with file
				return -2;
			}
			//add vertices to vector
			indices.push_back(ix);
			indices.push_back(iy);
			indices.push_back(iz);
		}	
	}

	return 0;
}

#endif
