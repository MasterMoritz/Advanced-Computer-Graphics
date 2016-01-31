#ifndef OBJ_READER_HXX_
#define OBJ_READER_HXX_

#include "math.hxx"
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

std::vector<Vec3f> readObjFile(const char* file_path) {
	
	std::ifstream obj_file(file_path);
	std::string line;
	std::string line_type;

	std::vector<Vec3f> vertices;
	float x, y, z;

	//read obj vertices from file
	while (std::getline(obj_file, line)) {
		std::istringstream iss(line);
		iss >> line_type;

		if (line_type.compare("v") == 0) {
			if (!(iss >> x >> y >> z)) {
				//error with file
				break;
			}
			//add vertices to vector
			vertices.push_back(Vec3f(x, y, z));
		}
	}

	for (int i = 0; i < vertices.size(); i++) {
		std::cout << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << "\n";
	}
	return vertices;
}

#endif
