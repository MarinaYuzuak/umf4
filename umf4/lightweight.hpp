#ifndef _LIGHTWEIGHT_HPP
#define _LIGHTWEIGHT_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

void printErrorFile(std::string file_name) {
	std::cerr << "File " << file_name << " is not open." << std::endl;
	std::exit(1);
}

template <typename T>
void readVector(std::vector<T>& v, std::ifstream& fin) {
	for (size_t i = 0; i < v.size(); i++)
		fin >> v[i];
}

void printVector(std::vector<double>& v) {
	for (size_t i = 0; i < v.size(); i++)
		std::cout << std::setprecision(15) << v[i] << std::endl;
}

#endif // !_LIGHTWEIGHT_HPP
