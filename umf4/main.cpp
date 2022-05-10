#include "Solver.hpp"

int main() {
	using ::std::filesystem::path;

	path dir_to_file = "files/test4";
	Solver s(dir_to_file);

	return 0;
}