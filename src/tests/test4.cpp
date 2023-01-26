#include <iostream>
#include <fstream>
#include <string>

#include "graphutils.hpp"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	std::cout << "Test 4 PI_Project" << std::endl;
	// Attempt to read a file and output the corresponding graph

	std::ifstream myfile("../../../scenarios/scenario6v2.txt");
	graph g = fstream_graph(myfile);

	std::cout << g << std::endl;

	std::cout << "Preprocessing on g" << std::endl;
	int vrtx_rmvd = remove_leafs(g); // It's equal to the number of edges AND vertices removed when getting rid of leafs
	int edgs_rmvd = 0;
	std::pair<int, int> _temp = optimize_degree_2(g); vrtx_rmvd+=_temp.first;  edgs_rmvd+=_temp.second;
	std::cout << g;
	std::cout << "(Removed " << vrtx_rmvd << " vertices and " << vrtx_rmvd+edgs_rmvd << " edges)" << std::endl;
	
	
	return 0;
}
