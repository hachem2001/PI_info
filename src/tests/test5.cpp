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

	//std::ifstream myfile("../../../scenarios/scenario6v2.txt");
	graph g = empty_graph();
	add_vertices(g, 7);
	add_edge(g, 0, 1, 9);
	add_edge(g, 0, 5, 11);
	add_edge(g, 0, 6, 10);
	add_edge(g, 4, 5, 10);
	add_edge(g, 4, 6, 12);
	add_edge(g, 2, 6, 8);
	add_edge(g, 3, 6, 11);
	add_edge(g, 1, 2, 9);

	std::cout << g << std::endl;
	
	std::map<int, std::map<int, double>> min_matrix = min_distance(g);
	for (std::map<int, std::map<int, double>>::iterator v=min_matrix.begin(); v!=min_matrix.end(); v++){
		std::cout << v->first << ":";
		for (std::map<int, double>::iterator u=v->second.begin(); u!=v->second.end(); u++){
			std::cout << "(" << u->first << "," << u->second << ")" << std::endl;
		}
		std::cout << std::endl;
	}
	return 0;
}
