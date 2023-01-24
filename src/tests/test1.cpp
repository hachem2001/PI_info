#include <iostream>
#include "graphutils.hpp"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	std::cout << "PI_Project" << std::endl;

	graph g = empty_graph();
	add_vertices(g, 6);
	add_edge(g, 0, 1, 2.4);
	add_edge(g, 0, 2, 2.3);
	add_edge(g, 4, 1, 2.1);
	add_edge(g, 1, 3, 1.1);
	add_edge(g, 3, 2, 1.2);
	std::cout << "Initial Graph" << std::endl;
	std::cout << g << std::endl;
	
	std::cout << "Adding the edge (3 ~ 4,1.2)" << std::endl;
	add_edge(g, 3, 4, 1.2);
	std::cout << g << std::endl;

	std::cout << "Disconnecting 3" << std::endl;
	disconnect_vertex(g, 3);
	std::cout << g << std::endl;

	std::cout << "Removing vertex 2" << std::endl;
	remove_vertex(g, 2);
	std::cout << g << std::endl;

	std::cout << "Removing vertex 3" << std::endl;
	remove_vertex(g, 3);
	std::cout << g << std::endl;

	std::cout << "FUN!" << std::endl << std::endl;

	std::cout << "Retrieving weights 0~1, 1~3, 4~1" << std::endl;
	std::cout << get_weight(g, 0, 0).second << std::endl;
	std::cout << get_weight(g, 3, 1).second << std::endl;
	std::cout << get_weight(g, 4, 1).second << std::endl;

	std::cout << "Adding vertices+edge (2 ~ 6, 0.3)" << std::endl;
	add_link(g, 2, 6, 0.3);
	std::cout << g << std::endl;
	return 0;
}
