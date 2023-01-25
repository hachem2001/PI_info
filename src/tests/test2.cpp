#include <iostream>
#include "graphutils.hpp"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	std::cout << "PI_Project" << std::endl;

	graph g = empty_graph();
	add_vertices(g, 6); // Let's make a star formation : one terminal only connected with leafs
	add_edge(g, 0, 1, 2.3);
	add_edge(g, 0, 2, 2.3);
	add_edge(g, 0, 3, 2.3);
	add_edge(g, 0, 4, 2.3);
	add_edge(g, 0, 5, 2.3);

	std::cout << "Initial Graph" << std::endl;
	std::cout << g << std::endl;
	
	std::cout << "Set 0 and 5 as a terminal" << std::endl;
	set_terminal(g, 0, true);
	set_terminal(g, 5, true);
	std::cout << g << std::endl;
	return 0;
}
