#include <iostream>
#include "graphutils.hpp"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	std::cout << "PI_Project" << std::endl;

	graph g = empty_graph();
	add_vertices(g,6);
	add_edge(g, 0, 1, 2.4);
	add_edge(g, 0, 2, 2.3);
	add_edge(g, 4, 1, 2.1);
	add_edge(g, 1, 3, 1.1);
	add_edge(g, 3, 2, 1.2);

	std::cout << g << std::endl;
	return 0;
}
