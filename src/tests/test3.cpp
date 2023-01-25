#include <iostream>
#include "graphutils.hpp"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	std::cout << "PI_Project" << std::endl;

	graph g = empty_graph();
	add_vertices(g, 9); // Let's make a turning star formation.

	add_edge(g, 1, 2, 2.3);
	add_edge(g, 2, 3, 2.3);
	add_edge(g, 3, 1, 2.3);

	add_edge(g, 4, 5, 2.3);
	add_edge(g, 5, 6, 2.3);
	add_edge(g, 6, 4, 2.3);

	add_edge(g, 0, 1, 2.3);
	add_edge(g, 0, 2, 2.3);

	add_edge(g, 0, 4, 2.3);
	add_edge(g, 0, 5, 2.3);

	add_edge(g, 3, 7, 2.3);
	add_edge(g, 5, 8, 2.3);
	std::cout << "Initial Graph" << std::endl;
	std::cout << g << std::endl;
	
	std::cout << "Set 0 as a terminal" << std::endl;
	set_terminal(g, 0, true);
	std::cout << g << std::endl;

	std::cout << "Remove leafs." << std::endl;
	remove_leafs(g); // Words as intended.
	std::cout << g << std::endl;

	std::cout << "Traverse 2-link edges" << std::endl;
	optimize_degree_2(g);
	std::cout << g << std::endl;
	return 0;
}
