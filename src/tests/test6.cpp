#include <iostream>
#include <fstream>
#include <string>

#include "graphutils.hpp"
#include "time.h"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	clock_t start_c, end_c;
	std::cout << "Test 4 PI_Project" << std::endl;
	std::cout << "Using scenario 6v2" << std::endl;
	// Attempt to read a file and output the corresponding graph

	std::ifstream myfile("../../../scenarios/scenario5.txt");
	graph g = fstream_graph(myfile);
	/*graph g = empty_graph(); //fstream_graph(myfile);
	add_vertices(g, 7);
	add_edge(g, 0, 1, 9);
	add_edge(g, 0, 5, 11);
	add_edge(g, 0, 6, 10);
	add_edge(g, 4, 5, 10);
	add_edge(g, 4, 6, 12);
	add_edge(g, 2, 6, 8);
	add_edge(g, 3, 6, 11);
	add_edge(g, 1, 2, 9);
	add_edge(g, 0, 4, 60);
	add_edge(g, 5, 6, 60);
	set_terminal(g, 1, true);
	set_terminal(g, 4, true);
	*/

	std::cout << g << std::endl;
	std::cout << "Preprocessing on g" << std::endl;
	start_c = clock();
	int vrtx_rmvd = remove_leafs(g); // It's equal to the number of edges AND vertices removed when getting rid of leafs
	end_c = clock();

	printf("It took %.2f µs\n", 1000000*double(end_c - start_c)/CLOCKS_PER_SEC);
	int edgs_rmvd = 0;
	start_c = clock();
	std::pair<int, int> _temp = optimize_degree_2(g); vrtx_rmvd+=_temp.first;  edgs_rmvd+=_temp.second;
	end_c = clock();


	printf("It took %.2f µs\n", 1000000*double(end_c - start_c)/CLOCKS_PER_SEC);
	std::cout << "(Removed " << vrtx_rmvd << " vertices and " << vrtx_rmvd+edgs_rmvd << " edges)" << std::endl;

	std::cout << "Creating min_distance matrix" << std::endl;
	start_c = clock();
	std::map<int, std::map<int, double>> min_matrix = min_distance(g);
	end_c = clock();
	printf("It took %.2f ms\n", 1000*double(end_c - start_c)/CLOCKS_PER_SEC);
	std::cout << "The minimum distance matrix :" << std::endl;
	//std::cout << min_matrix << std::endl;

	std::cout << "Removing edges longer than corresponding min_distance" << std::endl;
	start_c = clock();
	edgs_rmvd = optimize_edge_with_min_distance_matrix(g, min_matrix);
	end_c = clock();
	std::cout << "(Removed " << edgs_rmvd << " edges)" << std::endl;
	printf("It took %.2f ms\n", 1000*double(end_c - start_c)/CLOCKS_PER_SEC);
	//std::cout << g << std::endl;

	std::cout << "Transforming minimum matrix to graph" << std::endl;
	start_c = clock();
	graph min_dist_graph = min_distance_graph(g, min_matrix);
	end_c = clock();
	printf("It took %.2f ms\n", 1000*double(end_c - start_c)/CLOCKS_PER_SEC);
	//std::cout << min_dist_graph << std::endl;
	

	std::cout << "Minimum Spanning Tree for D" << std::endl;
	start_c = clock();
	graph mst_d = primm_mst(min_dist_graph);
	end_c = clock();
	printf("It took %.2f ms\n", 1000*double(end_c - start_c)/CLOCKS_PER_SEC);

	std::cout << "Calculating the minimum steiner tree in D" << std::endl;
	start_c = clock();
	std::pair<graph,graph> result = enumeration_steiner_tree(g, min_dist_graph);
	end_c = clock();
	printf("It took %.2f ms\n", 1000*double(end_c - start_c)/CLOCKS_PER_SEC);

	std::cout << result.second << std::endl;
	return 0;
}
