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
	// Attempt to read a file and output the corresponding graph

	std::ifstream myfile("../../../scenarios/scenario6v2.txt");
	graph g = fstream_graph(myfile);
	std::cout << "Working with scenario6v2" << std::endl;
	//std::cout << g << std::endl;

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

	//std::cout << g;
	std::cout << "(Removed " << vrtx_rmvd << " vertices and " << vrtx_rmvd+edgs_rmvd << " edges)" << std::endl;
	
	std::cout << "Creating min_distance matrix";
	start_c = clock();
	std::map<int, std::map<int, double>> min_matrix = min_distance(g);
	end_c = clock();
	printf("It took %.2f ms\n", 1000*double(end_c - start_c)/CLOCKS_PER_SEC);

	/*
	for (std::map<int, std::map<int, double>>::iterator v=min_matrix.begin(); v!=min_matrix.end(); v++){
		std::cout << v->first << ":";
		for (std::map<int, double>::iterator u=v->second.begin(); u!=v->second.end(); u++){
			std::cout << "(" << u->first << "," << u->second << ")";
		}
		std::cout << std::endl;
	}*/

	std::cout << "Finished calculating" << std::endl;
	

	return 0;
}
