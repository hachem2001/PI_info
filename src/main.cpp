#include <iostream>
#include "graphutils.hpp"

using namespace Graphutils;
// Main
int main(int argc, char *argv[])
{
	std::cout << "PI_Project" << std::endl;

	std::set<int> a = {1,2,3,5,6,7,8};
	std::vector<std::set<int>> b = Setutils::get_subsets(a, 3);

	for (auto v: b) {
		for (auto m: v) {
			std::cout << m << "-";
		}
		std::cout << std::endl;
	}
	return 0;
}
