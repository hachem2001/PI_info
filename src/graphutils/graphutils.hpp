#pragma once
// pragma once should work with almost all modern compilers.
#include <iostream>
#include <vector>
#include <list>

namespace Graphutils {
	struct graph {
		std::vector<std::list<int>> adjs;
	};
	std::ostream &operator<<(std::ostream &os, graph const &m);
}

