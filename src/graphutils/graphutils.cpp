#include "graphutils.hpp"

std::ostream &operator<<(std::ostream &os, std::list<int> adjs) {
  if (adjs.size() == 0) {
    os << "[]";
  } else if (adjs.size() == 1) {
    os << "[" << adjs.front() << "]";
  } else {
    auto before_last = --adjs.end();
    os << "[";
    for (auto it = adjs.begin(); it != before_last; it++) {
      int neighbour = *it;
      os << neighbour << ",";
    }
    os << adjs.back() << "]";
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, Graphutils::graph const &m) {

  for (int vert = 0; vert < m.adjs.size() - 1; vert++) {
    std::list<int> adj = m.adjs[vert];
    os << vert << ":" << adj << std::endl;
  }
  os << m.adjs.size() - 1 << ":" << m.adjs[m.adjs.size() - 1];
  return os;
}
