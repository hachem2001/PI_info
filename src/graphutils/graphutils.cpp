#include "graphutils.hpp"

// For now this print implementation does not show the label. This can be easily modified later.
std::ostream &operator<<(std::ostream &os, std::list<Graphutils::cell> neighs) {
  if (neighs.size() == 0) {
    os << "[]";
  } /*else if (neighs.size() == 1) {
    os << "[" << neighs.front() << "]";
  } */else {
    auto before_last = --neighs.end();
    os << "[";
    for (auto it = neighs.begin(); it != before_last; it++) {
      Graphutils::cell neighbour = *it;
      os << "(" << neighbour.vertex << "," << neighbour.weight << ")" << ";";
    }
    os << "(" << neighs.back().vertex << "," << neighs.back().weight << ")" << "]";
  }

  return os;
}

std::ostream &operator<<(std::ostream &os, Graphutils::graph const &m) {

  for (int vert = 0; vert < m.info.size() - 1; vert++) {
    std::list<Graphutils::cell> adj = m.info[vert].second;
    os << vert << ":" << adj << std::endl;
  }
  os << m.info.size() - 1 << ":" << m.info[m.info.size() - 1].second;
  return os;
}
