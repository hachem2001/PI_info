#include "graphutils.hpp"

// For now this print implementation does not show the label. This can be easily modified later.


namespace Graphutils
{
  std::ostream &operator<<(std::ostream &os, std::list<cell> neighs) {
    if (neighs.size() == 0) {
      os << "[]";
    } /*else if (neighs.size() == 1) {
      os << "[" << neighs.front() << "]";
    } */else {
      auto before_last = --neighs.end();
      os << "[";
      for (auto it = neighs.begin(); it != before_last; it++) {
        cell neighbour = *it;
        os << "(" << neighbour.vertex << ", " << neighbour.weight << ")" << ";";
      }
      os << "(" << neighs.back().vertex << ", " << neighs.back().weight << ")" << "]";
    }

    return os;
  }

  std::ostream &operator<<(std::ostream &os, graph const &m) {

    for (std::map<int, std::pair<label, std::list<cell>>>::const_iterator it=m.info.begin(); it!=m.info.end(); ++it){
      std::list<cell> adj = it->second.second;
      os << it->first << ":" << adj << std::endl;
    }

    return os;
  }

  graph empty_graph()
  {
    graph a;
    // Nothing special here really. Just creating a structure and sending it.
    // With a struct, we worry less about properly freeing memory upon deletion.
    return a;
  }
  
  int add_vertex(graph& g)
  {
    g._number_of_created_vertices++;
    g.number_of_vertices++;
    label l; l.index=g._number_of_created_vertices;
    std::list<cell> empty_adjacency_list;
    g.info.emplace(g.number_of_vertices-1, std::pair<label, std::list<cell>>(l, empty_adjacency_list)); //(std::pair<label, std::list<cell>>(l, empty_adjacency_list));
    return g.number_of_vertices;
  }
  
  int add_vertices(graph& g, int n)
  {
    int first_index = add_vertex(g);
    for (int i=0; i<n-1; i++){
      add_vertex(g);
    }
    return first_index;
  }
  
  void disconnect_vertex(graph& g, int vertex)
  {
    std::list<cell> &neighs = g.info.at(vertex).second; // List of all the neigbours to vertex
    while (neighs.size() > 0) {
      // Remove each connection from both the vertex's list and the list of the element in question
      cell a = neighs.front();
      std::list<cell>::iterator &corresponding_cell = a._p;

      //std::list<cell> &corresponding_neighs = g.info.at(a.vertex).second;
      std::cout << "Passed by here " << corresponding_cell->vertex;
      //corresponding_neighs.erase(corresponding_cell);
      neighs.pop_front();
    }
  }

  void remove_vertex(graph& g, int vertex)
  {
    disconnect_vertex(g, vertex);
    g.info.erase(vertex);
    g.number_of_vertices--;
  }
  
  void add_edge(graph& g, int vertex1, int vertex2, double weight)
  {
    // This function is optimised as in to insert the new cells and have them ordered by weight.
    // Since we create the graph once and only seem to create few points in the future : this seems like a better approach.
    cell c1, c2; // The two cells to add.
    c1.vertex = vertex2; c2.vertex = vertex1;
    c1.weight = c2.weight = weight;
    
    // Insert both.
    std::list<cell> &adj1 = g.info.at(vertex1).second, &adj2 = g.info.at(vertex2).second;

    std::list<cell>::iterator it1, it2;
    if (adj1.size() == 0){
      adj1.push_front(c1);
      it1 = adj1.begin();
    } else {
      it1 = adj1.begin();
      while (it1!=adj1.end() && it1->weight < weight){
        it1++;
      }
      it1 = adj1.insert(it1, c1);
    }
    c2._p = it1;

    if (adj2.size() == 0){
      adj2.emplace_front(c2);
      it2 = adj2.begin();
    } else {
      it2 = adj2.begin();
      while (it2!=adj2.end() && it2->weight < weight) {
        it2++;
      }
      it2 = adj2.insert(it2, c2);
    }
    it1->_p = it2;
  }
  
  std::pair<bool, std::list<cell>::iterator> get_edge(graph& g, int vertex1, int vertex2)
  {
    std::list<cell> &adj1 = g.info.at(vertex1).second;

    for (std::list<cell>::iterator it = adj1.begin(); it!=adj1.end(); it++) {
      if (it->vertex == vertex2) {
        return std::pair<bool, std::list<cell>::iterator>(true, it);
      }
    }
    return std::pair<bool, std::list<cell>::iterator>(false, adj1.begin());
  }
  
  void remove_edge(graph& g, int vertex1, int vertex2)
  {
    std::list<cell> &adj1 = g.info.at(vertex1).second;
    std::list<cell> &adj2 = g.info.at(vertex2).second;

    for (std::list<cell>::const_iterator it = adj1.begin(); it!=adj1.end(); it++) {
      if (it->vertex == vertex2) {
        
        return;
      }
    }
  }
}