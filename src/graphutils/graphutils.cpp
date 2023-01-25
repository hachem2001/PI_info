#include "graphutils.hpp"
#include "assert.h"

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
      os << it->first;
      if (it->second.first.terminal) {
        os << "*"; // Terminals have a star next to them
      }
      os << ":" << adj << std::endl;
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
  
  graph fstream_graph(std::ifstream& stream)
  {
    graph a;
    std::string k_str; stream >> k_str; // Reads k and n
    std::string n_str; stream >> n_str;
    std::string str; // Will be used to read cells consecutively.

    int k = static_cast<int>(stod(k_str)); // Convert k and n to ints.
    int n = static_cast<int>(stod(n_str)); // We must convert to double first then to int (because of the scientific notation)

    //std::cout << "Number of terminals : " << k << std::endl;
    //std::cout << "Number of nodes : " << n << std::endl;

    add_vertices(a, n); // Adds all the vertices at once.
    for (int i=0; i<k; i++){
      set_terminal(a, i, true);
    }
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
        stream >> str;
        if (i>j) { // We don't want to add the edges twice. The matrix is symmetric. However we do need to keep reading.
          continue;
        }
        double weight = stod(str); // Read the weight in position i,j in the adjacency matrix.
        if (weight>0){ // I'm not sure if the string 0.0000 will be converted to a perfect 0 float, and that this comparison will work... Let's hope.
          add_edge(a, i, j, weight);
        }
      }
    }

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
  
  int add_vertex(graph& g, int index)
  {
    if (g.info.count(index) == 0) {
      if (index >= g._number_of_created_vertices) {
        g._number_of_created_vertices = index+1;
        // Make sure any newly added vertices with the other functions have a unique index using this trick.
      }
      label l;
      std::list<cell> empty_adjacency_list;
      g.info.emplace(index, std::pair<label, std::list<cell>>(l, empty_adjacency_list));
    }
    return index;
  }

  int add_vertices(graph& g, int n)
  {
    // Add n different vertices (and return the index of the first one).
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

      std::list<cell> &corresponding_neighs = g.info.at(a.vertex).second;
      //std::cout << "Passed by here " << corresponding_cell->vertex;
      corresponding_neighs.erase(corresponding_cell);
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

    if (!(g.info.count(vertex1)==0 || g.info.count(vertex2)==0)) {
      std::list<cell> &adj1 = g.info.at(vertex1).second;

      for (std::list<cell>::iterator it = adj1.begin(); it!=adj1.end(); it++) {
        if (it->vertex == vertex2) {
          return std::pair<bool, std::list<cell>::iterator>(true, it);
        }
      }
    }
    return std::pair<bool, std::list<cell>::iterator>(false, std::list<cell>::iterator());
  }
  
  std::pair<bool, double> get_weight(graph& g, int vertex1, int vertex2)
  {
    std::pair<bool, std::list<cell>::iterator> s = get_edge(g, vertex1, vertex2);
    if (s.first == true) {
      return std::pair<bool, double>(true, s.second->weight);
    }
    return std::pair<bool, double>(false, -1);
  }
  
  void remove_edge(graph& g, int vertex1, int vertex2)
  {
    std::list<cell> &adj1 = g.info.at(vertex1).second;
    std::list<cell> &adj2 = g.info.at(vertex2).second;

    
    for (std::list<cell>::iterator it = adj1.begin(); it!=adj1.end(); it++) {
      if (it->vertex == vertex2) { // First first cell leading to vertex2 from vertex1
        adj2.erase(it->_p); // Erase the corresponding link from vertex 2's adjacency list
        adj1.erase(it); // Erase the link from vertex 1's adjacency list
        return; // Or break;
      }
    }
  }
  
  void add_link(graph& g, int vertex1, int vertex2, double weight)
  {
    add_vertex(g, vertex1); // If the vertex already exists : it won't do anything 
    add_vertex(g, vertex2); // Same thing.

    add_edge(g, vertex1, vertex2, weight);
  }
  
  void set_terminal(graph &g, int vertex, bool set=1)
  {
    std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.find(vertex);
    assert(it != g.info.end()); // Asserts the vertex exists. Throws an error otherwise
    
    it->second.first.terminal = set;
  }
  
  int remove_leafs(graph& g)
  {
    int count = 0; // Counts the number of removed leafs.
    std::map<int, std::pair<label, std::list<cell>>>::iterator it=g.info.begin();
    while (it!=g.info.end()) {
      if (!it->second.first.terminal && it->second.second.size() <= 1) { // If not terminal and length less than 1
        disconnect_vertex(g, it->first); // Disconnect instead of delete, so we can keep iterating.
        it = g.info.erase(it); // Should increment the iterator after erasing it.
        count++;
      } else {
        it++;
      }
    }
    return count;
  }
  
  std::pair<int, int> optimize_degree_2(graph& g)
  {
    int number_of_removed_edges = 0, number_of_removed_vertices = 0;
    std::map<int, std::pair<label, std::list<cell>>>::iterator it=g.info.begin();
    while (it!=g.info.end()) {
      if (!it->second.first.terminal && it->second.second.size() == 2) { // If not terminal and length less than 1
        std::list<cell> &adj = it->second.second;
        std::list<cell>::iterator it1 = adj.begin(),  it2 = ++adj.begin();
        double weight1 = it1->weight, weight2 = it2->weight;
        std::pair<bool, double> middle_weight = get_weight(g, it1->vertex, it2->vertex);
        if (middle_weight.first == true) {
          if (weight1 + weight2 > middle_weight.second) {
            disconnect_vertex(g, it->first); // Edge
            it = g.info.erase(it);
            number_of_removed_vertices += 1;
            number_of_removed_edges += 2;
          } else {
            // Remove mediating edge
            remove_edge(g, it1->vertex, it2->vertex);
            number_of_removed_edges += 1;
            it++;
          }
        } else {
          it++;
        }
      } else {
        it++;
      }
    }
    return std::pair<int, int>(number_of_removed_vertices, number_of_removed_edges);
  }



}