#include "graphutils.hpp"
#include "assert.h"

#include <queue>
#include <cmath>
#include <algorithm> // For std::set_union.
//#include <bitset>

// For now this print implementation does not show the label. This can be easily modified later.
namespace Setutils
{
  /*
  unsigned long int _bit_next_permutation(unsigned long int v) { // Binary Hack
    unsigned long int t = (v | (v - 1)) + 1;
    unsigned long int w = t | ((((t & -t) / (v & -v)) >> 1) - 1);
    return w; // If the input is 000111, the result is 001011. Then 001101 then 001110 ... Useful bit operation.
  }*/

  void get_subsets_rec(std::vector<int>& set, int k, int idx, std::vector<int>& current, std::vector<std::vector<int>>& solution){
    if (current.size() == k) {
      solution.emplace_back(std::vector<int>(current));
    }

    for (int i=idx; i<set.size(); i++) {
      current.emplace_back(set[i]);
      get_subsets_rec(set, k, i+1, current, solution);
      current.pop_back();
    }
  }

  std::vector<std::set<int>> get_subsets(std::set<int>& set, int k)
  {
    std::vector<int> set_as_vector(set.begin(), set.end());
    std::vector<int> current;
    std::vector<std::vector<int>> result;
    get_subsets_rec(set_as_vector, k, 0, current, result);
    std::vector<std::set<int>> result_as_vect_of_sets;
    for (auto vect:result){
      std::set<int> _temp(vect.begin(), vect.end());
      result_as_vect_of_sets.emplace_back(_temp);
    }
    /*for (int element: set_as_vector) {
      std::set<int> s; s.insert(element);
      results.emplace_back(s);
      previous.emplace_back(s);
    }*/

    /*
    for (int i=0; i<k; i++)
    { // Everytime, current has at most 3 bits.
      double second_counter = 0;
      for (std::set<int> s:previous) {
        int counter = 0;
        for (int element: set_as_vector) {
          if (counter >= second_counter+i) {
            std::set<int> set_to_add = std::set<int>(s);
            set_to_add.insert(element);
            current.emplace_back(set_to_add);
          }
          counter++;
        }

        second_counter+=1.0/(i+1.0);
      }
      for (std::set<int> s:current) {
        results.push_back(s);
      }
      previous = current;
      current = std::vector<std::set<int>>();
    }
    */
    return result_as_vect_of_sets;
  }

}

namespace Graphutils
{
  std::ostream& operator<<(std::ostream& os, std::list<cell> neighs) {
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

  std::ostream& operator<<(std::ostream& os, graph const& m) {

    for (std::map<int, std::pair<label, std::list<cell>>>::const_iterator it = m.info.begin(); it != m.info.end(); ++it) {
      std::list<cell> adj = it->second.second;
      os << it->first;
      if (it->second.first.terminal) {
        os << "*"; // Terminals have a star next to them
      }
      os << ":" << adj << std::endl;
    }

    return os;
  }

  std::ostream& operator<<(std::ostream& os, std::map<int, std::map<int, double>> const& min_matrix)
  {
    for (std::map<int, std::map<int, double>>::const_iterator v = min_matrix.begin(); v != min_matrix.end(); v++) {
      os << v->first << ":";
      for (std::map<int, double>::const_iterator u = v->second.begin(); u != v->second.end(); u++) {
        os << "(" << u->first << "," << u->second << ")" << ";";
      }
      os << std::endl;
    }
    os << std::endl;
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
    for (int i = 0; i < k; i++) {
      set_terminal(a, i, true);
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        stream >> str;
        if (i > j) { // We don't want to add the edges twice. The matrix is symmetric. However we do need to keep reading.
          continue;
        }
        double weight = stod(str); // Read the weight in position i,j in the adjacency matrix.
        if (weight > 0) { // I'm not sure if the string 0.0000 will be converted to a perfect 0 float, and that this comparison will work... Let's hope.
          add_edge(a, i, j, weight);
        }
      }
    }

    return a;
  }

  int add_vertex(graph& g)
  {
    label l; l.index = g._number_of_created_vertices;
    std::list<cell> empty_adjacency_list;
    g.info.emplace(g.number_of_vertices, std::pair<label, std::list<cell>>(l, empty_adjacency_list)); //(std::pair<label, std::list<cell>>(l, empty_adjacency_list));
    g._number_of_created_vertices++;
    g.number_of_vertices++;
    return g.number_of_vertices;
  }

  int add_vertex(graph& g, int index)
  {
    if (g.info.count(index) == 0) {
      if (index >= g._number_of_created_vertices) {
        g._number_of_created_vertices = index + 1;
        // Make sure any newly added vertices with the other functions have a unique index using this trick.
      }
      label l; l.index = g._number_of_created_vertices - 1;
      std::list<cell> empty_adjacency_list;
      g.info.emplace(index, std::pair<label, std::list<cell>>(l, empty_adjacency_list));
    }
    g.number_of_vertices++;
    return index;
  }

  int add_vertices(graph& g, int n)
  {
    // Add n different vertices (and return the index of the first one).
    int first_index = add_vertex(g);
    for (int i = 0; i < n - 1; i++) {
      add_vertex(g);
    }
    return first_index;
  }

  void disconnect_vertex(graph& g, int vertex)
  {
    std::list<cell>& neighs = g.info.at(vertex).second; // List of all the neigbours to vertex
    while (neighs.size() > 0) {
      // Remove each connection from both the vertex's list and the list of the element in question
      cell a = neighs.front();
      std::list<cell>::iterator& corresponding_cell = a._p;

      std::list<cell>& corresponding_neighs = g.info.at(a.vertex).second;
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
    std::list<cell>& adj1 = g.info.at(vertex1).second, & adj2 = g.info.at(vertex2).second;

    std::list<cell>::iterator it1, it2;
    if (adj1.size() == 0) {
      adj1.push_front(c1);
      it1 = adj1.begin();
    }
    else {
      it1 = adj1.begin();
      while (it1 != adj1.end() && it1->weight < weight) {
        it1++;
      }
      it1 = adj1.insert(it1, c1);
    }
    c2._p = it1;

    if (adj2.size() == 0) {
      adj2.emplace_front(c2);
      it2 = adj2.begin();
    }
    else {
      it2 = adj2.begin();
      while (it2 != adj2.end() && it2->weight < weight) {
        it2++;
      }
      it2 = adj2.insert(it2, c2);
    }
    it1->_p = it2;
    g.number_of_edges += 1;
    g.total_cost += weight;
  }

  std::pair<bool, std::list<cell>::iterator> get_edge(graph& g, int vertex1, int vertex2)
  {

    if (!(g.info.count(vertex1) == 0 || g.info.count(vertex2) == 0)) {
      std::list<cell>& adj1 = g.info.at(vertex1).second;

      for (std::list<cell>::iterator it = adj1.begin(); it != adj1.end(); it++) {
        if (it->vertex == vertex2) {
          return std::pair<bool, std::list<cell>::iterator>(true, it);
        }
      }
    }
    return std::pair<bool, std::list<cell>::iterator>(false, std::list<cell>::iterator());
  }

  void add_edge_without_repetition(graph& g, int vertex1, int vertex2, double weight)
  {
    std::list<cell>& adj1 = g.info.at(vertex1).second, & adj2 = g.info.at(vertex2).second;
    for (std::list<cell>::iterator it=adj1.begin(); it!=adj1.end(); it++) {
      if (it->vertex == vertex2) return;
    }
    add_edge(g, vertex1, vertex2, weight);
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
    std::list<cell>& adj1 = g.info.at(vertex1).second;
    std::list<cell>& adj2 = g.info.at(vertex2).second;


    for (std::list<cell>::iterator it = adj1.begin(); it != adj1.end(); it++) {
      if (it->vertex == vertex2) { // First first cell leading to vertex2 from vertex1
        g.number_of_edges--;
        g.total_cost -= it->weight;
        adj2.erase(it->_p); // Erase the corresponding link from vertex 2's adjacency list
        adj1.erase(it); // Erase the link from vertex 1's adjacency list
        return; // Or break;
      }
    }
  }

  void remove_edge(graph& g, std::list<cell>& adj, std::list<cell>::iterator& it) { // Removes edge using cell
    std::list<cell>& adj2 = g.info.at(it->vertex).second;
    adj.erase(it);
    adj2.erase(it->_p); // Efficiency of augmented adjacency lists exploited at it's best :)
  }

  void add_link(graph& g, int vertex1, int vertex2, double weight)
  {
    add_vertex(g, vertex1); // If the vertex already exists : it won't do anything 
    add_vertex(g, vertex2); // Same thing.

    add_edge(g, vertex1, vertex2, weight);
  }

    
  void add_link_without_repetition(graph& g, int vertex1, int vertex2, double weight)
  {
    add_vertex(g, vertex1);
    add_vertex(g, vertex2);

    add_edge_without_repetition(g, vertex1, vertex2, weight);
  }

  void set_terminal(graph& g, int vertex, bool set = 1)
  {
    std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.find(vertex);
    assert(it != g.info.end()); // Asserts the vertex exists. Throws an error otherwise

    it->second.first.terminal = set;
  }

  int remove_leafs(graph& g)
  {
    int count = 0; // Counts the number of removed leafs.
    std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.begin();
    while (it != g.info.end()) {
      if (!it->second.first.terminal && it->second.second.size() <= 1) { // If not terminal and length less than 1
        disconnect_vertex(g, it->first); // Disconnect instead of delete, so we can keep iterating.
        it = g.info.erase(it); // Should increment the iterator after erasing it.
        count++;
      }
      else {
        it++;
      }
    }
    return count;
  }

  std::pair<int, int> optimize_degree_2(graph& g)
  {
    int number_of_removed_edges = 0, number_of_removed_vertices = 0;
    std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.begin();
    while (it != g.info.end()) {
      if (!it->second.first.terminal && it->second.second.size() == 2) { // If not terminal and length less than 1
        std::list<cell>& adj = it->second.second;
        std::list<cell>::iterator it1 = adj.begin(), it2 = ++adj.begin();
        double weight1 = it1->weight, weight2 = it2->weight;
        std::pair<bool, double> middle_weight = get_weight(g, it1->vertex, it2->vertex);
        if (middle_weight.first == true) {
          if (weight1 + weight2 > middle_weight.second) {
            disconnect_vertex(g, it->first); // Edge
            it = g.info.erase(it);
            number_of_removed_vertices += 1;
            number_of_removed_edges += 2;
          }
          else {
            // Remove mediating edge
            remove_edge(g, it1->vertex, it2->vertex);
            number_of_removed_edges += 1;
            it++;
          }
        }
        else {
          it++;
        }
      }
      else {
        it++;
      }
    }
    return std::pair<int, int>(number_of_removed_vertices, number_of_removed_edges);
  }

  class _comparison_for_min_distance_to_source {
    bool reverse;
  public:
    _comparison_for_min_distance_to_source(const bool& revparam = false)
    {
      reverse = revparam;
    }
    bool operator() (const std::pair<int, double>& lhs, const std::pair<int, double>& rhs) const
    {
      if (!reverse) {
        if (lhs.second < 0) { return true; }
        if (rhs.second < 0) { return false; }
        return lhs.second > rhs.second;
      }
      else {
        if (lhs.second < 0) { return false; }
        if (rhs.second < 0) { return true; }
        return lhs.second < rhs.second;
      }
    }
  };

  std::map<int, double> min_distance_to_source(graph& g, int vertex)
    // Basically Dijkstra's algorithm
  {
    std::map<int, double> result; // -1 for infinity.
    result.emplace(vertex, 0);

    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, _comparison_for_min_distance_to_source> queue;

    for (auto v_info : g.info) {
      if (v_info.first != vertex) { // Go through all the vertex
        result.emplace(v_info.first, -1);
      }
      queue.push(std::pair<int, double>(v_info.first, result.at(v_info.first)));
    }

    while (!queue.empty()) {
      std::pair<int, double> u = queue.top();

      std::pair<label, std::list<cell>>& u_info = g.info.at(u.first);
      std::list<cell>& adj = u_info.second;

      for (std::list<cell>::iterator it = adj.begin(); it != adj.end(); it++) {
        double d = result.at(u.first) + it->weight;
        double dist_neigh = result.at(it->vertex);
        if (dist_neigh < 0 || d < dist_neigh) {
          if (result.count(it->vertex) == 0) {
            result.emplace(it->vertex, d);
          }
          else {
            result[it->vertex] = d;
          }
          queue.push(std::pair<int, double>(it->vertex, d));
        }
      }
      queue.pop();
    }

    return result;
  }

  std::map<int, std::map<int, double>> min_distance(graph& g)
  {
    std::map<int, std::map<int, double>> min_distance_matrix;

    for (auto v_info : g.info) {
      min_distance_matrix.emplace(v_info.first, min_distance_to_source(g, v_info.first));
    }

    return min_distance_matrix;
  }

  graph min_distance_graph(graph& g, std::map<int, std::map<int, double>> min_distance_matrix)
  {
    graph d = empty_graph();
    typedef std::map<int, std::map<int, double>>::iterator __map_iterator; // For less writing

    // Add all the vertices and eges
    for (__map_iterator it = min_distance_matrix.begin(); it != min_distance_matrix.end(); it++) {
      add_vertex(d, it->first);
      if (g.info.at(it->first).first.terminal) {
        set_terminal(d, it->first, true);
      }
      for (std::map<int, double>::iterator jt = it->second.begin(); jt != it->second.end(); jt++) {
        if (jt->first > it->first) { // In order not to go through the same edge twice.
          add_vertex(d, jt->first);
          add_edge(d, it->first, jt->first, jt->second); // Add the edge.
        }
      }
    }
    return d; // Return the graph D.
  }

  class _comparison_for_primm_algorithm {
    bool reverse;
  public:
    _comparison_for_primm_algorithm(const bool& revparam = false)
    {
      reverse = revparam;
    }
    bool operator() (const std::pair<int, double>& lhs, const std::pair<int, double>& rhs) const
    {
      if (!reverse) {
        if (lhs.second < 0) { return true; }
        if (rhs.second < 0) { return false; }
        return lhs.second > rhs.second;
      }
      else {
        if (lhs.second < 0) { return false; }
        if (rhs.second < 0) { return true; }
        return lhs.second < rhs.second;
      }
    }
  };

  graph primm_mst(graph& g)
  {

    graph mst;
    std::map<int, std::pair<int, double>> result;
    std::map<int, bool> visited; // Keeps track of all visited nodes.
    std::map<int, double> value; // Minimum weight for each source.
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, _comparison_for_primm_algorithm> queue;

    int first_vertex = g.info.begin()->first;
    add_vertex(mst, first_vertex);
    queue.push(std::pair<int, double>(first_vertex, -1));
    typedef std::map<int, std::pair<label, std::list<cell>>>::iterator __ginfoiterator;
    for (__ginfoiterator v_info = g.info.begin(); v_info != g.info.end(); v_info++) {
      value.emplace(v_info->first, -1); // As for Dijkstra's algorithm, I use -1 to represent infinities,
    }

    while (!queue.empty()) {
      std::pair<int, double> u = queue.top();queue.pop(); // Get the vertex connected by the smallest weight, and remove from queue
      visited[u.first] = true; // Mark as visited

      std::pair<label, std::list<cell>>& u_info = g.info.at(u.first);
      std::list<cell>& adj = u_info.second;

      for (std::list<cell>::iterator it = adj.begin(); it != adj.end(); it++) {
        int vertex = it->vertex;
        double weight = it->weight;

        if (!visited[vertex] && (value[vertex] < 0 || value[vertex] > weight)) { // If not visited and has lower weight : update the queue.
          value[vertex] = weight;
          result[vertex] = std::pair<int, double>(u.first, weight);
          queue.push(std::pair<int, double>(vertex, weight));
        }
      }
    }
    for (std::map<int, std::pair<int, double>>::iterator it = result.begin(); it != result.end(); it++) { // Create the minimum spanning tree here by going through the nodes in order.
      add_link(mst, it->first, it->second.first, it->second.second);
    }

    for (std::map<int, std::pair<label, std::list<cell>>>::iterator it = mst.info.begin(); it!=mst.info.end(); it++){
      if (g.info[it->first].first.terminal) {
        set_terminal(mst, it->first, true);
      }
    }
    return mst;
  }

  int optimize_edge_with_min_distance_matrix(graph& g, std::map<int, std::map<int, double>> min_dist_matrix) // This is for the question 6
  {
    int number_of_removed_edges = 0;
    std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.begin();
    while (it != g.info.end()) {
      std::list<cell>& adj = it->second.second;
      std::list<cell>::iterator it_v = adj.begin();
      int vertex1 = it->second.first.index, vertex2; // The indices of the vertices
      while (it_v != adj.end()) {
        vertex2 = it_v->vertex;
        if (it_v->weight > min_dist_matrix.at(vertex1).at(vertex2)) {
          std::list<cell>::iterator it = it_v; // Copy the iterator
          it_v++; // Advance the iterator before deleting the edge.
          remove_edge(g, adj, it);
          number_of_removed_edges++;
        }
        else it_v++; // Advance the current iterator.
      }
      it++;
    }
    return number_of_removed_edges;
  }

  graph restrict_graph_copy(graph& g, std::set<int>& points)
  {
    graph result = empty_graph();

    for (std::set<int>::iterator it = points.begin(); it != points.end(); it++) {
      add_vertex(result, *it); // Add the vertex
      if (g.info[*it].first.terminal) {
        set_terminal(result, *it, true);
      }
      std::list<cell> adj = g.info[*it].second;
      for (std::list<cell>::iterator it_adj = adj.begin(); it_adj != adj.end(); it_adj++) // Iterate through adjacent points, and add the vertices that exist in our set of points to keep 
      {
        if (*it > it_adj->vertex && points.find(it_adj->vertex) != points.end()) {
          // If the adjacent vertex exists in the points set, and make sure we don't add the same edge twice by verifying that the index of the second iterator is lower than the first.
          add_edge(result, *it, it_adj->vertex, it_adj->weight);
        }
      }
    }

    return result;
  }

  std::set<int> get_set_of_all_vertices(graph& g)
  {
    std::set<int> points; // Set of points
    for (std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.begin(); it != g.info.end(); it++) {
      points.insert(it->first); // Insert the index.
    }
    return points;
  }

  std::set<int> get_terminals(graph& g)
  {
    std::set<int> terminals; // Contains the terminals

    for (std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.begin(); it != g.info.end(); it++) {
      if (it->second.first.terminal) {
        terminals.insert(it->first); // Insert the index.
      }
    }

    return terminals;
  }

  std::pair<std::set<int>, std::set<int>> get_terminal_and_branching_set(graph& g)
  {
    std::set<int> terminals, branching_points; // Self explanatory

    for (std::map<int, std::pair<label, std::list<cell>>>::iterator it = g.info.begin(); it != g.info.end(); it++) {
      if (it->second.first.terminal) {
        terminals.insert(it->first); // Insert the index.
      }
      else if (it->second.second.size() >= 3) { // We don't include terminals that are branching points in here.
        branching_points.insert(it->first);
      }
    }

    return std::pair<std::set<int>, std::set<int>>(terminals, branching_points);
  }

  double get_graph_cost(graph& g)
  {
    return g.total_cost; // The total cost is already stored in this variable
  }

  std::list<std::pair<int, double>> shortest_path(graph& g, std::map<int, std::map<int, double>>& min_dist_matrix, int a, int b)
  {
    std::list<std::pair<int, double>> path;
    int vertex = a;
    while (vertex != b) {
      //We traverse the neighbours of vertex  
      for (auto it = g.info[vertex].second.begin();it != g.info[vertex].second.end();it++) {
        //We assume that the weights are bigger that 10E-5 so that we
        //can safely check the equality between doubles via this inequality
        if (abs(min_dist_matrix[vertex][it->vertex] + min_dist_matrix[it->vertex][b] - min_dist_matrix[vertex][b]) < 1E-10) {
          path.push_back(std::pair<int, double>(it->vertex, it->weight));
          vertex = it->vertex;
          break;
        }
      }
    }
    return path;
  };

  std::pair<graph, graph> enumeration_steiner_tree(graph& g, graph& d, std::map<int,std::map<int,double>>& min_dist_matrix)
  {
    std::pair<std::set<int>, std::set<int>> KunionS = get_terminal_and_branching_set(g); // S union K discussed in the program.
    std::set<int>& K = KunionS.first, & S = KunionS.second; // Get K and S seperately

    graph steiner_min_tree_in_g = empty_graph(), steiner_min_tree_in_d = empty_graph();
    double coust_min_span_tree_in_d = -1;


    int number_of_terminals = K.size(); // We could use g or d too. Same thing.
    int number_of_branching_points = S.size();

    for (int subset_size = 0; subset_size <= number_of_terminals - 2; subset_size++) {
      std::vector<std::set<int>> possible_sets = Setutils::get_subsets(S, subset_size);
      for (std::set<int> subs : possible_sets) { // Go through all the listed subsets of size subset_size
        /*
        for (auto v: subs) {
          for (auto m: v) {
            std::cout << m << "-";
          }
          std::cout << std::endl;
        }
        */
        // subs is a possible subset.
        std::set<int> KunionSubset;
        std::set<int> copy_of_K(K);
        KunionSubset.merge(copy_of_K);
        KunionSubset.merge(subs); // Only C++17 functionality used.
        std::set<int>::iterator it = KunionSubset.begin();
        graph restricted_copy = restrict_graph_copy(d, KunionSubset); // restricts D to the terminals and the subset of branching
        graph min_span_tree = primm_mst(restricted_copy);
        double cost_of_tree = get_graph_cost(min_span_tree);
        if ((coust_min_span_tree_in_d < 0) || (cost_of_tree < coust_min_span_tree_in_d)) {
          steiner_min_tree_in_d = min_span_tree;
          coust_min_span_tree_in_d = cost_of_tree;
        }
      }
    }

    for (std::map<int, std::pair<label, std::list<cell>>>::iterator it=steiner_min_tree_in_d.info.begin(); it!=steiner_min_tree_in_d.info.end(); it++){
      std::list<cell>& adj = it->second.second; // Adjacency list of the point
      if (it->second.first.terminal) {
        add_vertex(steiner_min_tree_in_g, it->first);
        set_terminal(steiner_min_tree_in_g, it->first);
      }
      for (std::list<cell>::iterator vt=adj.begin(); vt!=adj.end(); vt++) {
        std::list<std::pair<int, double>> path = shortest_path(g, min_dist_matrix, it->first, vt->vertex);
        int current_vertex = it->first;
        add_vertex(steiner_min_tree_in_g, current_vertex);
        for (std::pair<int, double>next : path) { // Add the intermediate paths
          add_vertex(steiner_min_tree_in_g, next.first);
          add_edge_without_repetition(steiner_min_tree_in_g, current_vertex, next.first, next.second);
          current_vertex = next.first;
        }
      }
    }

    // TODO: Add the part where I calculate the inverse tree in G from the tree in D. Missing.
    return std::pair<graph, graph>(steiner_min_tree_in_g, steiner_min_tree_in_d);
  }

  graph shortest_heuristic_path_algorithm(graph& g, std::map<int, std::map<int, double>>& min_dist_matrix)
  {
    std::set<int> K = get_terminals(g);
    graph T = empty_graph();
    add_vertex(T, *K.begin());
    K.erase(*K.begin());
    while (!K.empty()) {
      int closest_terminal = *K.begin();
      int closest_vertex = T.info.begin()->first;
      double closest_dist = min_dist_matrix[closest_terminal][closest_vertex];
      for (auto it = K.begin();it != K.end();it++) {
        for (auto jt = T.info.begin(); jt != T.info.end(); jt++) {
          if (min_dist_matrix[jt->first][*it] < closest_dist) {
            closest_vertex = jt->first;
            closest_terminal = *it;
          }
        }
      }
      K.erase(closest_terminal);
      std::set<int> T_vertices = get_set_of_all_vertices(T);
      std::list<std::pair<int, double>> path = shortest_path(g, min_dist_matrix, closest_vertex, closest_terminal);
      for (auto v = path.begin();v != path.end(); v++) {
        T_vertices.insert(v->first);
      }

      graph restriction = restrict_graph_copy(g, T_vertices);
      T = primm_mst(restriction);
      remove_leafs(T);
    }
    return T;
  }

  graph distance_network_heuristic_algorithm(graph& g, graph& min_dist_graph, std::map<int, std::map<int, double>>& min_dist_matrix){
    std::set<int> K = get_terminals(g);
    graph restriction = restrict_graph_copy(min_dist_graph,K);
    graph T_D = primm_mst(restriction);
    std::set<int> set_of_vertices = get_set_of_all_vertices(T_D);
    std::list<std::pair<int,double>> path;
    for(std::map<int, std::pair<label, std::list<cell>>>::iterator it = T_D.info.begin();it!= T_D.info.end();it++){
      for(std::list<cell>::iterator jt=it->second.second.begin();jt!=it->second.second.end();jt++){
        path = shortest_path(g,min_dist_matrix,it->first,jt->vertex);
        for(auto v = path.begin();v!=path.end();v++){
          set_of_vertices.insert(v->first);
        }
        set_of_vertices.insert(it->first);
      }
    }
    restriction = restrict_graph_copy(g,set_of_vertices);
    graph T = primm_mst(restriction);
    remove_leafs(T);
    return T;
  }

}