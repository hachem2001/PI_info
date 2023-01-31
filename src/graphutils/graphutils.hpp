#pragma once // pragma once should work with almost all modern compilers.
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <fstream>

namespace Setutils {
	/// @brief Finds all subsets of size equal to k of set s given its size n.
	/// @param set the set from which we'll get the subsets
	/// @param n The number of elemenst in s
	/// @param k The maximum size of the subsets to find
	/// @return A vector containing all subsets of set, of size k.
	std::vector<std::set<int>> get_subsets(std::set<int>& set, int k); 
}

namespace Graphutils {

	/// @brief This is an intermediary definition used to make working with graph adjacency lists easier.
	struct cell {
		int vertex;
		double weight;
		std::list<cell>::iterator _p; // Points to the cell from vertex's adjacency lists that points to this cell as well.
	};

	/// @brief This is a label attached to a vertex. It could be the name of a city, a number, ...
	struct label {
		int index; // For the moment we just keep a number, and identify it to the index of creation of the graph, which should be unique
		bool terminal = 0; // Indicates if the label in question is a terminal or not.
	};

	/// This is an augmented adjacency list representattion - every element {y,w,p} of an adjacency list to a vertex x contains : 
	/// The index of the adjacent point, the weight connecting them, and a pointer to the triplet {x,w,p'} from y's adjacency list.
	/// The lists being double-linked, this would mean our deleting a link would be -slightly- faster.
	struct graph {
		int _number_of_created_vertices = 0;

		int number_of_vertices = 0;
		int number_of_edges = 0;
		double total_cost = 0; // Contains the sum of all weights in the graph.
		std::map<int, std::pair<label, std::list<cell>>> info;
		// For _number_of_created_vertices: this number goes up and never down, and does up when a new vertex is created.
		// It is currently used to modify the label to keep it unique for each vertex.
	};


	/// @brief Creates an empty graph
	/// @return Empty Graph.
	graph empty_graph();

	/// @brief Attemps to create a graph from an ifstream (typically a file that is read for example). This correspond's to the scenarios format
	/// @param stream Stream to transform into graph
	/// @return Transformed graph
	graph fstream_graph(std::ifstream& stream);

	/// @brief Adds a vertex. Returns the number corresponding to the vertex
	/// @param g Graph
	int add_vertex(graph& g);

	/// @brief Adds a vertex with the specified index. Returns the same number of the index.
	/// @param g Graph
	/// @param index Index to create
	/// @return 
	int add_vertex(graph& g, int index);

	/// @brief Overload : adds n vertices to the graph g
	/// @param g Graph
	/// @param n number of vertices to add
	/// @return index of the first added vertex.
	int add_vertices(graph& g, int n);

	/// @brief Disconnects a vertex from the graph. It will therefore become connected to NO other vertex.
	/// @param g Graph
	/// @param vertex Vertex number (from 0 to number_of_vertices-1)
	void disconnect_vertex(graph& g, int vertex);

	/// @brief Removes an index and all the weights connecting to it. 
	/// @param g Graph
	/// @param vertex Vertex number (from 0 to number_of_vertices-1)
	void remove_vertex(graph& g, int vertex);

	/// @brief Adds an edge to the Graph g linking vertex1 and vertex2 with the cost weight 
	/// @param g Graph
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	/// @param weight Weight between both vertices.
	void add_edge(graph& g, int vertex1, int vertex2, double weight);

	/// @brief Adds the edge if it does not exist yet.
	/// @param g Graph g
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	/// @param weight Total weight
	void add_edge_without_repetition(graph& g, int vertex1, int vertex2, double weight);

	/// @brief Gets the cell link from vertex1 to vertex2
	/// @param g Graph
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	/// @return <true, cell linking 1 to 2> if an edge exists, <false, empty cell> otherwise
	std::pair<bool, std::list<cell>::iterator> get_edge(graph& g, int vertex1, int vertex2);

	/// @brief Gets the weight from vertex 1 to vertex 2. 
	/// @param g Graph
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	/// @return Returns <true, weight> if edge exists. <false, -1> otherwise.
	std::pair<bool, double> get_weight(graph& g, int vertex1, int vertex2);

	/// @brief Removes the edge (vertex1,vertex2) from the graph
	/// @param g Graph
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	void remove_edge(graph& g, int vertex1, int vertex2);

	/// @brief Removes the edge corresponding to the iterator it. Runs in O(1) compared to the previous function
	/// @param g Graph.
	/// @param adj Adjacency list of the vertex containing the iterator
	/// @param it The iterator in question
	void remove_edge(graph& g, std::list<cell>& adj, std::list<cell>::iterator& it);

	/// @brief Adds the vertices if they are missing, and connect them with the cost weight.
	/// @param g Graph
	/// @param vertex1 Vertex 1
	/// @param vertex2 Vertex 2
	/// @param weight Weight
	void add_link(graph& g, int vertex1, int vertex2, double weight);

	/// @brief Adds the vertices if missing, and the edge if not existant.
	/// @param g Graph
	/// @param vertex1 Vertex 1
	/// @param vertex2 Vertex 2
	/// @param weight Weight
	void add_link_without_repetition(graph& g, int vertex1, int vertex2, double weight);

	/// @brief Sets the vertex as a terminal (or non terminal) (part of K) in the graph
	/// @param g Graph
	/// @param vertex The vertex to be 
	/// @param set Boolean, true or false;
	void set_terminal(graph&, int vertex, bool set);

	// HERE WE START WORKING ON THE QUESTIONS ASKED IN THE PI PDF PROPERLY.

	/// @brief Deletes leafs from the graph that are not labeled terminals
	/// @param g Graph
	/// @return Returns the number of leafs removed.
	int remove_leafs(graph& g);

	/// @brief Remove third edges if degree 2 non terminal vertices if sum of weights smaller than third edge, Or remove edge if otherwise. If non existant third edge : do nothing.
	/// @param g Graph
	std::pair<int, int> optimize_degree_2(graph& g);

	/// @brief Minimum distance array from vertex to other vertices
	/// @param g Graph g
	/// @param vertex Source verterx
	/// @return Vector containing distance to vertex
	std::map<int, double> min_distance_to_source(graph& g, int vertex);

	/// @brief Minimum distance matrix 
	/// @param g 
	/// @param vertex1 
	/// @param vertex2 
	/// @return 
	std::map<int, std::map<int, double>> min_distance(graph& g);

	/// @brief Transforms the minimum distance matrix into a graph.
	/// @param min_distance_matrix The minimum distance matrix
	/// @return Minimum distance graph D
	graph min_distance_graph(graph& g, std::map<int, std::map<int, double>> min_distance_matrix);

	/// @brief Calculates the minimum spanning tree of G using primm's algorithm
	/// @param g Connected input graph. (Must be connected!)
	/// @return MST graph.
	graph primm_mst(graph& g);

	/// @brief Remove edges that are longer than the minimum distance linking their two edges
	/// @param g Graph to pre-process
	/// @param min_dist_matrix the minimum distance matrix calculated with min_distance. 
	int optimize_edge_with_min_distance_matrix(graph& g, std::map<int, std::map<int, double>> min_dist_matrix);

	/// @brief Restricts the graph to the corresponding points.
	/// @param g Graph to restrict
	/// @param points A set containing the indices of the vertices to keep.
	/// @return A graph copy containing only the indices.
	graph restrict_graph_copy(graph& g, std::set<int>& points);

	/// @brief Takes graph as input, gives out a set of all the indices (of all the vertices)
	/// @param g Input graph
	/// @return std::set<int> of all the vertices as integers.
	std::set<int> get_set_of_all_vertices(graph& g);

	/// @brief Returns a set of the indices corresponding to the terminals
	/// @param g Graph g.
	/// @return std::set<int> with the indices
	std::set<int> get_terminals(graph& g);

	/// @brief Gives a set containing the indices of the terminals and branching points in g.
	/// @param g Graph to go through
	/// @return Returns std::pair<set of terminals, set of branching points>, each set containing the indices.
	std::pair<std::set<int>, std::set<int>> get_terminal_and_branching_set(graph& g);

	/// @brief Returns the sum of all weights in the graph.
	/// @param g Input graph
	/// @return Graph cost.
	double get_graph_cost(graph& g);

	/// @brief Given the minimum distance matrix, it gives out the shortest path corresponding to it.
	/// @param min_dist_matrix The minimum distance matrix
	/// @param G the original graph G
	/// @param a The first vertex
	/// @param b The second vertex
	/// @return Shortest path between a and b
	std::list<std::pair<int,double>> shortest_path(graph& g, std::map<int,std::map<int,double>>& min_dist_matrix, int a,int b);

	/// @brief Returns the steiner minimal tree in G and D, D being the minimum distance graph of G.
	/// @param g Original Graph
	/// @param d Minimum distance grance
	/// @return std::pair<graph1, graph2> with graph1 being the steiner minimal tree in g, and graph2 being the steiner minimal tree in d.
	std::pair<graph, graph> enumeration_steiner_tree(graph& g, graph& d, std::map<int,std::map<int,double>>& min_dist_matrix);

	/// @brief Returns the graph corresponding Steiner tree in G given the one in the minimum distance matrix
	/// @param g Original graph
	/// @param steiner_in_d The steiner tree graph in D.
	/// @param min_dist_matrix Minimum distance matrix
	/// @return The steiner tree in G
	graph min_in_d_to_g(graph& g, graph& steiner_in_d, std::map<int,std::map<int,double>>& min_dist_matrix);

	/// @brief Uses the heuristic method to calculate the steiner Tree.
	/// @param g Graph g
	/// @param min_dist_matrix The min distance matrix corresponding to g.
	/// @return A graph that is a steiner tree, approximating the Steiner minimal tree.
	graph shortest_heuristic_path_algorithm(graph &g, std::map<int,std::map<int,double>>& min_dist_matrix);

	/// @brief Uses the first heuristic algorithm mentioned in the PDF
	/// @param g Graph g
	/// @param min_dist_graph Corresponding minimum distance graph.
	/// @param min_dist_matrix Corresponding minimum distance matrix
	/// @return Approximation of the Steiner minimal Tree
	graph distance_network_heuristic_algorithm(graph& g, graph& min_dist_graph, std::map<int, std::map<int, double>>& min_dist_matrix);

	/// @brief Dreyfus_Wagner algorithm
	/// @param g Graph g
	/// @param min_dist_matrix Corresponding minimum distance matrix
	/// @return Exact computation of the Steiner minimal Tree	
	double dreyfus_wagner_algorithm(graph& g, std::map<int,std::map<int,double>>& min_dist_matrix);

	/// @brief Allows to sort integers from 0 to 2^(n-1) sorted along the number of 1's
	/// @param n integer
	/// @return set 
	std::vector<std::pair<int,int>> sorted_set_int(int n);

	/// @brief Allows to convert an integer to a set
	/// @param super_set the bigger set sorted
	/// @param number number representation
	/// @return set representation
	std::set<int> decode_set(std::set<int> super_set,int number);

	/// @brief Allows to convert a set to an integer
	/// @param set subset
	/// @param k number of elements in the superset
	/// @return integer representation.
	int encode_set(std::set<int> set, int k);

	/// @brief Allows to print the graph in a readable format (for std::cout for example)
	/// @param os Output stream
	/// @param m The graph to print
	/// @return Output stream.
	std::ostream& operator<<(std::ostream& os, graph const& m);

	/// @brief Allows to print the mininal distance matrix in a readable format.
	/// @param os Output stream
	/// @param min_matrix The matrix in question
	/// @return Output stream.
	std::ostream& operator<<(std::ostream& os, std::map<int, std::map<int, double>> const& min_matrix);

}
