#pragma once // pragma once should work with almost all modern compilers.
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <fstream>

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

	/// @brief Adds the vertices if they are missing, and connect them with the cost weight.
	/// @param g Graph
	/// @param vertex1 Vertex 1
	/// @param vertex2 Vertex 2
	/// @param weight Weight
	void add_link(graph& g, int vertex1, int vertex2, double weight);

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

	/// @brief Allows to print the graph in a readable format (for std::cout for example)
	/// @param os Output stream
	/// @param m The graph to print
	/// @return Output stream.
	std::ostream& operator<<(std::ostream& os, graph const& m);
}
