#pragma once // pragma once should work with almost all modern compilers.
#include <iostream>
#include <vector>
#include <list>



namespace Graphutils {

	/// @brief This is an intermediary definition used to make working with graph adjacency lists easier.
	struct cell  {
		int vertex;
		double weight;
		void* _p; // Points to the cell from vertex's adjacency lists that points to this cell as well.
	};
	
	/// @brief This is a label attached to a vertex. It could be the name of a city, a number, ...
	struct label {
		int index; // For the moment we just keep a number, and identify it to the index of creation of the graph, which should be unique
	};

	/// This is an augmented adjacency list representattion - every element {y,w,p} of an adjacency list to a vertex x contains : 
	/// The index of the adjacent point, the weight connecting them, and a pointer to the triplet {x,w,p'} from y's adjacency list.
	/// The lists being double-linked, this would mean our deleting a link would be -slightly- faster.
	struct graph {
		int _number_of_created_vertices = 0; 
		int number_of_vertices = 0;
		int number_of_edges = 0;
		std::vector<std::pair<label, std::list<cell>>> info; // For each vertex, we associate a list of (connected_vertex, weight).
		// I use void* because a recursive definition of the types would be much more messy.
		// We therefore try to be careful when verifying the integrity of the code.
		// For _number_of_created_vertices: this number goes up and never down, and does up when a new vertex is created.
		// It is currently used to modify the label to keep it unique for each vertex.
	};

	/// @brief Creates an empty graph
	/// @return Empty Graph.
	graph empty_graph();

	/// @brief Adds a vertex. Returns the number corresponding to the vertex
	/// @param g 
	int add_vertex(graph& g);

	/// @brief Removes an index and all the weights connecting to it. 
	/// @param g Graph
	/// @param vertex Vertex number (from 0 to number_of_vertices-1)
	void remove_vertex(graph& g, int vertex);
	
	/// @brief Disconnects a vertex from the graph. It will therefore become connected to NO other vertex.
	/// @param g Graph
	/// @param vertex Vertex number (from 0 to number_of_vertices-1)
	void disconnect_vertex(graph& g, int vertex);

	/// @brief Adds an edge to the Graph g linking vertex1 and vertex2 with the cost weight 
	/// @param g Graph
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	/// @param weight Weight between both vertices.
	void add_edge(graph& g, int vertex1, int vertex2, double weight);

	/// @brief Removes the edge (vertex1,vertex2) from the graph
	/// @param g Graph
	/// @param vertex1 First vertex
	/// @param vertex2 Second vertex
	void remove_edge(graph& g, int vertex1, int vertex2);

	/// @brief Adds the vertices that are missing. 
	/// @param g 
	/// @param vertex1 
	/// @param vertex2 
	/// @param weight 
	void add_link(graph& g, int vertex1, int vertex2, double weight);
	/// @brief Allows to print the graph in a readable format (for std::cout for example)
	/// @param os Output stream
	/// @param m The graph to print
	/// @return Output stream.
	std::ostream &operator<<(std::ostream &os, graph const &m);
}

