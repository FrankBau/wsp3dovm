#include "common.h"

// write tetrahedralization in vtk format
void write_vtk(const Mesh &mesh, std::string filename);

// write singel source shortest path tree in vtk format
void write_shortest_path_tree_vtk
(
	const Graph &graph,
	std::vector<GraphNode_descriptor>& predecessors,
	std::vector<double>& distances,
	std::string filename
);

// write a longest shortest path in vtk format
void write_shortest_path_to_vtk
(
	const Graph &graph,
	GraphNode_descriptor t,
	std::vector<GraphNode_descriptor>& predecessors,
	std::vector<double>& distances,
	std::string filename
);

// write complete graph, mainly for debugging small examples
void write_graph_vtk
(
	const Graph &graph,
	std::string filename
);

