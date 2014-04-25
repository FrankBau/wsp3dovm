#include "common.h"

// write tetrahedralization in vtk format
void write_vtk( Mesh &mesh, std::string filename);

// write singel source shortest path tree in vtk format
void write_shortest_path_tree_vtk
(
Graph &graph,
std::vector<GraphNode_descriptor>& predecessors,
std::vector<double>& distances,
std::string filename
);

// write a longest shortest path in vtk format
void write_shortest_path_max_vtk
(
Graph &graph,
std::vector<GraphNode_descriptor>& predecessors,
std::vector<double>& distances,
std::string filename
);
