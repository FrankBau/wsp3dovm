#include "write_vtk.h"

void write_shortest_path_tree_vtk
(
	Graph &graph,
	std::vector<GraphNode_descriptor>& predecessors,
	std::vector<double>& distances,
	std::string filename
)
{
	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	int n = boost::num_vertices(graph);

	file <<
		"# vtk DataFile Version 2.0\n"
		"shortest path tree\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n";
	file << "POINTS " << n << " double\n";

	std::vector<int> id(num_vertices(graph));
	int i = 0;
	Graph::vertex_iterator vertexIt, vertexEnd;
	for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph); vertexIt != vertexEnd; ++vertexIt)
	{
		GraphNode_descriptor v = *vertexIt;
		id[v] = i++;
		file << graph[v].point << "\n";
	}

	file << "CELLS " << n << " " << 3 * n << "\n";
	for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph); vertexIt != vertexEnd; ++vertexIt)
	{
		GraphNode_descriptor u = *vertexIt;
		GraphNode_descriptor v = predecessors[u];
		file << "2 " << id[u] << " " << id[v] << "\n";
	}

	file << "CELL_TYPES " << n << "\n";
	for (int i = 0; i < n; ++i)
	{
		// vtk cell type 3 is line
		file << "3" "\n";
	}

	file
		<< "POINT_DATA " << n << "\n"
		<< "SCALARS distance double 1\n"
		<< "LOOKUP_TABLE default\n";

	for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph); vertexIt != vertexEnd; ++vertexIt)
	{
		GraphNode_descriptor u = *vertexIt;
		file << distances[u] << "\n";
	}
	file << "\n";
}

int hops(std::vector<GraphNode_descriptor> &predecessors, GraphNode_descriptor t)
{
	int h = 0;
	while (t != predecessors[t])
	{
		t = predecessors[t];
		++h;
	}
	assert(t == 0);
	return h;
}

void write_shortest_path_max_vtk
(
	Graph &graph,
	std::vector<GraphNode_descriptor>& predecessors,
	std::vector<double>& distances,
	std::string filename
)
{
	int s = 0; // source node index
	int t = s; // graph node having max. distance, initially t==s
	double max_distance = distances[t];
	assert(max_distance == 0.0);

	for (int i = 0; i < distances.size(); ++i)
	{
		if (distances[i] > max_distance)
		{
			t = i;
			max_distance = distances[i];
		}
	}
	std::cout << "s = " << graph[s].cell.idx() << "; t = " << graph[t].cell.idx() << std::endl;

	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	int n = boost::num_vertices(graph);

	file <<
		"# vtk DataFile Version 2.0\n"
		"longest shortest path\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n";
	file << "POINTS " << n << " double\n";

	std::vector<int> id(num_vertices(graph));
	int i = 0;
	Graph::vertex_iterator vertexIt, vertexEnd;
	for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph); vertexIt != vertexEnd; ++vertexIt)
	{
		GraphNode_descriptor u = *vertexIt;
		id[u] = i++;
		file << graph[u].point << "\n";
	}

	int h = hops(predecessors, t);

	file << "CELLS " << h << " " << 3 * h << "\n";
	for (GraphNode_descriptor r = t; r != s; r = predecessors[r])
	{
		GraphNode_descriptor v = predecessors[r];
		file << "2 " << id[r] << " " << id[v] << "\n";
	}

	file << "CELL_TYPES " << h << "\n";
	for (int i = 0; i < h; ++i)
	{
		// vtk cell type 3 is line
		file << "3" "\n";
	}

	file
		<< "POINT_DATA " << n << "\n"
		<< "SCALARS distance double 1\n"
		<< "LOOKUP_TABLE default\n";

	for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph); vertexIt != vertexEnd; ++vertexIt)
	{
		GraphNode_descriptor u = *vertexIt;
		file << distances[u] << "\n";
	}
	file << "\n";
}

void write_vtk(Mesh &mesh, std::string filename)
{
	std::cout << "write_vtk " << filename << std::endl;
	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}


	int number_of_vertices = mesh.n_vertices();

	file <<
		"# vtk DataFile Version 2.0\n"
		"tetrahedralization\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n";

	file << "POINTS " << number_of_vertices << " double\n";

	std::map<VertexHandle, int> V;
	int inum = 0;

	for (auto vit = mesh.vertices_begin(), end = mesh.vertices_end();
		vit != end;
		++vit
	)
	{
		Point point = mesh.vertex(*vit);

		file
			<< point[0] << " "
			<< point[1] << " "
			<< point[2] << " "
			<< "\n";

		V[*vit] = inum++;
	}

	int number_of_cells = mesh.n_cells();

	// for each cell 5 values are given: numPoints (4) and the 4 vertex indices
	file << "CELLS " << number_of_cells << " " << 5 * number_of_cells << "\n";

	for (
		auto cit = mesh.cells_begin(), end = mesh.cells_end();
		cit != end;
		++cit)
	{
		CellHandle ch = *cit;

		file << "4 ";
		for (auto vit = mesh.cv_iter(ch); vit; ++vit)
			file << vit->idx() << " ";
		
		file << "\n";
	}

	file << "CELL_TYPES " << number_of_cells << "\n";
	for (int i = 0; i < number_of_cells; ++i)
	{
		// vtk cell type 10 is tetrahedron
		file << "10" "\n";
	}

	file
		<< "CELL_DATA " << number_of_cells << "\n"
		<< "SCALARS weight double 1\n"
		<< "LOOKUP_TABLE default\n";

	for (
		auto cit = mesh.cells_begin(), end = mesh.cells_end();
		cit != end;
	++cit)
	{
		file << mesh.weight(*cit);
		file << "\n";
	}

	file.close();
}
