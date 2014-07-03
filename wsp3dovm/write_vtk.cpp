#include "write_vtk.h"

void write_shortest_path_tree_vtk
(
	const Graph &graph,
	GraphNode_descriptor s,
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

	size_t n = boost::num_vertices(graph);

	file <<
		"# vtk DataFile Version 2.0\n"
		"shortest paths tree with root node " << graph[s].vh.idx() << "\n"
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
	for (size_t i = 0; i < n; ++i)
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

int hops(std::vector<GraphNode_descriptor> &predecessors, GraphNode_descriptor s, GraphNode_descriptor t)
{
	int h = 0;
	while (s != t)
	{
		t = predecessors[t];
		++h;
	}
	return h;
}

void write_shortest_path_from_to_vtk
(
	const Graph &graph,
	GraphNode_descriptor s, 
	GraphNode_descriptor t,
	std::vector<GraphNode_descriptor>& predecessors,
	std::vector<double>& distances,
	std::string filename
)
{
	if (distances[t] == std::numeric_limits<double>::infinity())
	{
		std::cout << "write_shortest_path_from_to_vtk: target node unreachable" << std::endl;
	}

	int h = hops(predecessors, s, t);

	std::cout 
	<< "from s=" << graph[s].vh.idx() 
	<< " to t=" << graph[t].vh.idx() 
	<< " distance=" << distances[t]
	<< " #hops = " << h 
	<< std::endl;

	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	size_t n = boost::num_vertices(graph);

	file <<
		"# vtk DataFile Version 2.0\n"
		"shortest path from " << graph[s].vh.idx() << " to " << graph[t].vh.idx() << "\n"
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

void write_vtk(const Mesh &mesh, std::string filename)
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

void write_graph_vtk
(
const Graph &graph,
std::string filename
)
{
	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	file <<
		"# vtk DataFile Version 2.0\n"
		"steiner graph\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n";
	file << "POINTS " << num_vertices(graph) << " double\n";

	std::vector<int> id(num_vertices(graph));
	int i = 0;
	Graph::vertex_iterator vertexIt, vertexEnd;
	for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph); vertexIt != vertexEnd; ++vertexIt)
	//for (auto v : graph.m_vertices )
	{
		GraphNode_descriptor u = *vertexIt;
		id[u] = i++;
		file << graph[u].point << "\n";
	}

	file << "CELLS " << num_edges(graph) << " " << 3 * num_edges(graph) << "\n";
	for (auto e : graph.m_edges)
	{
		GraphNode_descriptor u = e.m_source;
		GraphNode_descriptor v = e.m_target;
		file << "2 " << id[u] << " " << id[v] << "\n";
	}

	file << "CELL_TYPES " << num_edges(graph) << "\n";
	for (size_t i = 0; i < num_edges(graph); ++i)
	{
		// vtk cell type 3 is line
		file << "3" "\n";
	}

	// dump weight (not cost!) of edges to check that the adjacent  fetures have been calculated correctly
	file
		<< "CELL_DATA " << num_edges(graph) << "\n"
		<< "SCALARS edge_weight double 1\n"
		<< "LOOKUP_TABLE default\n";

	Graph::edge_iterator edgeIt, edgeEnd;
	for (boost::tie(edgeIt, edgeEnd) = boost::edges(graph); edgeIt != edgeEnd; ++edgeIt)
	{
		GraphEdge_descriptor edge = *edgeIt;
		file << graph[edge].weight << "\n";
	}
	file << "\n";

	file.close();
}

std::set<CellHandle> cells_from_graph_nodes
(
	const Graph& graph, 
	const Mesh &mesh,
	GraphNode_descriptor s,
	GraphNode_descriptor t,
	std::vector<GraphNode_descriptor>& predecessors
)
{
	std::set<CellHandle> cells;

	GraphNode_descriptor r = t;
	for (;;)
	{
		if (graph[r].vh != OpenVolumeMesh::TopologyKernel::InvalidVertexHandle)
		{
			VertexHandle vh = graph[r].vh;

			std::cout << "vrtx node " << vh.idx() << std::endl;
			for (auto vc_iter = mesh.vc_iter(vh); vc_iter.valid(); ++vc_iter )
			{
				CellHandle ch = *vc_iter;
				if (ch.is_valid())
				{
					std::cout << "  incident cell: " << ch.idx() << std::endl;
					cells.insert(ch);
				}
			}
		}
		else if (graph[r].eh != OpenVolumeMesh::TopologyKernel::InvalidEdgeHandle)
		{
			EdgeHandle eh = graph[r].eh;
			std::cout << "edge node " << eh.idx() << std::endl;

			// we can choose either halfegde, both yield the same set (tested)
			HalfEdgeHandle heh0 = mesh.halfedge_handle(eh, 0);
			for (auto hec_iter = mesh.hec_iter(heh0); hec_iter.valid(); ++hec_iter)
			{
				CellHandle ch = *hec_iter;
				if (ch.is_valid())
				{
					std::cout << "  incident cell: " << ch.idx() << std::endl;
					cells.insert(ch);
				}
			}
		}
		else if (graph[r].fh != OpenVolumeMesh::TopologyKernel::InvalidFaceHandle)
		{
			FaceHandle fh = graph[r].fh;
			std::cout << "face node " << fh.idx() << std::endl;
			
			HalfFaceHandle hfh0 = mesh.halfface_handle(fh, 0);
			CellHandle ch0 = mesh.incident_cell(hfh0);
			if (ch0.is_valid())
			{
				cells.insert(ch0);
			}

			HalfFaceHandle hfh1 = mesh.halfface_handle(fh, 1);
			CellHandle ch1 = mesh.incident_cell(hfh1);
			if (ch1.is_valid())
			{
				cells.insert(ch1);
			}
		}
		else
		{
			assert(0 && "graph node should be vrtx, edge or face node.");
		}

		if (r == s)
			break; 
		r = predecessors[r];
	}

	std::cout << cells.size() << " cells out of " << mesh.n_cells() << " included in shortest path subcomplex" << std::endl;

	return cells;
}

void write_shortest_path_cells_from_to_vtk
(
	const Graph &graph,
	const Mesh &mesh,
	GraphNode_descriptor s,
	GraphNode_descriptor t,
	std::vector<GraphNode_descriptor>& predecessors,
	std::vector<double>& distances,
	std::string filename
)
{
	if (distances[t] == std::numeric_limits<double>::infinity())
	{
		std::cout << "write_shortest_path_cells_from_to_vtk: target node unreachable" << std::endl;
	}

	std::set<CellHandle> cells = cells_from_graph_nodes(graph, mesh, s, t, predecessors);

	int h = hops(predecessors, s, t);

	std::cout
		<< "from s=" << graph[s].vh.idx()
		<< " to t=" << graph[t].vh.idx()
		<< " distance=" << distances[t]
		<< " #hops = " << h
		<< std::endl;

	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	size_t n = boost::num_vertices(graph);

	file <<
		"# vtk DataFile Version 2.0\n"
		"cells along shortest path from " << graph[s].vh.idx() << " to " << graph[t].vh.idx() << "\n"
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

	file << "CELLS " << cells.size() << " " << 5 * cells.size() << "\n";
	for (
		auto cit =cells.begin(), end = cells.end();
		cit != end;
		++cit)
	{
		CellHandle ch = *cit;

		file << "4 ";
		for (auto vit = mesh.cv_iter(ch); vit; ++vit)
			file << vit->idx() << " ";

		file << "\n";
	}

	file << "CELL_TYPES " << cells.size() << "\n";
	for (int i = 0; i < cells.size(); ++i)
	{
		// vtk cell type 10 is tetra
		file << "10" "\n";
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
