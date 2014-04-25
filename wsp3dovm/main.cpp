#include "common.h"

#include "read_tet.h"
#include "write_vtk.h"
#include "create_steinerpoints.h"

void dump(Mesh &mesh)
{
	std::cout << "Vertices " << mesh.n_vertices() << std::endl;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		std::cout << "  Position of vertex " << v_it->idx() << ": " << mesh.vertex(*v_it) << std::endl;
	}

	std::cout << "Edges " << mesh.n_edges() << std::endl;
	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		std::cout << "  Vertices of edge " << e_it->idx() << ": " << mesh.edge(*e_it) << "  weight of edge " << mesh.weight(*e_it) << std::endl;
		// std::cout << "  Vertices of edge " << e_it->idx() << ": " << mesh.edge(*e_it) << std::endl;
	}

	std::cout << "Facets " << mesh.n_faces() << std::endl;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::cout << "  Halfedges of face " << f_it->idx() << ": " << mesh.face(*f_it) << "  weight of face " << mesh.weight(*f_it) << std::endl;
		//std::cout << "  Halfedges of face " << f_it->idx() << ": " << mesh.face(*f_it) << std::endl;
	}

	std::cout << "Cells " << mesh.n_cells() << std::endl;
	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		std::cout << "  Halffaces of cell " << c_it->idx() << ": " << mesh.cell(*c_it) << "  weight of cell " << mesh.weight(*c_it) << std::endl;
		// std::cout << "  Halffaces of cell " << c_it->idx() << ": " << mesh.cell(*c_it) << std::endl;
	}
}

void set_cell_weights(Mesh &mesh)
{
	boost::mt19937 rng;
	boost::uniform_real<> unity(0.0, 1.0);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > get_random_weight(rng, unity);
	int i = 1;

	mesh._cellWeight.resize(mesh.n_cells());

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;
		mesh.weight(ch) = get_random_weight();
		//mesh.weight(ch) = i++;
	}
}

void calc_face_weights(Mesh &mesh)
{
	mesh._faceWeight.resize(mesh.n_faces());

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		OpenVolumeMesh::FaceHandle fh = *f_it;

		OpenVolumeMesh::HalfFaceHandle hfh0 = Kernel::halfface_handle(fh, 0);
		OpenVolumeMesh::HalfFaceHandle hfh1 = Kernel::halfface_handle(fh, 1);

		CellHandle cell0 = mesh.incident_cell(hfh0);
		CellHandle cell1 = mesh.incident_cell(hfh1);

		double w0 = cell0.is_valid() ? mesh.weight(cell0) : max_weight;
		double w1 = cell1.is_valid() ? mesh.weight(cell1) : max_weight;

		mesh.weight(fh) = std::min(w0, w1);
	}
}

void calc_edge_weights(Mesh &mesh)
{
	mesh._edgeWeight.resize(mesh.n_edges());

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		OpenVolumeMesh::EdgeHandle eh = *e_it;

		//std::cout << "processing edge " << e_it->idx() << " " << mesh.edge(*e_it) << "of valence " << mesh.valence(*e_it) << std::endl;

		OpenVolumeMesh::HalfEdgeHandle heh0 = Kernel::halfedge_handle(eh, 0);
		//OpenVolumeMesh::HalfEdgeHandle heh1 = Kernel::halfedge_handle(eh, 1);

		double weight = max_weight;

		for( auto hec0 = mesh.hec_iter(heh0); hec0; ++hec0 )
		{
			//std::cout << "  hec0 incident cell " << hec0->idx() << std::endl;
			if (hec0->is_valid())
				weight = std::min(weight, mesh.weight(*hec0));
		}

		// this yields exactly the same cells in reverse order

		//for (auto hec1 = mesh.hec_iter(heh1); hec1; ++hec1)
		//{
		//	std::cout << "  hec1 incident cell " << hec1->idx() << std::endl;
		//if (hec1->is_valid())
		//		weight = std::min(weight, mesh.weight(*hec1));
		//}

		mesh.weight(eh) = weight;
	}
}

int main(int _argc, char** _argv)
{
	Mesh mesh;

	//boost::filesystem::path inputfilename("../../Meshes/holmes_off/geometry/tetrahedron.1");
	boost::filesystem::path inputfilename("../..//Meshes/misc/cube.1");
	//boost::filesystem::path inputfilename("../..//Meshes/misc/moomoo.1");

	//boost::filesystem::path inputfilename("../../Meshes/GEN_IV_plenum/NewFineMesh.1");
	// output:
	//reading 51944 nodes... done.
	//reading 172284 tetras... done.
	//Detected non-three-manifold configuration!
	//Connectivity probably won't work.
	//Detected non-three-manifold configuration!
	//Connectivity probably won't work.
	//read_tet finished.
	//write_vtk NewFineMesh.vtk

	read_tet(mesh, inputfilename.string() );

	set_cell_weights(mesh);
	calc_face_weights(mesh);
	calc_edge_weights(mesh);

	//dump(mesh);

	write_vtk(mesh, inputfilename.filename().replace_extension(".vtk").string() );

	Graph graph;

	create_steiner_points(graph, mesh);

	// the distances are temporary, so we choose an external property for that
	std::vector<double> distances(num_vertices(graph));
	std::vector<GraphNode_descriptor> predecessors(num_vertices(graph));

	boost::dijkstra_shortest_paths
	(
		graph,
		*vertices(graph).first,
		boost::weight_map(get(&GraphEdge::weight, graph)).
		distance_map(boost::make_iterator_property_map(distances.begin(), get(boost::vertex_index, graph))).
		predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(boost::vertex_index, graph)))
	);

	write_shortest_path_tree_vtk(graph, predecessors, distances, inputfilename.filename().replace_extension("_wsp_tree.vtk").string());

	write_shortest_path_max_vtk(graph, predecessors, distances, inputfilename.filename().replace_extension("_wsp_max.vtk").string());

	// write_graph_dot("graph.dot", graph);

	std::cout << "This is the end..." << std::endl;

	return EXIT_SUCCESS;
}
