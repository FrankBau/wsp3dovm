#include "common.h"

#include "read_tet.h"
#include "write_vtk.h"
#include "create_steinerpoints.h"
#include "statistics.h"

using namespace std;
using namespace boost;
using namespace boost::chrono;

template< class Clock >
class timer
{
  typename Clock::time_point start;
public:
  timer() : start( Clock::now() ) {}
  typename Clock::duration elapsed() const
  {
    return Clock::now() - start;
  }
  double seconds() const
  {
    return elapsed().count() * ((double)Clock::period::num/Clock::period::den);
  }
};


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
		//mesh.weight(ch) = get_random_weight();
		mesh.weight(ch) = 1;
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

		Weight weight = max_weight;

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

int main(int argc, char** argv)
{
	timer<high_resolution_clock> total_time;

	std::cout << "sizeof an int:   " << sizeof(int) << std::endl;
	std::cout << "sizeof a void*:  " << sizeof(void*) << std::endl;
	std::cout << "sizeof a Handle: " << sizeof(EdgeHandle) << std::endl;

	int start_vertex; // for single source shortest paths (Dijkstra)
	int termination_vertex; // an optional termination vertex for which the shortest path will be reported
	double stretch; // spaner graph stretch factor

	program_options::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("start_vertex,s", program_options::value<int>(&start_vertex)->default_value(0), "shortest path start vertex number")
		("termination_vertex,t", program_options::value<int>(&termination_vertex)->default_value(-1), "shortest path termination vertex number (-1==none)")
		("spanner_stretch,x", program_options::value<double>(&stretch)->default_value(0.0), "spanner graph stretch factor")
		("input-mesh", program_options::value<std::string>(), "set input filename (tetgen 3D mesh files wo extension)")
		;

	program_options::positional_options_description positional_options;
	positional_options.add("input-mesh", 1);

	program_options::variables_map vm;
	program_options::store(
		program_options::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), 
		vm
	);
	program_options::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << endl;
		return EXIT_FAILURE;
	}

	if (vm.count("input-mesh"))
	{
		cout << "Input mesh: " << vm["input-mesh"].as<string>() << endl;
	}
	else
	{
		cout << "no input-mesh specified, exit." << endl;
		std::cout << desc << endl;
		return EXIT_FAILURE;
	}

	boost::filesystem::path inputfilename(vm["input-mesh"].as<string>());

	Mesh mesh;

	timer<high_resolution_clock> t;
	read_tet(mesh, inputfilename.string() );
	std::cout << "read_tet: " << t.seconds() << " s" << std::endl;

	set_cell_weights(mesh);
	calc_face_weights(mesh);
	calc_edge_weights(mesh);

	mesh.print_memory_statistics();
	print_mesh_statistics(mesh);

	{
		timer<high_resolution_clock> t;
		write_vtk(mesh, inputfilename.filename().replace_extension(".vtk").string() );
		std::cout << "write_vtk: " << t.seconds() << " s" << std::endl;
	}
	
	Graph graph;

	{
		timer<high_resolution_clock> t;
		
		//create_barycentric_steiner_points(graph, mesh);
		//std::cout << "create_barycentric_steiner_points: " << t.seconds() << " s" << std::endl;

		if (stretch = 0.0)
		{
			create_surface_steiner_points(graph, mesh);
			std::cout << "create_surface_steiner_points: " << t.seconds() << " s" << std::endl;
		}
		else
		{
			create_steiner_graph_improved_spanner(graph, mesh);
			std::cout << "create_steiner_graph_improved_spanner: " << t.seconds() << " s" << std::endl;
		}
	}
	
	std::cout << "graph nodes: " << graph.m_vertices.size() << std::endl;
	std::cout << "graph edges: " << graph.m_edges.size() << std::endl;

	{
		timer<high_resolution_clock> t;
		write_graph_vtk(graph, inputfilename.filename().replace_extension("_steiner_graph.vtk").string());
		std::cout << "write_graph_vtk: " << t.seconds() << " s" << std::endl;
	}

	// the distances are temporary, so we choose an external property for that
	std::vector<double> distances(num_vertices(graph));
	std::vector<GraphNode_descriptor> predecessors(num_vertices(graph));

	{
		timer<high_resolution_clock> t;

		boost::dijkstra_shortest_paths
		(
			graph,
			start_vertex,
			boost::weight_map(get(&GraphEdge::weight, graph)).
			distance_map(boost::make_iterator_property_map(distances.begin(), get(boost::vertex_index, graph))).
			predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(boost::vertex_index, graph))).
			distance_inf(std::numeric_limits<double>::infinity())
		);
		std::cout << "dijkstra_shortest_paths: " << t.seconds() << " s" << std::endl;
	}
	
	{
		timer<high_resolution_clock> t;
		write_shortest_path_tree_vtk(graph, predecessors, distances, inputfilename.filename().replace_extension("_wsp_tree.vtk").string());
		std::cout << "write_shortest_path_tree_vtk: " << t.seconds() << " s" << std::endl;
	}

	if (termination_vertex > 0)
	{
		timer<high_resolution_clock> t;

		write_shortest_path_to_vtk(
			graph, 
			termination_vertex,
			predecessors, 
			distances, 
			inputfilename.filename().replace_extension("_wsp_max.vtk").string()
		);
		std::cout << "write_shortest_path_to_vtk: " << t.seconds() << " s" << std::endl;
	}

	// write_graph_dot("graph.dot", graph);

	std::cout << "This is the end, total time: " << total_time.seconds() << " s" << std::endl;

	return EXIT_SUCCESS;
}
