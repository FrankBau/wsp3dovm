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

	int start_vertex;		// for single source shortest paths (Dijkstra)
	int termination_vertex; // an optional termination vertex for which the shortest path will be reported

	int num_start_vertices; // number of start vertices to be used (with incremented indices starting at start_vertex)
	int num_term_vertices;	// number of termination vertices to be used (with incremented indices starting at term_vertex)

	double stretch;		// spaner graph stretch factor	
	double yardstick;	// max. size of edge for edge subdivisions

	program_options::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("start_vertex,s", program_options::value<int>(&start_vertex)->default_value(0), "shortest path start vertex number")
		("num_start_vertices,S", program_options::value<int>(&num_start_vertices)->default_value(1), "number of start vertices to be used (with incremented indices starting at start_vertex)")
		("termination_vertex,t", program_options::value<int>(&termination_vertex)->default_value(-1), "shortest path termination vertex number (-1==none)")
		("num_termination_vertices,T", program_options::value<int>(&num_term_vertices)->default_value(1), "number of termination vertices to be used (with incremented indices starting at term_vertex)")
		("spanner_stretch,x", program_options::value<double>(&stretch)->default_value(0.0), "spanner graph stretch factor")
		("yardstick,y", program_options::value<double>(&yardstick)->default_value(0.0), "interval length for interval scheme (0: do not subdivide edges)")
		("input-mesh", program_options::value<std::string>(), "set input filename (tetgen 3D mesh files wo extension)")
		;

	program_options::positional_options_description positional_options;
	positional_options.add("input-mesh", 1);

	program_options::variables_map vm;
	try
	{
		program_options::store(
			program_options::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), 
			vm
		);
	} 
	catch (...)
	{
		// exceptions in cmd line parsing are aweful, but hapenned easily
		std::cerr << "exception while parsing command line, exit." << std::endl;
		std::cerr << desc << endl;
		return EXIT_FAILURE;
	}

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

		if (stretch < 0.0)
		{
			create_surface_steiner_points(graph, mesh);
			std::cout << "create_surface_steiner_points: " << t.seconds() << " s" << std::endl;
		}
		else
		{
			create_steiner_graph_improved_spanner(graph, mesh, stretch, yardstick );
			std::cout << "create_steiner_graph_improved_spanner with stretch t=: " << stretch << " : " << t.seconds() << " s" << std::endl;
		}
	}
	
	print_steiner_point_statistics(mesh);

	std::cout << "graph nodes: " << graph.m_vertices.size() << std::endl;
	std::cout << "graph edges: " << graph.m_edges.size() << std::endl;

	if(0)
	{
		timer<high_resolution_clock> t;
		write_graph_vtk(graph, inputfilename.filename().replace_extension("_steiner_graph.vtk").string());
		std::cout << "write_graph_vtk: " << t.seconds() << " s" << std::endl;
	}

	std::cout << "running shortest paths calculations for " << num_start_vertices << " start and  " << num_term_vertices << " termination vertices" << std::endl;

	{
		double min_approx_ratio = std::numeric_limits<double>::max();
		double max_approx_ratio = std::numeric_limits<double>::min();
		double sum_approx_ratio = 0;

		timer<high_resolution_clock> t;

		for (int s = start_vertex; s < start_vertex + num_start_vertices; ++s)
		{
			// the distances are temporary, so we choose an external property for that
			std::vector<double> distance(num_vertices(graph));
			std::vector<GraphNode_descriptor> predecessor(num_vertices(graph));

			boost::dijkstra_shortest_paths
				(
				graph,
				s,
				boost::weight_map(get(&GraphEdge::weight, graph)).
				distance_map(boost::make_iterator_property_map(distance.begin(), get(boost::vertex_index, graph))).
				predecessor_map(boost::make_iterator_property_map(predecessor.begin(), get(boost::vertex_index, graph))).
				distance_inf(std::numeric_limits<double>::infinity())
				);

			if (termination_vertex > 0)
			{
				for (int t = termination_vertex; t < termination_vertex + num_term_vertices; ++t)
				{
					double euclidean_distance = norm(graph[s].point, graph[t].point);
					double approx_distance = distance[t];
					double approx_ratio = approx_distance / euclidean_distance;

					//std::cout << "approx. distance from " << s << " to " << t << " is " << 100 * approx_ratio << "% of true euclidean distance" << std::endl;

					min_approx_ratio = min(approx_ratio, min_approx_ratio);
					max_approx_ratio = max(approx_ratio, max_approx_ratio);
					sum_approx_ratio += approx_ratio;
				}
			}
		}

		double avg_approx_ratio = sum_approx_ratio / (num_start_vertices*num_term_vertices);

		std::cout << "total time for " << num_start_vertices << " dijkstra_shortest_paths: " << t.seconds() << " s" << std::endl;
		std::cout << "best  shortest path approximation ratio: " << min_approx_ratio << std::endl;
		std::cout << "avg.  shortest path approximation ratio: " << avg_approx_ratio << std::endl;
		std::cout << "worst shortest path approximation ratio: " << max_approx_ratio << std::endl;
	}
	
#if 0
	{
		timer<high_resolution_clock> t;
		write_shortest_path_tree_vtk(graph, predecessor, distance, inputfilename.filename().replace_extension("_wsp_tree.vtk").string());
		std::cout << "write_shortest_path_tree_vtk: " << t.seconds() << " s" << std::endl;
	}

	if (termination_vertex > 0)
	{
		timer<high_resolution_clock> t;

		write_shortest_path_to_vtk(
			graph, 
			termination_vertex,
			predecessor, 
			distance, 
			inputfilename.filename().replace_extension("_wsp_max.vtk").string()
		);
		std::cout << "write_shortest_path_to_vtk: " << t.seconds() << " s" << std::endl;
	}
#endif

	// write_graph_dot("graph.dot", graph);

	std::cout << "This is the end, total time: " << total_time.seconds() << " s" << std::endl;

	return EXIT_SUCCESS;
}
