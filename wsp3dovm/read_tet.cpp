#include "read_tet.h"

#include "MeshGenerator.hpp"

#include <fstream>

// starting with tetgen 1.5.1 beta, there are additional new file formats containing 
// top-down topological information
// We also tried to save (cache) the topology in .ovm file format and usethat later.
// However, even this did not shortcut running time because bottom-up data structures
// have to be generated in memory anyway.

// internally we start at index 0 so we need correction info
// the start value depends on the software, there is no standard
int first_point_number = 0;
int first_edge_number = 0;
int first_facet_number = 0;
int first_tetra_number = 0;

// make nodes from point input, add them to mesh and keep 0-based indices in vector
void read_nodes( MeshGenerator &meshGenerator, std::string filename)
{
	int number_of_nodes;
	int dimension;
	int number_of_attributes;
	bool has_boundary_marker;

	// file format see http://tetgen.berlios.de/fformats.node.html
	std::ifstream node_file(filename);
	if (!node_file.is_open())
	{
		std::cerr << "read_tet: failed to open .node file \"" << filename << "\", exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	node_file >> number_of_nodes >> dimension >> number_of_attributes >> has_boundary_marker;
	if (dimension != 3)
	{
		std::cerr << "read_tet: wrong dimension found in .node file \"" << filename << "\", exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "  reading " << number_of_nodes << " nodes..." << std::endl;
	for (int i = 0; i < number_of_nodes; ++i)
	{
		int j;
		double x;
		double y;
		double z;

		if (!(node_file >> j >> x >> y >> z))
		{
			std::cerr << "read_tet: failed to read node " << i << std::endl;
			exit(EXIT_FAILURE);
		}

		if (i == 0)
		{
			first_point_number = j;
		}
		else
		{
			assert(i + first_point_number == j);
		}

		meshGenerator.add_vertex_component(x);
		meshGenerator.add_vertex_component(y);
		meshGenerator.add_vertex_component(z);

		// attributes, ignored
		for (int k = 0; k < number_of_attributes; ++ k) 
			node_file >> j;

		// boundary_marker, ignored
		if (has_boundary_marker)
			node_file >> j;
	}

	node_file.close();
}

void read_tetras( MeshGenerator &meshGenerator, std::string filename)
{
	int number_of_tetras;
	int number_of_points;
	bool has_boundary_marker;

	// file format see 
	std::ifstream tetra_file(filename);
	if (!tetra_file.is_open())
	{
		std::cerr << "read_tet: failed to open .ele file \"" << filename << "\", exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	tetra_file >> number_of_tetras >> number_of_points >> has_boundary_marker;
	if (number_of_points != 4)
	{
		std::cerr << "number of points 4 expected, found " << number_of_points << ", exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "  reading " << number_of_tetras << " tetras..." << std::endl;
	for (int i = 0; i < number_of_tetras; ++i)
	{
		int j;
		int u;
		int v;
		int w;
		int x;
		if (!(tetra_file >> j >> u >> v >> w >> x))
		{
			std::cerr << "failed to read tetra " << i << std::endl;
			exit(EXIT_FAILURE);
		}
	
		if (i == 0)
		{
			first_tetra_number = j;
		}
		else
		{
			assert(i + first_tetra_number == j);
		}

		u -= first_point_number;
		v -= first_point_number;
		w -= first_point_number;
		x -= first_point_number;

		meshGenerator.add_cell_vertex(u + 1);
		meshGenerator.add_cell_vertex(v + 1);
		meshGenerator.add_cell_vertex(w + 1);
		meshGenerator.add_cell_vertex(x + 1);

		// weight
		double weight;
		if (has_boundary_marker)
		{
			tetra_file >> weight;
		}
		else
		{
			weight = 1.0;
		}
		meshGenerator.mesh()._cellWeight.push_back(weight);
	}
	tetra_file.close();
}

void read_tet(Mesh &mesh, std::string filename)
{
	boost::filesystem::path input_pathname = filename;
	
	std::cout << "read_tet: slurping tetraheder soup from " << filename << " ..." << std::endl;
	MeshGenerator generator(mesh);

	// all nodes
	std::cout << "  read_nodes..." << std::endl;
	read_nodes(generator, input_pathname.string() + ".node");

	// these are only the boundary edges
	// read_edges ( edges, "C:/Carleton/CGAL-4.4/demo/Polyhedron/data/ellipsoid.1.edge");

	// these are only the boundary facets
	// read_facets(facets, "C:/Carleton/CGAL-4.4/demo/Polyhedron/data/ellipsoid.1.face");

	// all tetrahedra

	std::cout << "  enable_bottom_up_incidences(false)..." << std::endl;
	mesh.enable_bottom_up_incidences(false);

	std::cout << "  enable_vertex_bottom_up_incidences(true)..." << std::endl;
	// otherwise add_edge makes an expensive quadratic search
	mesh.enable_vertex_bottom_up_incidences(true);

	std::cout << "  read_tetras..." << std::endl;
	read_tetras(generator, input_pathname.string() + ".ele");

	// now we can free some memory which is needed for other data structures:
	// this caused to performance penalty (tested on NewFineMesh.1 -p)
	// no longer true, needed for vc_iter
	//std::cout << "  enable_vertex_bottom_up_incidences(false)..." << std::endl;
	//mesh.enable_vertex_bottom_up_incidences(false);

	// otherwise would break later calls to some iterators (which?)
	std::cout << "  enable_edge_bottom_up_incidences(true)..." << std::endl;
	mesh.enable_edge_bottom_up_incidences(true);

	// otherwise would break later calls to some iterators (which?)
	std::cout << "  enable_face_bottom_up_incidences(true)..." << std::endl;
	mesh.enable_face_bottom_up_incidences(true);

	std::cout << "read_tet finished." << std::endl;
}
