#include "read_tet.h"

#include "MeshGenerator.hpp"

#include <fstream>

// internally we start at index 0 so we need correction info
// the start value depends on the software, there is no standard
int first_point_number = 0;
int first_edge_number = 0;
int first_facet_number = 0;
int first_tetra_number = 0;

// make vertices from point input, add them to mesh and keep 0-based indices in vector
void read_vertices( MeshGenerator &meshGenerator, std::string filename)
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
	std::cout << "reading " << number_of_nodes << " nodes...";
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

	std::cout << " done." << std::endl;
	node_file.close();
}

void read_tetras( MeshGenerator &meshGenerator, std::string filename)
{
	int number_of_tetras;
	int number_of_points;
	bool has_boundary_markers;

	// file format see 
	std::ifstream file(filename);
	file >> number_of_tetras >> number_of_points >> has_boundary_markers;
	if (number_of_points != 4)
	{
		std::cerr << "number of points 4 expected, found " << number_of_points << ", exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "reading " << number_of_tetras << " tetras...";
	for (int i = 0; i < number_of_tetras; ++i)
	{
		int j;
		int u;
		int v;
		int w;
		int x;
		if (!(file >> j >> u >> v >> w >> x))
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

		// ignored
		if (has_boundary_markers)
			file >> j;
	}
	std::cout << " done." << std::endl;
	file.close();
}


void read_tet(Mesh &mesh, std::string filename)
{
	boost::filesystem::path input_pathname = filename;
	
	std::cout << "slurping tetraheder soup from " << filename << " ..." << std::endl;
	MeshGenerator generator(mesh);

	// all nodes
	read_vertices(generator, input_pathname.string() + ".node");

	// these are only the boundary edges
	// read_edges ( edges, "C:/Carleton/CGAL-4.4/demo/Polyhedron/data/ellipsoid.1.edge");

	// these are only the boundary facets
	// read_facets(facets, "C:/Carleton/CGAL-4.4/demo/Polyhedron/data/ellipsoid.1.face");

	// all tetrahedra

	mesh.enable_bottom_up_incidences(false);
	mesh.enable_vertex_bottom_up_incidences(true); // otherwise add_edge makes an expensive quadratic search
	read_tetras(generator, input_pathname.string() + ".ele");
	mesh.enable_bottom_up_incidences(true); // rebuild them at once is much cheaper

	std::cout << "read_tet finished." << std::endl;
}
