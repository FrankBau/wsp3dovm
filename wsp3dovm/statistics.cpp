#include "statistics.h"

void dump_mesh(Mesh &mesh)
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

void print_edge_statistics(Mesh &mesh)
{
	double min_lenght = std::numeric_limits<double>::max();
	double max_lenght = 0;

	EdgeHandle min_eh = Kernel::InvalidEdgeHandle;
	EdgeHandle max_eh = Kernel::InvalidEdgeHandle;

	int num_edges = 0;

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		EdgeHandle eh = *e_it;

		double lenght = mesh.length(eh);
		if (lenght < min_lenght)
		{
			min_lenght = lenght;
			min_eh = eh;
		}
		if (lenght > max_lenght)
		{
			max_lenght = lenght;
			max_eh = eh;
		}

		++num_edges;
	}

	std::cout << "min edge length: " << min_lenght << " edge from/to vertex: " << mesh.edge(min_eh) << std::endl;
	std::cout << "max edge length: " << max_lenght << " edge from/to vertex: " << mesh.edge(max_eh) << std::endl;

	const int num_bins = 10;
	int histo[num_bins];
	for (int bin = 0; bin < num_bins; ++bin)
		histo[bin] = 0;

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		double lenght = mesh.length(*e_it);
		int bin = (int)(num_bins * (lenght - min_lenght) / (max_lenght - min_lenght));
		if (bin == num_bins)
			--bin;
		++histo[bin];
	}

	for (int bin = 0; bin < num_bins; ++bin)
	{
		std::cout << "edge histogram bin " << bin << " : " << histo[bin] << std::endl;
	}
}

//    1      |ax bx cx dx|
//V = - * det|ay by cy dy|
//    6      |az bz cz dz|
//           | 1  1  1  1|
double tetrahedral_volume(Point& a, Point& b, Point& c, Point& d )
{
	a -= d;
	b -= d;
	c -= d;
	double det
		= a[0] * b[1] * c[2]
		+ b[0] * c[1] * a[2]
		+ c[0] * a[1] * b[2]
		- c[0] * b[1] * a[2]
		- b[0] * a[1] * c[2]
		- a[0] * c[1] * b[2];

	return fabs(det / 6.0);
}

double tetrahedral_volume(Mesh &mesh, CellHandle ch)
{
	auto v_it = mesh.cv_iter(ch);
	Point p1 = mesh.vertex(*v_it); ++v_it;
	Point p2 = mesh.vertex(*v_it); ++v_it;
	Point p3 = mesh.vertex(*v_it); ++v_it;
	Point p4 = mesh.vertex(*v_it); ++v_it;
	assert(!v_it); // tetrahedron

	double volume =  tetrahedral_volume(p1, p2, p3, p4);

	return volume;
}

void print_volume_statistics(Mesh &mesh)
{
	double min_volume = std::numeric_limits<double>::max();
	double max_volume = std::numeric_limits<double>::min();

	CellHandle min_ch = Kernel::InvalidCellHandle;
	CellHandle max_ch = Kernel::InvalidCellHandle;

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;

		double volume = tetrahedral_volume( mesh, ch );
	
		if (volume < min_volume)
		{
			min_volume = volume;
			min_ch = ch;
		}
		if (volume > max_volume)
		{
			max_volume = volume;
			max_ch = ch;
		}
	}

	std::cout << "min cell volume: " << min_volume << " at cell: " << min_ch.idx() << std::endl;
	std::cout << "max cell volume: " << max_volume << " at cell: " << max_ch.idx() << std::endl;

	const int num_bins = 10;
	int histo[num_bins];
	for (int bin = 0; bin < num_bins; ++bin)
		histo[bin] = 0;

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;

		double volume = tetrahedral_volume(mesh, ch);
		int bin = (int)(num_bins * (volume - min_volume) / (max_volume - min_volume));
		if (bin == num_bins)
			--bin;
		++histo[bin];
	}

	for (int bin = 0; bin < num_bins; ++bin)
	{
		std::cout << "cell volume histogram bin " << bin << " : " << histo[bin] << std::endl;
	}
}

void print_general_statistics(Mesh &mesh)
{
	std::cout << "vertices: " << mesh.n_vertices() << std::endl;
	std::cout << "edges   : " << mesh.n_edges() << std::endl;
	std::cout << "facets  : " << mesh.n_faces() << std::endl;
	std::cout << "cells   : " << mesh.n_cells() << std::endl;
	std::cout << "genus   : " << mesh.genus() << std::endl;
} 

void print_mesh_statistics(Mesh &mesh)
{
	print_general_statistics(mesh);
	print_edge_statistics(mesh);
	print_volume_statistics(mesh);
}