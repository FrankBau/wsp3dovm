#ifndef COMMON_H
#define COMMON_H

// some VC++2013 internally include stuff needed this
#define _SCL_SECURE_NO_WARNINGS

// prevent min/max define in "Windows Kits\8.1\Include\shared\minwindef.h"
#define NOMINMAX

// C++ includes
#include <iostream>
#include <vector>
#include <map>
#include <tuple>

#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/chrono.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/small_world_generator.hpp>

#include <boost/property_map/property_map.hpp>

// Include vector classes
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshProperty.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>

// Make some typedefs to facilitate your life
typedef double                                      Real;
typedef OpenVolumeMesh::Geometry::Vec3d             Vector;
typedef OpenVolumeMesh::Geometry::Vec3d             Point;

typedef OpenVolumeMesh::VertexHandle            VertexHandle;
typedef OpenVolumeMesh::EdgeHandle              EdgeHandle;
typedef OpenVolumeMesh::FaceHandle              FaceHandle;
typedef OpenVolumeMesh::CellHandle              CellHandle;

typedef OpenVolumeMesh::HalfFaceHandle          HalfFaceHandle;

typedef OpenVolumeMesh::TopologyKernel          Kernel;

typedef Kernel::Cell                            Cell;
typedef Kernel::Face                            Face;

const double max_weight = std::numeric_limits<double>::max();

///////////////////////////// boost Graph /////////////////////////////

struct GraphNode;
struct GraphEdge;

typedef boost::adjacency_list <
	boost::vecS,
	boost::vecS,
	boost::undirectedS,
	GraphNode,
	GraphEdge
> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor GraphNode_descriptor;
typedef boost::graph_traits<Graph>::edge_descriptor   GraphEdge_descriptor;

///////////////////////////// boost Graph details /////////////////////////////

struct GraphNode
{
public:
	GraphNode()
	{
		// std::cout << "GraphNode()" << std::endl;
	}

	CellHandle cell = 0;   // the cell it belongs to
	int i = -1;             // index i (face related node) 
	int j = -1;             // rsp. i,j (edge related node)
	Point point;
};

struct GraphEdge
{
	// we could have done this with internal properties too
	double weight;
};

typedef float Weight;

struct Mesh : public OpenVolumeMesh::GeometricPolyhedralMeshV3d
{
	std::vector<double> _cellWeight;
	std::vector<double> _faceWeight;
	std::vector<double> _edgeWeight;

	double& weight(CellHandle ch) { return _cellWeight[ch.idx()]; }
	double& weight(FaceHandle fh) { return _faceWeight[fh.idx()]; }
	double& weight(EdgeHandle eh) { return _edgeWeight[eh.idx()]; }

	void print_memory_statistics()
	{
		std::cout << "verts: " << n_vertices() << std::endl;
		std::cout << "edges: " << n_edges() << std::endl;
		std::cout << "faces: " << n_faces() << std::endl;
		std::cout << "cells: " << n_cells() << std::endl;

		size_t total_edges_v = 2 * edges_.size();
		std::cout << "total_edges_v:   " << total_edges_v << std::endl;

		size_t total_faces_heh = 0;
		for (int i = 0; i < faces_.size(); ++i)
			total_faces_heh += faces_[i].halfedges().size();
		std::cout << "total_faces_heh: " << total_faces_heh << std::endl;

		size_t total_cells_hfh = 0;
		for (int i = 0; i < cells_.size(); ++i)
			total_cells_hfh += cells_[i].halffaces().size();
		std::cout << "total_cells_hfh: " << total_cells_hfh << std::endl;

		if (this->has_vertex_bottom_up_incidences())
		{
			size_t total_outgoing_hes_per_vertex = 0;
			for (int i = 0; i < outgoing_hes_per_vertex_.size(); ++i)
				total_outgoing_hes_per_vertex += outgoing_hes_per_vertex_[i].size();
			std::cout << "total_outgoing_hes_per_vertex: " << total_outgoing_hes_per_vertex << std::endl;
		}

		if (this->has_edge_bottom_up_incidences())
		{
			size_t total_incident_hfs_per_he = 0;
			for (int i = 0; i < incident_hfs_per_he_.size(); ++i)
				total_incident_hfs_per_he += incident_hfs_per_he_[i].size();
			std::cout << "total_incident_hfs_per_he: " << total_incident_hfs_per_he << std::endl;
		}

		if (this->has_face_bottom_up_incidences())
		{
			size_t total_incident_cell_per_hf = incident_cell_per_hf_.size();
			std::cout << "total_incident_cell_per_hf: " << total_incident_cell_per_hf << std::endl;
		}
	}
};

#endif
