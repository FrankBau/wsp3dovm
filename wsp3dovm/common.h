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
#include <set>
#include <tuple>

#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/chrono.hpp>
#include <boost/program_options.hpp>
// I had to remove libboost_program_options... from the linker input list to avoid duplicate symbols

#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/small_world_generator.hpp>

#include <boost/property_map/property_map.hpp>

#define INCLUDE_TEMPLATES
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshProperty.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#undef INCLUDE_TEMPLATES

// Make some typedefs to facilitate your life
typedef double                                      Real;
typedef OpenVolumeMesh::Geometry::Vec3d             Vector;
typedef OpenVolumeMesh::Geometry::Vec3d             Point;

typedef OpenVolumeMesh::VertexHandle            VertexHandle;
typedef OpenVolumeMesh::EdgeHandle              EdgeHandle;
typedef OpenVolumeMesh::FaceHandle              FaceHandle;
typedef OpenVolumeMesh::CellHandle              CellHandle;
typedef OpenVolumeMesh::OpenVolumeMeshHandle    MeshHandle;

typedef OpenVolumeMesh::HalfFaceHandle          HalfFaceHandle;
typedef OpenVolumeMesh::HalfEdgeHandle          HalfEdgeHandle;

typedef OpenVolumeMesh::TopologyKernel          Kernel;

typedef Kernel::Cell                            Cell;
typedef Kernel::Face                            Face;
typedef Kernel::Edge                            Edge;

typedef double Weight;

const Weight max_weight = std::numeric_limits<double>::max();

const double epsilon = 1E-8; // smaller differences are to be considered 0

///////////////////////////// boost Graph /////////////////////////////

struct GraphNode;
struct GraphEdge;

typedef 
		boost::adjacency_list 
		<
			boost::vecS,
			boost::vecS,
			boost::undirectedS,
			GraphNode,
			GraphEdge
		>
Graph;

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

	// the original vertex in tetrahedralization
	VertexHandle vertex = OpenVolumeMesh::TopologyKernel::InvalidCellHandle;
	Point point;	// geometric location
};

struct GraphEdge
{
	// we could have done this with internal properties too
	Weight weight;
};

struct Mesh : public OpenVolumeMesh::GeometricPolyhedralMeshV3d
{
	std::vector<double> _cellWeight;
	std::vector<double> _faceWeight;
	std::vector<double> _edgeWeight;
	
	// graph nodes interior to mesh features
	// we dont have cell interior nodes for now
	std::vector<std::vector<GraphNode_descriptor>> _faceNodes;
	std::vector<std::vector<GraphNode_descriptor>> _edgeNodes;
	std::vector<GraphNode_descriptor> _vertexNode; // only one per vertex

	double& weight(CellHandle ch) { return _cellWeight[ch.idx()]; }
	double& weight(FaceHandle fh) { return _faceWeight[fh.idx()]; }
	double& weight(EdgeHandle eh) { return _edgeWeight[eh.idx()]; }

	std::vector<GraphNode_descriptor>& f_nodes(FaceHandle fh) { return _faceNodes[fh.idx()]; }
	std::vector<GraphNode_descriptor>& e_nodes(EdgeHandle eh) { return _edgeNodes[eh.idx()]; }
	GraphNode_descriptor& v_node(VertexHandle vh) { return _vertexNode[vh.idx()]; }

	void print_memory_statistics()
	{
		std::cout << "verts: " << n_vertices() << std::endl;
		std::cout << "edges: " << n_edges() << std::endl;
		std::cout << "faces: " << n_faces() << std::endl;
		std::cout << "cells: " << n_cells() << std::endl;

		size_t total_edges_v = 2 * edges_.size();
		std::cout << "total_edges_v:   " << total_edges_v << std::endl;

		size_t total_faces_heh = 0;
		for (size_t i = 0; i < faces_.size(); ++i)
			total_faces_heh += faces_[i].halfedges().size();
		std::cout << "total_faces_heh: " << total_faces_heh << std::endl;

		size_t total_cells_hfh = 0;
		for (size_t i = 0; i < cells_.size(); ++i)
			total_cells_hfh += cells_[i].halffaces().size();
		std::cout << "total_cells_hfh: " << total_cells_hfh << std::endl;

		if (this->has_vertex_bottom_up_incidences())
		{
			size_t total_outgoing_hes_per_vertex = 0;
			for (size_t i = 0; i < outgoing_hes_per_vertex_.size(); ++i)
				total_outgoing_hes_per_vertex += outgoing_hes_per_vertex_[i].size();
			std::cout << "total_outgoing_hes_per_vertex: " << total_outgoing_hes_per_vertex << std::endl;
		}

		if (this->has_edge_bottom_up_incidences())
		{
			size_t total_incident_hfs_per_he = 0;
			for (size_t i = 0; i < incident_hfs_per_he_.size(); ++i)
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

static inline double norm(const Vector& v)
{
	return sqrt(v.sqrnorm());
}

static inline double norm(const Point& p, const Point& q)
{
	return norm(p-q);
}

static inline double length(const Mesh& mesh, const Edge& e)
{
	// hu, ugly that we need the mesh here...
	return norm(mesh.vertex(e.from_vertex()), mesh.vertex(e.to_vertex()));
}

#endif
