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
};

#endif
