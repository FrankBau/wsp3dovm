#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#if 0

// small test cuboid
int x_min = 0;	// number of cubes in x direction (left to right)
int y_min = 0;		// number of cubes in y direction
int z_min = 0;		// number of cubes in z direction (down)

int x_max = 2;	// number of cubes in x direction (left to right)
int y_max = 3;		// number of cubes in y direction
int z_max = 4;		// number of cubes in z direction (down)

#else

int x_min = -3*8;	// number of cubes in x direction (left to right)
int y_min = -8*3;		// number of cubes in y direction
int z_min = -144/6;		// number of cubes in z direction (down)

int x_max = 3*8;	// number of cubes in x direction (left to right)
int y_max = 8*3;		// number of cubes in y direction
int z_max = 144/6;		// number of cubes in z direction (down)

#endif

// cube edge length
double dx = 1.0/8;
double dy = 1.0/3;
double dz = 1.0*6;

string filename = "cuboid.poly";

// linearized index of a 3D node
int idx(int x, int y, int z)
{
	return (((x-x_min) * (y_max - y_min + 1)) + (y-y_min)) * (z_max - z_min + 1) + (z-z_min) + 1;
}

int weight(int x, int y, int z)
{
	if (z > z_min+(z_max-z_min)/2)
		return 4;
	else
		return 1;
}

int main1(int argc, char * argv)
{

	ofstream file(filename, ios::trunc);
	if (!file.is_open())
	{
		cerr << "failed to open file " << filename << endl;
		return EXIT_FAILURE;
	}

	// Part 1 - node list
	// node count, 3 dim, no attribute, no boundary marker

	file << ((x_max-x_min) + 1)*((y_max-y_min) + 1)*((z_max-z_min) + 1) << " 3 0 0" << endl;

	for (int xx = x_min; xx <= x_max; ++xx)
	{
		for (int yy = y_min; yy <= y_max; ++yy)
		{
			for (int zz = z_min; zz <= z_max; ++zz)
			{
				// Node index, node coordinates
				int node_index = idx(xx, yy, zz);
				file << node_index << " " << dx*xx << " " << dy*yy << " " << dz*zz << endl; // z axis goes down
			}
		}
	}

	// Part 2 - facet list
	// facet count, no boundary marker
	file << ((x_max - x_min) + 1)* (y_max - y_min) * (z_max - z_min) + (x_max - x_min) *((y_max - y_min) + 1)* (z_max - z_min) + (x_max - x_min) * (y_max - y_min) *((z_max - z_min) + 1) << " 0" << endl;

	// facets
	// perpendicular to x axis
	for (int xx = x_min; xx <= x_max; ++xx)
	{
		for (int yy = y_min; yy < y_max; ++yy)
		{
			for (int zz = z_min; zz < z_max; ++zz)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = idx(xx, yy + 0, zz + 0);
				int node2 = idx(xx, yy + 1, zz + 0);
				int node3 = idx(xx, yy + 1, zz + 1);
				int node4 = idx(xx, yy + 0, zz + 1);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// perpendicular to y axis
	for (int xx = x_min; xx < x_max; ++xx)
	{
		for (int yy = y_min; yy <= y_max; ++yy)
		{
			for (int zz = z_min; zz < z_max; ++zz)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = idx(xx + 0, yy, zz + 0);
				int node2 = idx(xx + 1, yy, zz + 0);
				int node3 = idx(xx + 1, yy, zz + 1);
				int node4 = idx(xx + 0, yy, zz + 1);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// perpendicular to z axis
	for (int xx = x_min; xx < x_max; ++xx)
	{
		for (int yy = y_min; yy < y_max; ++yy)
		{
			for (int zz = z_min; zz <= z_max; ++zz)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = idx(xx + 0, yy + 0, zz);
				int node2 = idx(xx + 1, yy + 0, zz);
				int node3 = idx(xx + 1, yy + 1, zz);
				int node4 = idx(xx + 0, yy + 1, zz);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// Part 3 - hole list
	file << "0" << endl; // no holes

	// Part 3 - region list
	// number of regions, each center of a cube defines its own region
	file << (x_max-x_min) * (y_max-y_min) * (z_max-z_min) << endl; // no holes

	for (int xx = x_min; xx < x_max; ++xx)
	{
		for (int yy = y_min; yy < y_max; ++yy)
		{
			for (int zz = z_min; zz < z_max; ++zz)
			{
				int node_index = idx(xx, yy, zz);

				// center coord
				file << node_index << " " << dx * (xx+0.5) << " " << dy * (yy+0.5) << " " << dz * (zz+0.5) << " " << weight(xx,yy,zz) << endl; // z axis goes down
			}
		}
	}

	return EXIT_SUCCESS;
}

/////////////////////////////// generation with open volume mesh ? //////////////////////

#define INCLUDE_TEMPLATES
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshProperty.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>

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

using namespace OpenVolumeMesh;
using namespace OpenVolumeMesh::Geometry;

struct Mesh : GeometricHexahedralMeshV3f
{
	CellPropertyT<float> weight = request_cell_property<float>("weight");

	void write_poly(const std::string filename)
	{
		std::ofstream file(filename, std::ios::trunc);
		if (!file.is_open())
		{
			std::cerr << "failed to open file " << filename << std::endl;
			return;
		}

		file << "# poly file generated by build_cuboid" << endl;

		{
			file << n_vertices() << " 3 0 0" << endl;

			for (auto v : vertices())
			{
				auto point = vertex(v);
				file << std::setw(6) << 1+v.idx() << " " << point[0] << " " << point[1] << " " << point[2] << "\n";
			}
		}

		// select halffaces for output
		StatusAttrib statusAttrib(*this);

		int k = 0;

		for (auto hf : halffaces()) {
			statusAttrib[hf].set_selected(true);
			k++;
		}

		for (auto hf : halffaces()) {
			if (statusAttrib[hf].selected() ) {
				statusAttrib[ opposite_halfface_handle(hf) ].set_selected(false);
				k--;
			}
		}

		{
			file << k << " 0\n";
			for (auto hf : halffaces()) {
				if (statusAttrib[hf].selected()) {
					file << "1" << endl;  // 1 polygon
					file << "4 ";	// we know they are quads
					for (auto v : halfface_vertices(hf)) {
						file << " " << 1 + v.idx();
					}
					file << endl;
				}
			}
		}

		// no holes
		file << " 0" << endl;

		// regions
		{
			file << n_cells() << endl;

			for (auto c : cells()) {
				Vec3f pt = barycenter(c);
				file << c.idx() + 1 << " " << pt[0] << " " << pt[1] << " " << pt[2] << " " << weight[c] << endl;
			}
		}

		file.close();
	}

	void write_nodes(const std::string filename)
	{
		std::ofstream file(filename, std::ios::trunc);
		if (!file.is_open())
		{
			std::cerr << "failed to open file " << filename << std::endl;
			return;
		}

		file << n_vertices() << " 3 0 0\n";
		int i = 1;

		for (auto v : vertices() )
		{
			auto point = vertex(v);
			file << std::setw(6) << i++ << " " << point[0] << " " << point[1] << " " << point[2] << "\n";
		}
		file << "# node file generated by wsp3dovm\n";
		file.close();
	}

	void write_ele(const std::string filename)
	{
		std::ofstream file(filename, std::ios::trunc);
		if (!file.is_open())
		{
			std::cerr << "failed to open file " << filename << std::endl;
			return;
		}

		int i = 1;

		file << n_cells() << "  8  0\n";
		for (auto c : cells()) {
			cout << "n_vertices_in_cell=" << n_vertices_in_cell(c) << endl; // should be 8 
			barycenter(c);
			file << std::setw(5) << i++ << " ";
			for (auto v : cell_vertices(c) )
				file << std::setw(5) << 1 + v.idx() << " "; // tet format is 1-based
			file << "\n";
		}

		file << "# ele file generated by wsp3dovm\n";
		file.close();
	}

	// check for and avoid duplicates
	VertexHandle add_vertex(Vec3f& point) {
		for (auto v : vertices() ) {
			if (point == vertex(v)) {
				return v;
			}
		}
		return OpenVolumeMesh::GeometricHexahedralMeshV3f::add_vertex(point);
	}

	// x increases left -> right
	// y increases top -> bottom
	// z increases into the screen (front --> rear)
	void add_cuboid(int x, int y, int z, int dx, int dy, int dz, float weight=1.0 )
	{
		if (dx < 0) { x += dx; dx = -dx; }
		if (dy < 0) { y += dy; dy = -dy; }
		if (dz < 0) { z += dz; dz = -dz; }

		VertexHandle ltf = add_vertex(Vec3f(x, y, z));
		VertexHandle ltr = add_vertex(Vec3f(x, y, z + dz));
		VertexHandle lbf = add_vertex(Vec3f(x, y + dy, z));
		VertexHandle lbr = add_vertex(Vec3f(x, y + dy, z + dz));
		VertexHandle rtf = add_vertex(Vec3f(x + dx, y, z));
		VertexHandle rtr = add_vertex(Vec3f(x + dx, y, z + dz));
		VertexHandle rbf = add_vertex(Vec3f(x + dx, y + dy, z));
		VertexHandle rbr = add_vertex(Vec3f(x + dx, y + dy, z + dz));

		auto ch = add_cell({
			{ lbf, rbf, rtf, ltf, lbr, ltr, rtr, rbr }
		}, true );
		if (ch == InvalidCellHandle) {
			cerr << "failed to create cuboid, exit." << endl;
			exit(-1);
		}
		this->weight[ch] = weight;
	}
};


/////////////////////////////////////////// d=3 n=2 ////////////////////////////////////
int main1()
{
	Mesh mesh;
	
	int w = 1;

	// center - weight 1
	mesh.add_cuboid( 0, 0, 0, 1,  6,  144, 16 );
	mesh.add_cuboid(-1, 0, 0, 1,  6,  144, 16 );
	mesh.add_cuboid( 0, 0, 0, 1,  6, -144, 16 );
	mesh.add_cuboid(-1, 0, 0, 1,  6, -144, 16 );
	mesh.add_cuboid( 0, 0, 0, 1, -6,  144, 16 );
	mesh.add_cuboid(-1, 0, 0, 1, -6,  144, 16 );
	mesh.add_cuboid( 0, 0, 0, 1, -6, -144, 16 );
	mesh.add_cuboid(-1, 0, 0, 1, -6, -144, 16 );

	// green: top-right and bottom-left
	mesh.add_cuboid( 1, 0, 0,  2, -6,  144, 4 );
	mesh.add_cuboid( 1, 0, 0,  2, -6, -144, 4 );
	mesh.add_cuboid(-1, 0, 0, -2,  6,  144, 4 );
	mesh.add_cuboid(-1, 0, 0, -2,  6, -144, 4 );

	// gold: front-left and rear-right
	mesh.add_cuboid(-1, 8, 0, 1, -2, 144, 1 );
	mesh.add_cuboid( 0, 8, 0, 1, -2, 144, 1 );

	mesh.add_cuboid(1, -8, 0, -1, 2, -144, 1 );
	mesh.add_cuboid(0, -8, 0, -1, 2, -144, 1 );


	cout << " n_cells=" << mesh.n_cells();
	cout << " n_faces=" << mesh.n_faces();
	cout << " n_edges=" << mesh.n_edges();
	cout << " n_vertices=" << mesh.n_vertices();
	cout << " genus=" << mesh.genus();
	cout << endl;

	mesh.write_poly("test.poly");

	// preview (coarse): ..\..\tetgen.exe -pqA test.poly
	// finer: ..\..\tetgen.exe -pqAa0.1 test.poly
	// finest: ..\..\tetgen.exe -pq1.3Aa0.01 test.poly (Mesh tetrahedra: 3810112)

	// generates: .node .ele. face .edge

	// reads .node and .ele, ignores others
	// ..\x64\RelWithDebInfo\wsp3dovm.exe --input-mesh test.1  --write_mesh_vtk 1 --start_vertex 1 --termination_vertex 12

	return 0;
}


int main()
{
	Mesh mesh;

	int n = 4;
	const int wx = 1;
	const int wy = 3 * n;
	const int wz = 4 * 3 * wy * n;

#if 1
	// weight=16: center
	// these are n large slabs of height 1 and weight 1
	// they are only partitioned into smaller parts to match satellite cuboid sizes
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			for (int k=0; k < n; ++k) {
				mesh.add_cuboid( i*wx, j*wy, k*wz, wx, wy, wz, 16 );
			}
		}
	}
#endif

#if 1
	// weight=4: top-right and bottom-left
	for (int j = 0; j < n; ++j) {
		for (int k = 0; k < n; ++k) {
			if (j % 2) {
				mesh.add_cuboid( n, j * wy, k * wz, 1, wy, wz, 4);
			}
			else {
				mesh.add_cuboid( 0, j * wy, k * wz, -1, wy, wz, 4);
			}
		}
	}
#endif

	for (int i = 0; i < n; ++i) {
		for (int k = 0; k < n; ++k) {
			if (k % 2) {
				mesh.add_cuboid( i*wx, n * wy, k * wz, wx, 1, wz, 1);
			}
			else {
				mesh.add_cuboid( i*wx, 0, k * wz, wx, -1, wz, 1);
			}
		}
	}

	// weight=1: front-left and rear-right

#if 0
	mesh.add_cuboid(0, 0, 0, 1, 6, 144, 16);
	mesh.add_cuboid(-1, 0, 0, 1, 6, 144, 16);
	mesh.add_cuboid(0, 0, 0, 1, 6, -144, 16);
	mesh.add_cuboid(-1, 0, 0, 1, 6, -144, 16);
	mesh.add_cuboid(0, 0, 0, 1, -6, 144, 16);
	mesh.add_cuboid(-1, 0, 0, 1, -6, 144, 16);
	mesh.add_cuboid(0, 0, 0, 1, -6, -144, 16);
	mesh.add_cuboid(-1, 0, 0, 1, -6, -144, 16);

	// green: top-right and bottom-left
	mesh.add_cuboid(1, 0, 0, 2, -6, 144, 4);
	mesh.add_cuboid(1, 0, 0, 2, -6, -144, 4);
	mesh.add_cuboid(-1, 0, 0, -2, 6, 144, 4);
	mesh.add_cuboid(-1, 0, 0, -2, 6, -144, 4);

	// gold: front-left and rear-right
	mesh.add_cuboid(-1, 8, 0, 1, -2, 144, 1);
	mesh.add_cuboid(0, 8, 0, 1, -2, 144, 1);

	mesh.add_cuboid(1, -8, 0, -1, 2, -144, 1);
	mesh.add_cuboid(0, -8, 0, -1, 2, -144, 1);
#endif

	cout << " n_cells=" << mesh.n_cells();
	cout << " n_faces=" << mesh.n_faces();
	cout << " n_edges=" << mesh.n_edges();
	cout << " n_vertices=" << mesh.n_vertices();
	cout << " genus=" << mesh.genus();
	cout << endl;

	mesh.write_poly("test.poly");

	// ..\..\tetgen.exe -pqA test.poly (Mesh tetrahedra: 2106337)

	// generates: .node .ele. face .edge

	// reads .node and .ele, ignores others
	// ..\x64\RelWithDebInfo\wsp3dovm.exe --input-mesh test.1  --write_mesh_vtk 1 --start_vertex 60 --termination_vertex 64

	return 0;
}
