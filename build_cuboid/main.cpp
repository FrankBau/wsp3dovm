#include <string>
#include <cstdlib>
#include <iostream>
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

int x_min = -3;	// number of cubes in x direction (left to right)
int y_min = -8;		// number of cubes in y direction
int z_min = -144;		// number of cubes in z direction (down)

int x_max = 3;	// number of cubes in x direction (left to right)
int y_max = 8;		// number of cubes in y direction
int z_max = 144;		// number of cubes in z direction (down)

#endif

// cube edge length
double dx = 1.0;
double dy = 1.0;
double dz = 1.0;

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

int main(int argc, char * argv)
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
