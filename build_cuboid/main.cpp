#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char * argv)
{
	int x = 200;	// number of cubes in x direction (left to right)
	int y = 1;	// number of cubes in y direction
	int z = 50;	// number of cubes in z direction (down)

	// cube edge length
	double dx = 1.0;
	double dy = 1.0;
	double dz = 1.0;

	string filename = "cuboid.poly";

	ofstream file(filename, ios::trunc);
	if (!file.is_open())
	{
		cerr << "failed to open file " << filename << endl;
		return EXIT_FAILURE;
	}

	// Part 1 - node list
	// node count, 3 dim, no attribute, no boundary marker

	file << (x + 1)*(y + 1)*(z + 1) << " 3 0 0" << endl;

	for (int i = 0; i <= x; ++i)
	{
		for (int j = 0; j <= y; ++j)
		{
			for (int k = 0; k <= z; ++k)
			{
				// Node index, node coordinates
				int node_index = 1 + i*(y + 1)*(z + 1) + j*(z + 1) + k;
				file << node_index << " " << dx*i << " " << dy*j << " " << dz*k << endl;
			}
		}
	}

	// Part 2 - facet list
	// facet count, no boundary marker
	file << (x + 1)*y*z + x*(y + 1)*z + x*y*(z + 1) << " 0" << endl;

	// facets
	// perpendicular to x axis
	for (int i = 0; i <= x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			for (int k = 0; k < z; ++k)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = 1 + i*(y + 1)*(z + 1) + (j + 0)*(z + 1) + (k + 0);
				int node2 = 1 + i*(y + 1)*(z + 1) + (j + 1)*(z + 1) + (k + 0);
				int node3 = 1 + i*(y + 1)*(z + 1) + (j + 1)*(z + 1) + (k + 1);
				int node4 = 1 + i*(y + 1)*(z + 1) + (j + 0)*(z + 1) + (k + 1);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// perpendicular to y axis
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j <= y; ++j)
		{
			for (int k = 0; k < z; ++k)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = 1 + (i + 0)*(y + 1)*(z + 1) + j*(z + 1) + (k + 0);
				int node2 = 1 + (i + 1)*(y + 1)*(z + 1) + j*(z + 1) + (k + 0);
				int node3 = 1 + (i + 1)*(y + 1)*(z + 1) + j*(z + 1) + (k + 1);
				int node4 = 1 + (i + 0)*(y + 1)*(z + 1) + j*(z + 1) + (k + 1);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// perpendicular to z axis
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			for (int k = 0; k <= z; ++k)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = 1 + (i + 0)*(y + 1)*(z + 1) + (j + 0)*(z + 1) + k;
				int node2 = 1 + (i + 1)*(y + 1)*(z + 1) + (j + 0)*(z + 1) + k;
				int node3 = 1 + (i + 1)*(y + 1)*(z + 1) + (j + 1)*(z + 1) + k;
				int node4 = 1 + (i + 0)*(y + 1)*(z + 1) + (j + 1)*(z + 1) + k;

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// Part 3 - hole list
	file << "0" << endl; // no holes

	// Part 3 - region list
	file << "0" << endl; // no regions (yet)

	return EXIT_SUCCESS;
}
