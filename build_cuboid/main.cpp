#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char * argv)
{
	int x = 100;	// number of cubes in x direction (left to right)
	int y = 1;		// number of cubes in y direction
	int z = 31;		// number of cubes in z direction (down)

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

	for (int xx = 0; xx <= x; ++xx)
	{
		for (int yy = 0; yy <= y; ++yy)
		{
			for (int zz = 0; zz <= z; ++zz)
			{
				// Node index, node coordinates
				int node_index = 1 + xx*(y + 1)*(z + 1) + yy*(z + 1) + zz;
				file << node_index << " " << dx*xx << " " << dy*yy << " " << -(dz*zz) << endl; // z axis goes down
			}
		}
	}

	// Part 2 - facet list
	// facet count, no boundary marker
	file << (x + 1)*y*z + x*(y + 1)*z + x*y*(z + 1) << " 0" << endl;

	// facets
	// perpendicular to x axis
	for (int xx = 0; xx <= x; ++xx)
	{
		for (int yy = 0; yy < y; ++yy)
		{
			for (int zz = 0; zz < z; ++zz)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = 1 + xx*(y + 1)*(z + 1) + (yy + 0)*(z + 1) + (zz + 0);
				int node2 = 1 + xx*(y + 1)*(z + 1) + (yy + 1)*(z + 1) + (zz + 0);
				int node3 = 1 + xx*(y + 1)*(z + 1) + (yy + 1)*(z + 1) + (zz + 1);
				int node4 = 1 + xx*(y + 1)*(z + 1) + (yy + 0)*(z + 1) + (zz + 1);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// perpendicular to y axis
	for (int xx = 0; xx < x; ++xx)
	{
		for (int yy = 0; yy <= y; ++yy)
		{
			for (int zz = 0; zz < z; ++zz)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = 1 + (xx + 0)*(y + 1)*(z + 1) + yy*(z + 1) + (zz + 0);
				int node2 = 1 + (xx + 1)*(y + 1)*(z + 1) + yy*(z + 1) + (zz + 0);
				int node3 = 1 + (xx + 1)*(y + 1)*(z + 1) + yy*(z + 1) + (zz + 1);
				int node4 = 1 + (xx + 0)*(y + 1)*(z + 1) + yy*(z + 1) + (zz + 1);

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// perpendicular to z axis
	for (int xx = 0; xx < x; ++xx)
	{
		for (int yy = 0; yy < y; ++yy)
		{
			for (int zz = 0; zz <= z; ++zz)
			{
				// 1 polygon, no hole, no boundary marker
				file << "1" << endl;

				int node1 = 1 + (xx + 0)*(y + 1)*(z + 1) + (yy + 0)*(z + 1) + zz;
				int node2 = 1 + (xx + 1)*(y + 1)*(z + 1) + (yy + 0)*(z + 1) + zz;
				int node3 = 1 + (xx + 1)*(y + 1)*(z + 1) + (yy + 1)*(z + 1) + zz;
				int node4 = 1 + (xx + 0)*(y + 1)*(z + 1) + (yy + 1)*(z + 1) + zz;

				file << "4 " << node1 << " " << node2 << " " << node3 << " " << node4 << endl;
			}
		}
	}

	// Part 3 - hole list
	file << "0" << endl; // no holes

	// Part 3 - region list
#if 0
	// number of regions
	file << "0" << endl; // no regions (yet)
#elif defined(ZURICH_1)
	// number of regions
	file << x*y*z << endl;

	for (int xx = 0; xx < x; ++xx)
	{
		for (int yy = 0; yy < y; ++yy)
		{
			for (int zz = 0; zz < z; ++zz)
			{
				int region_index = 1 + xx*y*z + yy*z + zz;
				
				unsigned int weight; // tetgen conversion needs integral values
				if (zz >= 30)
					weight = 10000 / 5000;
				else if (zz >= 10)
					weight = 10000 / 2000;
				else
					weight = 10000 / 500;

				file << region_index << " " << dx*xx+dx/2 << " " << dy*yy+dy/2 << " " << -(dz*zz+dz/2) << " " << weight << endl;
			}
		}
	}
#else // ZURICH_2
	// number of regions
	file << x*y*z << endl;

	for (int xx = 0; xx < x; ++xx)
	{
		for (int yy = 0; yy < y; ++yy)
		{
			for (int zz = 0; zz < z; ++zz)
			{
				int region_index = 1 + xx*y*z + yy*z + zz;

				unsigned int weight; // tetgen conversion needs integral values
				if (zz >= 30)
					weight = 10000 / (5000 + 200*(zz-30));	// v3 = 5000 m/s (dv/dz = 200 m/s)
				else if (zz >= 10)
					weight = 10000 / (2000 + 100*(zz-10));	// v2 = 2000 m/s (dv/dz = 100 m/s),
				else
					weight = 10000 / ( 500 + 100*(zz- 0));	// v1 = 500 m/s (dv/dz = 100 m/s),

				file << region_index << " " << dx*xx + dx / 2 << " " << dy*yy + dy / 2 << " " << -(dz*zz + dz / 2) << " " << weight << endl;
			}
		}
	}
#endif

	return EXIT_SUCCESS;
}
