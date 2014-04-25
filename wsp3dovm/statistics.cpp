#include "statistics.h"

void print_mesh_statistics(Mesh &mesh)
{
	double min_lenght = std::numeric_limits<double>::max();
	double max_lenght = 0;

	int num_edges = 0;

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		double lenght = mesh.length(*e_it);
		if (lenght < min_lenght)
			min_lenght = lenght;
		if (lenght > max_lenght)
			max_lenght = lenght;

		++num_edges;
	}

	std::cout << "min edge length: " << min_lenght << std::endl;
	std::cout << "max edge length: " << max_lenght << std::endl;

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
