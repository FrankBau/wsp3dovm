#include "write_tet.h"

// see http://wias-berlin.de/software/tetgen/fformats.html

void write_nodes_tet(const Mesh& mesh, const std::string filename)
{
	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	file << mesh.n_vertices() << " 3 0 0\n";
	int i = 1;

	for (auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
	{
		Point point = mesh.vertex(*vit);
		file << std::setw(6) << i++ << " " << point[0] << " " << point[1] << " " << point[2] << "\n";
	}
	file << "# node file generated by wsp3dovm\n";
	file.close();
}

void write_ele_tet(const Mesh& mesh, const std::set<CellHandle>& cells, const std::string filename)
{
	std::ofstream file(filename, std::ios::trunc);
	if (!file.is_open())
	{
		std::cerr << "failed to open file " << filename << std::endl;
		return;
	}

	int i = 1;

	file << cells.size() << " 4 0\n";
	for (auto ch : cells)
	{
		file << std::setw(6) << i++ << " ";
		for (auto vit = mesh.cv_iter(ch); vit; ++vit)
			file << vit->idx() << " ";
		file << "\n";
	}

	file << "# ele file generated by wsp3dovm\n";
	file.close();
}

void write_cells_tet
(
	const Mesh &mesh,
	const std::set<CellHandle>& cells,
	const boost::filesystem::path basename
)
{
	boost::filesystem::path node_filename = basename.filename().replace_extension(".2.node");
	write_nodes_tet(mesh, node_filename.string());

	boost::filesystem::path ele_filename = basename.filename().replace_extension(".2.ele");
	write_ele_tet(mesh, cells, ele_filename.string());
}