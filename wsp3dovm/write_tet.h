#ifndef WRITE_TET_H
#define WRITE_TET_H

// write in tetgen format

#include "common.h"

// write a subset of mesh's cells 
// a .node file and a .ele file will be generated
// the (mesh-)original node numbering is retained
void write_cells_tet
(
	const Mesh &mesh,
	const std::set<CellHandle>& cells,
	const boost::filesystem::path basename
);

#endif
