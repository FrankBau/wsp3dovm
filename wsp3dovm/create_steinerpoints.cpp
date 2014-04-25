#include "create_steinerpoints.h"

// utmost simple scheme: one point per cell, the barycenter
static void create_steiner_points_for_cell(Graph &graph, std::vector<GraphNode_descriptor> &steiner_points, Mesh &mesh, CellHandle ch)
{
	Point center = mesh.barycenter(ch);

	GraphNode_descriptor node = boost::add_vertex(graph);
	graph[node].cell = ch;
	graph[node].point = center;

	steiner_points[ch.idx()] = node;
}

static void connect_steiner_points_for_cell(Graph &graph, std::vector<GraphNode_descriptor> &steiner_points, Mesh &mesh, CellHandle ch)
{
	GraphNode_descriptor u = steiner_points[ch.idx()];
	Point pu = graph[u].point;

	std::vector<HalfFaceHandle> halffaces = mesh.cell(ch).halffaces();

	for (auto hfh = halffaces.begin(); hfh != halffaces.end(); ++hfh )
	{
		//auto v_it = mesh.hfv_iter(hfh);
		//Point p0 = mesh.vertex(*v_it); ++v_it;
		//Point p1 = mesh.vertex(*v_it); ++v_it;
		//Point p2 = mesh.vertex(*v_it); ++v_it;
		//assert(!v_it); // face is a triangle

		CellHandle och = mesh.incident_cell( mesh.opposite_halfface_handle(*hfh) );

		if (och.is_valid())
		{
			GraphNode_descriptor v = steiner_points[och.idx()];
			Point pv = graph[v].point;

			Weight edge_weight = (pu - pv).norm(); // TODO: weighted length !

			//Weight w0 = sqrt(SquaredDistPoint3Triangle3( pu, p0, p1, p2));
			//Weight w1 = sqrt(SquaredDistPoint3Triangle3( pv, p0, p1, p2));
		
			GraphEdge_descriptor edge;
			bool inserted;
			boost::tie(edge, inserted) = boost::add_edge(u, v, graph);
			assert(inserted);
			graph[edge].weight = edge_weight;
		}
	}
}

void create_steiner_points(Graph &graph, Mesh &mesh)
{
	std::vector<GraphNode_descriptor> steiner_points(mesh.n_cells());

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;
		if (ch.is_valid())
		{
			create_steiner_points_for_cell(graph, steiner_points, mesh, ch);
		}
	}

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;
		if (ch.is_valid())
		{
			connect_steiner_points_for_cell(graph, steiner_points, mesh, ch);
		}
	}
}
