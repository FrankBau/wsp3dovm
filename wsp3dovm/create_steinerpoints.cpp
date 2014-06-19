#include "create_steinerpoints.h"

// utmost simple scheme: one point per cell, the barycenter
static void create_barycentric_steiner_points_for_cell(Graph &graph, std::vector<GraphNode_descriptor> &steiner_points, Mesh &mesh, CellHandle ch)
{
	Point center = mesh.barycenter(ch);

	GraphNode_descriptor node = boost::add_vertex(graph);
	graph[node].point = center;

	steiner_points[ch.idx()] = node;
}

static void connect_barycentric_steiner_points_for_cell(Graph &graph, std::vector<GraphNode_descriptor> &steiner_points, Mesh &mesh, CellHandle ch)
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

			Weight edge_weight = static_cast<Weight>((pu - pv).norm()); // TODO: weighted length !

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

void create_barycentric_steiner_points(Graph &graph, Mesh &mesh)
{
	std::vector<GraphNode_descriptor> steiner_points(mesh.n_cells());

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;
		if (ch.is_valid())
		{
			create_barycentric_steiner_points_for_cell(graph, steiner_points, mesh, ch);
		}
	}

	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
	{
		CellHandle ch = *c_it;
		if (ch.is_valid())
		{
			connect_barycentric_steiner_points_for_cell(graph, steiner_points, mesh, ch);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
// second strategy: use all vertices as graph nodes, create steiner points for faces and edges

void connect_surface_steiner_points_for_cell(Graph &graph, Mesh& mesh, CellHandle ch)
{
	for (auto hfh1 : mesh.cell(ch).halffaces())
	{
		FaceHandle fh1 = mesh.face_handle(hfh1);

		for (auto hfh2 : mesh.cell(ch).halffaces())
		{
			FaceHandle fh2 = mesh.face_handle(hfh2);

			if (fh1 < fh2) // do respect (an arbitrary) face order to avoid multiple edges
			{
				// connect every fh1 node to every fh2 node
				for (auto f1_node : mesh.f_nodes(fh1))
				{
					for (auto f2_node : mesh.f_nodes(fh2))
					{
						GraphEdge_descriptor edge;
						bool inserted;
						boost::tie(edge, inserted) = boost::add_edge(f1_node, f2_node, graph);
						double length = norm(graph[f1_node].point, graph[f2_node].point);
						graph[edge].weight = length * mesh.weight(ch);
						assert(inserted);
					}
				}
			}
		}
	}
}

void connect_surface_steiner_points_for_face(Graph &graph, Mesh& mesh, FaceHandle fh)
{
	// incident edge-edge connections
	for (auto heh1 : mesh.face(fh).halfedges())
	{
		EdgeHandle eh1 = mesh.edge_handle(heh1);

		for (auto heh2 : mesh.face(fh).halfedges())
		{
			EdgeHandle eh2 = mesh.edge_handle(heh2);

			if (eh1 < eh2) // enforce arbitrary order to avoid duplicates
			{
				// connect every fh1 node to every fh2 node
				for (auto e1_node : mesh.e_nodes(eh1))
				{
					for (auto e2_node : mesh.e_nodes(eh2))
					{
						GraphEdge_descriptor edge;
						bool inserted;
						boost::tie(edge, inserted) = boost::add_edge(e1_node, e2_node, graph);
						double length = norm(graph[e1_node].point, graph[e2_node].point );
						graph[edge].weight = length * mesh.weight(fh);
						assert(inserted);
					}
				}
			}
		}
	}

	// face interior node(s) to incident edges
	for (auto heh : mesh.face(fh).halfedges())
	{
		EdgeHandle eh = mesh.edge_handle(heh);

		for (auto e_node : mesh.e_nodes(eh))
		{
			for (auto f_node : mesh.f_nodes(fh))
			{
				GraphEdge_descriptor edge;
				bool inserted;
				boost::tie(edge, inserted) = boost::add_edge( e_node, f_node, graph);
				double length = norm( graph[e_node].point, graph[f_node].point );
				graph[edge].weight = length * mesh.weight(fh);
				assert(inserted);
			}
		}
	}
}

void connect_surface_steiner_points_for_edge(Graph& graph, Mesh& mesh, EdgeHandle eh)
{
	Edge e = mesh.edge(eh);
	VertexHandle v1 = e.from_vertex();
	VertexHandle v2 = e.to_vertex();
	GraphNode_descriptor node1 = mesh.v_node(v1);
	GraphNode_descriptor node2 = mesh.v_node(v2);

	// this is correct fo 1 edge node  only, otherwise, all nodes should be sorted
	assert(mesh.e_nodes(eh).size()==1);
	for (auto node : mesh.e_nodes(eh))
	{
		GraphEdge_descriptor edge1;
		GraphEdge_descriptor edge2;

		bool inserted;
		boost::tie(edge1, inserted) = boost::add_edge(node1, node, graph);
		assert(inserted);
		boost::tie(edge2, inserted) = boost::add_edge(node, node2, graph);
		assert(inserted);
		
		double length1 = norm(graph[node1].point, graph[node].point);
		double length2 = norm(graph[node].point, graph[node2].point);

		graph[edge1].weight = length1 * mesh.weight(eh);
		graph[edge2].weight = length2 * mesh.weight(eh);
	}
}

void create_surface_steiner_points(Graph &graph, Mesh &mesh)
{
	mesh._vertexNode.resize(mesh.n_vertices());

	// create a graph node for each mesh vertex
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it)
	{
		VertexHandle vh = *it;
		if (vh.is_valid())
		{
			GraphNode_descriptor node = boost::add_vertex(graph);
			graph[node].vertex = vh;
			graph[node].point = mesh.vertex(vh);
			mesh.v_node(vh) = node;
		}
	}

	mesh._edgeNodes.resize(mesh.n_edges());

	// create a graph node for each mesh edge
	for (auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it)
	{
		EdgeHandle eh = *it;
		if (eh.is_valid())
		{
			GraphNode_descriptor node = boost::add_vertex(graph);
			graph[node].point = mesh.barycenter(eh);
			mesh.e_nodes(eh).push_back(node);
		}
	}

	mesh._faceNodes.resize(mesh.n_faces());
	// create a graph node for each mesh face
	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it)
	{
		FaceHandle fh = *it;
		if (fh.is_valid())
		{
			GraphNode_descriptor node = boost::add_vertex(graph);
			graph[node].point = mesh.barycenter(fh);
			mesh.f_nodes(fh).push_back(node);
		}
	}

	// we do not add cell interior graph nodes because shortest paths wont bend in the interior of a cell

	// connect graph nodes within edges
	for (auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it)
	{
		EdgeHandle eh = *it;
		if (eh.is_valid())
		{
			connect_surface_steiner_points_for_edge(graph, mesh, eh);
		}
	}

	// connect graph nodes within faces
	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it)
	{
		FaceHandle fh = *it;
		if (fh.is_valid())
		{
			connect_surface_steiner_points_for_face(graph, mesh, fh);
		}
	}

	// connect graph nodes between cell faces 
	for (auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it)
	{
		CellHandle ch = *it;
		if (ch.is_valid())
		{
			connect_surface_steiner_points_for_cell(graph, mesh, ch);
		}
	}
}
