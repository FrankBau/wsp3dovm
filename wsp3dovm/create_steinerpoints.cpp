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

/////////////////////////////////////////////////////////////////////////////////////////////
// third strategy: ctreate steine nodes and run the improved greedy alg. to calulate a t-spanner
// see http://people.scs.carleton.ca/~michiel/greedyspanner.pdf
// and http://cg.scs.carleton.ca/~mfarshi/pub/ESA05.pdf


// this will add an edge (u,v) only if
// 1 no such edge exists in graph, or
// 2 existing edge is more expensive (because it runs across a face or mesh edge of a more expensive cell)
// in the latter case, the original edge will be replaced by the cheaper one
void add_edge(Graph &graph, GraphNode_descriptor u, GraphNode_descriptor v, Weight cellcost)
{
	Weight weight = cellcost * norm(graph[u].point, graph[v].point);

	// this pair stuff is counter intuitive. first is the edge, second the bool
	std::pair<Graph::edge_descriptor, bool> retrievedEdge = boost::edge(u, v, graph);
	if (retrievedEdge.second)
	{
		graph[retrievedEdge.first].weight = std::min(weight, graph[retrievedEdge.first].weight);
		//std::cout << "updating edge weigth of edge " << u << ", " << v << std::endl;
	}
	else
	{
		std::pair<Graph::edge_descriptor, bool> newEdge = boost::add_edge(u, v, graph);
		graph[newEdge.first].weight = weight;
		//std::cout << "adding new edge " << u << ", " << v << std::endl;
	}
}

void create_steiner_graph_nodes(Graph &graph, Mesh &mesh)
{
	const int nodes_per_vertex = 1;		// 0 or 1
	const int nodes_per_edge = 1;		// 0..k equidistant
	const int nodes_per_face = 1;		// 0..1 (barycenter) ..k (random)
	const int nodes_per_cell = 0;		// only 0 because shortetst paths wont bend within a (constant weight!) cell

	if (nodes_per_vertex == 1)
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
	}

	if (nodes_per_edge > 0)
	{
		mesh._edgeNodes.resize(nodes_per_edge * mesh.n_edges());

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
	}

	if (nodes_per_face > 0)
	{
		mesh._faceNodes.resize(nodes_per_face * mesh.n_faces());
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
	}
}

std::vector<GraphNode_descriptor> cell_nodes(Graph &graph, Mesh &mesh, CellHandle ch)
{
	std::set<EdgeHandle> edges;			// edges of cell
	std::set<VertexHandle> vertices;	// vertices of cell

	std::vector<GraphNode_descriptor> all_nodes;
	// collect all graph nodes belonging to that cell

	for (auto hfh : mesh.cell(ch).halffaces())
	{
		FaceHandle fh = mesh.face_handle(hfh);
		all_nodes.insert(all_nodes.end(), mesh.f_nodes(fh).begin(), mesh.f_nodes(fh).end());

		Face face = mesh.face(fh);

		for (auto heh : face.halfedges())
		{
			EdgeHandle eh = mesh.edge_handle(heh);
			if (edges.find(eh) == edges.end())
			{
				all_nodes.insert(all_nodes.end(), mesh.e_nodes(eh).begin(), mesh.e_nodes(eh).end());
				edges.insert(eh);

				Edge edge = mesh.edge(eh);

				VertexHandle vh1 = edge.from_vertex();
				if (vertices.find(vh1) == vertices.end())
				{
					all_nodes.insert(all_nodes.end(), mesh.v_node(vh1));
					vertices.insert(vh1);
				}

				VertexHandle vh2 = edge.to_vertex();
				if (vertices.find(vh2) == vertices.end())
				{
					all_nodes.insert(all_nodes.end(), mesh.v_node(vh2));
					vertices.insert(vh2);
				}
			}
		}
	}
	return all_nodes;
}

// after fighting with boost::subgraph for a while, I decided to re-invent the wheel

struct SpannerGraphNode;
struct SpannerGraphEdge;

typedef
boost::adjacency_list
<
	boost::vecS,
	boost::vecS,
	boost::undirectedS,
	SpannerGraphNode,
	SpannerGraphEdge
>
SpannerGraph;

typedef boost::graph_traits<SpannerGraph>::vertex_descriptor SpannerGraphNode_descriptor;
typedef boost::graph_traits<SpannerGraph>::edge_descriptor   SpannerGraphEdge_descriptor;

struct SpannerGraphNode
{
	GraphNode_descriptor original_node;
};

struct SpannerGraphEdge
{
	SpannerGraphEdge() : length(-1) {}

	SpannerGraphEdge(SpannerGraphNode_descriptor u, SpannerGraphNode_descriptor v, double length)
	: u(u), v(v), length(length)
	{
	}

	bool operator < (const SpannerGraphEdge &rhs) const { return length < rhs.length; }

	SpannerGraphNode_descriptor u;
	SpannerGraphNode_descriptor v;
	double length;
};

void create_steiner_graph_improved_spanner(Graph &graph, Mesh &mesh, double stretch)
{
	create_steiner_graph_nodes(graph, mesh);

	for (auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it)
	{
		CellHandle ch = *it;

		std::vector<GraphNode_descriptor> nodes = cell_nodes(graph, mesh, ch);

		SpannerGraph spanner;

		for (auto node : nodes)
		{
			SpannerGraphNode_descriptor u = boost::add_vertex(spanner);
			spanner[u].original_node = node;
		}

		// now we build a sorted list of potential edges
		std::vector<SpannerGraphEdge> potential_edges;

		auto uit = vertices(spanner);
		for (auto u = uit.first; u != uit.second; ++u)
		{
			auto vit = vertices(spanner);
			for (auto v = u+1; v != vit.second; ++v)
			{
				double length = norm(graph[spanner[*u].original_node].point, graph[spanner[*v].original_node].point);
				assert(length > 0);
				potential_edges.push_back(SpannerGraphEdge(*u, *v, length));
			}
		}

		// sort edges by length
		std::sort(potential_edges.begin(), potential_edges.end());

		// determine shortest path arleady in subgraph
		for (auto potential_edge : potential_edges)
		{
			std::vector<double> distances(num_vertices(spanner));
			std::vector<SpannerGraphNode_descriptor> predecessors(num_vertices(spanner));

			SpannerGraphNode_descriptor u = potential_edge.u;
			SpannerGraphNode_descriptor v = potential_edge.v;

			// this is the "unimproved" version with cubic runtime
			boost::dijkstra_shortest_paths
			(
				spanner,
				u,
				boost::weight_map(get(&SpannerGraphEdge::length, spanner)).
				distance_map(boost::make_iterator_property_map(distances.begin(), get(boost::vertex_index, graph))).
				predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(boost::vertex_index, graph))).
				distance_inf(std::numeric_limits<double>::infinity())
			);

			GraphNode_descriptor ou = spanner[u].original_node;
			GraphNode_descriptor ov = spanner[v].original_node;

			Point pu = graph[ou].point;
			Point pv = graph[ov].point;

			// add potential edge to graph iff its shorter that limit
			if (distances[v] > (1 + stretch)*norm(pu,pv))
			{
				// add edge to spanner and to big graph
				SpannerGraphEdge_descriptor spanner_edge = boost::add_edge(u, v, spanner ).first;
				spanner[spanner_edge].length = potential_edge.length;
				spanner[spanner_edge].u = potential_edge.u;
				spanner[spanner_edge].v = potential_edge.v;

				add_edge(graph, ou, ov, mesh.weight(ch) );
			}
		}
	}
}
