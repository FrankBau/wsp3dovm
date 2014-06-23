#include "common.h"

void create_barycentric_steiner_points(Graph &graph, Mesh &mesh);

void create_surface_steiner_points(Graph &graph, Mesh &mesh);

void create_steiner_graph_improved_spanner(Graph &graph, Mesh &mesh, double stretch = 0, double yardstick=0);
