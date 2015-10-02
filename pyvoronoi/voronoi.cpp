#pragma warning(disable : 4503)
#include "voronoi.hpp"

VoronoiDiagram::VoronoiDiagram() {

}

void VoronoiDiagram::AddPoint(Point p) {
	points.push_back(p);
}

void VoronoiDiagram::AddSegment(Segment s) {
	segments.push_back(s);
}

void VoronoiDiagram::Construct() {
	construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);
}

std::vector<c_Edge> VoronoiDiagram::GetEdges() {
	std::vector<c_Edge> edges;
	std::set<const voronoi_edge<double> *> visited;

	for (voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin(); 	it != vd.edges().end(); ++it) {
		if (visited.find(&(*it)) != visited.end())
			continue;
		
		const voronoi_vertex<double> *start = it->vertex0();
		const voronoi_vertex<double> *end = it->vertex1();
		
		c_Vertex startVertex;
		if (start != 0)
			startVertex = c_Vertex(start->x(), start->y());
		
		c_Vertex endVertex;
		if (end != 0)
			endVertex = c_Vertex(end->x(), end->y());
		
		size_t firstIndex = it->cell()->source_index();

		const voronoi_edge<double> *twin = it->twin();
		visited.insert(twin);

		size_t secondIndex = twin->cell()->source_index();

		edges.push_back(c_Edge(start != 0, end != 0, startVertex, endVertex, it->is_primary(), firstIndex, secondIndex));
	}
	
	return edges;
}

std::vector<Point> VoronoiDiagram::GetPoints() {
	return points;
}

std::vector<Segment> VoronoiDiagram::GetSegments() {
	return segments;
}