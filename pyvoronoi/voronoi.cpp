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

void VoronoiDiagram::GetEdges(std::vector<c_Vertex> &vertices, std::vector<c_Edge> &edges) {
	std::set<const voronoi_edge<double> *> visited;
	std::map<const voronoi_vertex<double> *, long long> vertexMap;

	for (voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin(); 	it != vd.edges().end(); ++it) {
		if (visited.find(&(*it)) != visited.end())
			continue;
		
		const voronoi_vertex<double> *start = it->vertex0();
		const voronoi_vertex<double> *end = it->vertex1();
		
		long long startIndex = -1;
		if (start != 0) {
			std::map<const voronoi_vertex<double> *, long long>::iterator it = vertexMap.find(start);

			if (it == vertexMap.end())
			{
				c_Vertex startVertex = c_Vertex(start->x(), start->y());
				startIndex = (long long) vertices.size();
				vertices.push_back(startVertex);
				vertexMap[start] = startIndex;
			}
			else {
				startIndex = it->second;
			}
		}
		
		long long endIndex = -1;
		if (end != 0) {
			std::map<const voronoi_vertex<double> *, long long>::iterator it = vertexMap.find(end);

			if (it == vertexMap.end())
			{
				c_Vertex endVertex = c_Vertex(end->x(), end->y());
				endIndex = (long long) vertices.size();
				vertices.push_back(endVertex);
				vertexMap[end] = endIndex;
			}
			else {
				endIndex = it->second;
			}
		}
		
		size_t firstIndex = it->cell()->source_index();

		const voronoi_edge<double> *twin = it->twin();
		visited.insert(twin);

		size_t secondIndex = twin->cell()->source_index();

		edges.push_back(c_Edge(startIndex, endIndex, it->is_primary(), firstIndex, secondIndex));
	}
}

std::vector<Point> VoronoiDiagram::GetPoints() {
	return points;
}

std::vector<Segment> VoronoiDiagram::GetSegments() {
	return segments;
}