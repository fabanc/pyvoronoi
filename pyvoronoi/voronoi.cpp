#pragma warning(disable : 4503)
#include "voronoi.hpp"
#include "map"

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

std::vector<Point> VoronoiDiagram::GetPoints() {
	return points;
}

std::vector<Segment> VoronoiDiagram::GetSegments() {
	return segments;
}

long long VoronoiDiagram::CountVertices(){
	return vd.num_vertices();
}

long long VoronoiDiagram::CountEdges(){
	return vd.num_edges();
}

long long VoronoiDiagram::CountCells(){
	return vd.num_cells();
}

void VoronoiDiagram::MapVertexIndexes(){
		long long index = 0;
		for (voronoi_diagram<double>::const_vertex_iterator it = vd.vertices().begin(); it != vd.vertices().end(); ++it) {
			const voronoi_diagram<double>::vertex_type* vertex = &(*it);
			map_indexes_to_vertices.insert(index_to_vertex(index, vertex));
			map_vertices_to_indexes.insert(vertex_to_index(vertex, index));
			index++;
		}
}

void VoronoiDiagram::MapEdgeIndexes(){
	long long index = 0;
	for (voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin(); it != vd.edges().end(); ++it) {
		const voronoi_diagram<double>::edge_type* edge = &(*it);
		map_indexes_to_edges.insert(index_to_edge(index, edge));
		map_edges_to_indexes.insert(edge_to_index(edge, index));
		index++;
	}
}

void VoronoiDiagram::MapCellIndexes(){
	long long index = 0;
	for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin(); it != vd.cells().end(); ++it) {
		const voronoi_diagram<double>::cell_type* cell = &(*it);
		map_indexes_to_cells.insert(index_to_cell(index, cell));
		map_cells_to_indexes.insert(cell_to_index(cell, index));
		index++;
	}
}

c_Vertex VoronoiDiagram::GetVertex(long long index){
	const voronoi_diagram<double>::vertex_type* vertex = map_indexes_to_vertices[index];
	double x = vertex->x();
	double y = vertex->y();
	return c_Vertex(x, y);
}

c_Edge VoronoiDiagram::GetEdge(long long index)
{
	const voronoi_diagram<double>::edge_type* edge = map_indexes_to_edges[index];

	//Find vertex references
	long long edge_start = -1;
    long long edge_end = -1;

    if (edge->vertex0() != NULL){
        	edge_start = map_vertices_to_indexes[edge->vertex0()];
    }

    if (edge->vertex1() != NULL){
        	edge_end = map_vertices_to_indexes[edge->vertex1()];
	}

	//Find the twin reference using the segment object
	const voronoi_diagram<double>::edge_type * twin = edge->twin();
	long long twinIndex = -1;
	if (edge != NULL){
		twinIndex = map_edges_to_indexes[twin];
	}

	//Find the cell reference using ther cell object
	long long cellIndex = map_cells_to_indexes[edge->cell()];

	c_Edge c_edge = c_Edge(
		edge_start,
		edge_end,
		edge->is_primary(),
		edge->is_linear(),
		cellIndex,
		twinIndex
	);

	//Return the object
	return c_edge;
}

c_Cell VoronoiDiagram::GetCell(long long index)
{
	//std::map<long long, const voronoi_diagram<double>::cell_type *>::iterator cellMapIterator = cellMap2.find(index);
	const voronoi_diagram<double>::cell_type* cell = map_indexes_to_cells[index];
	std::vector<long long> edge_identifiers;
	std::vector<long long> vertex_identifiers;

	bool is_open = false;

	//Identify the source type
	int source_category = -1;
	if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_SINGLE_POINT){
		source_category = 0;
	}
	else if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_SEGMENT_START_POINT){
		source_category = 1;
	}
	else if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_SEGMENT_END_POINT){
		source_category = 2;
	}
	else if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_INITIAL_SEGMENT){
		source_category = 3;
	}
	else if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_REVERSE_SEGMENT){
		source_category = 4;
	}
	else if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_GEOMETRY_SHIFT){
		source_category = 5;
	}
	else if (cell->source_category() == boost::polygon::SOURCE_CATEGORY_BITMASK){
		source_category = 6;
	}


	const voronoi_diagram<double>::edge_type* edge = cell->incident_edge();
	if (edge != NULL)
	{
		do {
			//Get the edge index
			long long edge_index = map_edges_to_indexes[edge];
			edge_identifiers.push_back(edge_index);

			if (edge->vertex0() == NULL || edge->vertex1() == NULL)
				is_open = true;

			long long edge_start = -1;

			if (edge->vertex0() != NULL){
				edge_start = map_vertices_to_indexes[edge->vertex0()];
			}

			long vertices_count = vertex_identifiers.size();
			if (vertices_count == 0){
				vertex_identifiers.push_back(edge_start);
			}
			else{
				if ( vertex_identifiers.back() != edge_start){
					vertex_identifiers.push_back(edge_start);
				}
			}

			//Move to the next edge
			edge = edge->next();

		} while (edge != cell->incident_edge());
	}

	c_Cell c_cell = c_Cell(
		index,
		cell->source_index(),
		cell->contains_point(),
		cell->contains_segment(),
		is_open,
		source_category
	);

	c_cell.is_degenerate = cell->is_degenerate();
	c_cell.vertices = vertex_identifiers;
	c_cell.edges = edge_identifiers;

	return c_cell;
}
