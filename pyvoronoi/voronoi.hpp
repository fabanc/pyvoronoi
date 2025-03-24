// Voronoi.cpp : Defines the entry point for the console application.
//
#define BOOST_POLYGON_NO_DEPS
#define BOOST_NO_USER_CONFIG
#define BOOST_NO_COMPILER_CONFIG
#define BOOST_NO_STDLIB_CONFIG
#define BOOST_NO_PLATFORM_CONFIG
#define BOOST_HAS_STDINT_H

#undef __GLIBC__

#include "boost/polygon/voronoi.hpp"
#include "map"

struct Point {
	int X;
	int Y;
	Point(int x = 0, int y = 0) : X(x), Y(y) {}
};

struct Segment {
	Point p0;
	Point p1;
	Segment(Point a = Point(), Point b = Point()) : p0(a.X, a.Y), p1(b.X, b.Y) {}

    //https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

    // Given three collinear points p, q, r, the function checks if
    // point q lies on line segment 'pr'
    bool onSegment(Point p, Point q, Point r)
    {
        // If the segments just equals an endpoint, returns false. It touches, but does not intersect
        if (q.X == p.X && q.Y == p.Y)
            return false;
        if (q.X == r.X && q.Y == r.Y)
            return false;
        // Otherwise, it is in the pr space, returns true
        if (q.X <= std::max(p.X, r.X) && q.X >= std::min(p.X, r.X) && q.Y <= std::max(p.Y, r.Y) && q.Y >= std::min(p.Y, r.Y))
            return true;

        return false;
    }

    // To find orientation of ordered triplet (p, q, r).
    // The function returns following values
    // 0 --> p, q and r are collinear
    // 1 --> Clockwise
    // 2 --> Counterclockwise
    int orientation(Point p, Point q, Point r)
    {
        // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
        // for details of below formula.
        int val = (q.Y - p.Y) * (r.X - q.X) -
                (q.X - p.X) * (r.Y - q.Y);

        if (val == 0) return 0; // collinear

        return (val > 0)? 1: 2; // clock or counterclock wise
    }


    bool intersects(Segment otherSegment){
        // Find the four orientations needed for general and
        // special cases
        int o1 = orientation(p0, p1, otherSegment.p0);
        int o2 = orientation(p0, p1, otherSegment.p1);
        int o3 = orientation(otherSegment.p0, otherSegment.p1, p0);
        int o4 = orientation(otherSegment.p0, otherSegment.p1, p1);

        // General case
        if (o1 != o2 && o3 != o4)
            return true;

        // Special Cases
        // p0, p1 and otherSegment.p0 are collinear and otherSegment.p0 lies on segment p0p1
        if (o1 == 0 && onSegment(p0, otherSegment.p0, p1)) return true;

        // p0, p1 and otherSegment.p1 are collinear and otherSegment.p1 lies on segment p0p1
        if (o2 == 0 && onSegment(p0, otherSegment.p1, p1)) return true;

        // otherSegment.p0, otherSegment.p1 and p0 are collinear and p0 lies on segment otherSegment.p0otherSegment.p1
        if (o3 == 0 && onSegment(otherSegment.p0, p0, otherSegment.p1)) return true;

        // otherSegment.p0, otherSegment.p1 and p1 are collinear and p1 lies on segment otherSegment.p0otherSegment.p1
        if (o4 == 0 && onSegment(otherSegment.p0, p1, otherSegment.p1)) return true;

        return false; // Doesn't fall in any of the above cases
    }
};


namespace boost {
	namespace polygon {
		template <>
		struct geometry_concept<Point> { typedef point_concept type; };

		template <>
		struct point_traits<Point> {
			typedef int coordinate_type;

			static inline coordinate_type get(const Point& point, orientation_2d orient) {
				return (orient == HORIZONTAL) ? point.X : point.Y;
			}
		};

		template <>
		struct geometry_concept<Segment> { typedef segment_concept type; };

		template <>
		struct segment_traits<Segment> {
			typedef Segment segment_type;
			typedef Point point_type;
			typedef int coordinate_type;

			static point_type get(const segment_type& segment, direction_1d dir) {
				return dir.to_int() ? segment.p1 : segment.p0;
			}
		};
	}
}

struct c_Vertex {
	double X;
	double Y;

	c_Vertex(double x = 0, double y = 0) : X(x), Y(y) {}

	bool operator==(c_Vertex& other){
	    return (X == other.X) && (Y==other.Y);
	}
};


struct c_Edge {
	long long start;
	long long end;

	bool isPrimary;
	bool isLinear;

    long long cell;
    long long twin;

	c_Edge(long long start = -1, long long end = -1, bool isPrimary = false, bool isLinear = false, long long cell = -1, long long twin = -1) {
		this->start = start;
		this->end = end;
		this->isPrimary = isPrimary;
		this->isLinear = isLinear;
		this->cell = cell;
		this->twin = twin;
	}
};

//A structure to identify a segment as part of one cell only
struct c_Cell{
    long long cell_identifier;
	long long site;
	bool contains_point;
	bool contains_segment;
	bool is_open;
	bool is_degenerate;
	std::vector<long long> vertices;
	std::vector<long long> edges;

	int source_category;

	c_Cell(size_t cell_identifier = -1, size_t site = -1, bool contains_point = false, bool contains_segment = false, bool is_open = false, int source_category = -1){
        this->cell_identifier = cell_identifier;
		this->site = site;
		this->contains_point = contains_point;
		this->contains_segment = contains_segment;
		this->is_open = is_open;
		this->source_category = source_category;
	}
};

using namespace boost::polygon;

class VoronoiDiagram {
public:
	VoronoiDiagram();
	void AddPoint(Point p);
	void AddSegment(Segment s);
	void Construct();

	std::vector<Point> GetPoints();
	std::vector<Segment> GetSegments();

	long long CountVertices();
	long long CountEdges();
	long long CountCells();

    std::vector<int> GetIntersectingSegments();

	//Map index to vertex
	typedef std::pair<long long, const voronoi_diagram<double>::vertex_type*> index_to_vertex;
	std::map<long long, const voronoi_diagram<double>::vertex_type*> map_indexes_to_vertices;

	//Map vertex to index
	std::map <const voronoi_diagram<double>::vertex_type*, long long> map_vertices_to_indexes;
	typedef std::pair<const voronoi_diagram<double>::vertex_type*, long long> vertex_to_index;

	//Map index to edge
	typedef std::pair<long long, const voronoi_diagram<double>::edge_type*> index_to_edge;
	std::map<long long, const voronoi_diagram<double>::edge_type*> map_indexes_to_edges;

	//Map edge to index
	std::map <const voronoi_diagram<double>::edge_type*, long long> map_edges_to_indexes;
	typedef std::pair<const voronoi_diagram<double>::edge_type*, long long> edge_to_index;

	//Map index to cell
	typedef std::pair<long long, const voronoi_diagram<double>::cell_type*> index_to_cell;
	std::map<long long, const voronoi_diagram<double>::cell_type*> map_indexes_to_cells;

	//Map cell to index
	std::map <const voronoi_diagram<double>::cell_type*, long long> map_cells_to_indexes;
	typedef std::pair<const voronoi_diagram<double>::cell_type*, long long> cell_to_index;

	void MapVertexIndexes();
	void MapEdgeIndexes();
	void MapCellIndexes();

	c_Vertex GetVertex(long long index);
	c_Edge GetEdge(long long index);
	c_Cell GetCell(long long index);


private:
	std::vector<Point> points;
	std::vector<Segment> segments;
	voronoi_diagram<double> vd;

};
