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
#include <cmath>

struct IntersectionPoint{
	double X;
	double Y;
	IntersectionPoint(double x = 0, double y = 0) : X(x), Y(y) {}
};

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
    bool onSegment(Point p, IntersectionPoint q, Point r)
    {
        // Otherwise, it is in the pr space, returns true
        if (q.X <= std::max(p.X, r.X) && q.X >= std::min(p.X, r.X) && q.Y <= std::max(p.Y, r.Y) && q.Y >= std::min(p.Y, r.Y))
            return true;

        return false;
    }

    bool onSegment(Point p, Point q, Point r)
    {
        if (q.X <= std::max(p.X, r.X) && q.X >= std::min(p.X, r.X) && q.Y <= std::max(p.Y, r.Y) && q.Y >= std::min(p.Y, r.Y))
            return true;
        return false;
    }

    bool onSegment2(Point p, Point q, Point r)
    {
        if (q.X < std::max(p.X, r.X) && q.X > std::min(p.X, r.X) && q.Y < std::max(p.Y, r.Y) && q.Y > std::min(p.Y, r.Y))
            return true;
        return false;
    }

    bool onEndpoint(Point q){
        if (q.X == p0.X && q.Y == p0.Y)
            return true;
        if (q.X == p1.X && q.Y == p1.Y)
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

    //  Returns Point of intersection if do intersect otherwise default Point (null)
    bool findIntersection(Segment otherSegment)
    {
        double tolerance = 0.01;
        long long x1 = p0.X, y1 = p0.Y;
        long long x2 = p1.X, y2 = p1.Y;

        long long x3 = otherSegment.p0.X, y3 = otherSegment.p0.Y;
        long long x4 = otherSegment.p1.X, y4 = otherSegment.p1.Y;

        // equations of the form x=c (two vertical lines) with overlapping
        if (abs(x1 - x2) < tolerance && abs(x3 - x4) < tolerance && abs(x1 - x3) < tolerance)
        {
            //throw new Exception("Both lines overlap vertically, ambiguous intersection points.");
            if (onEndpoint(otherSegment.p0) || onEndpoint(otherSegment.p1))
                return false;
            else{
                return true;
            }
        }

        //equations of the form y=c (two horizontal lines) with overlapping
        if (abs(y1 - y2) < tolerance && abs(y3 - y4) < tolerance && abs(y1 - y3) < tolerance)
        {
            //throw new Exception("Both lines overlap horizontally, ambiguous intersection points.");
            if (onEndpoint(otherSegment.p0) || onEndpoint(otherSegment.p1))
                return false;
            return true;
        }

        //equations of the form x=c (two vertical parallel lines)
        if (abs(x1 - x2) < tolerance && abs(x3 - x4) < tolerance)
        {
            //return default (no intersection)
            if (onEndpoint(otherSegment.p0) || onEndpoint(otherSegment.p1))
                return false;
            return true;
        }

        //equations of the form y=c (two horizontal parallel lines)
        if (abs(y1 - y2) < tolerance && abs(y3 - y4) < tolerance)
        {
            //return default (no intersection)
            if (onEndpoint(otherSegment.p0) || onEndpoint(otherSegment.p1))
                return false;
            return true;
        }


        if (orientation(p0, otherSegment.p0, p1) == 0 && orientation(p0, otherSegment.p1, p1) == 0){
            if (
                onSegment2(p0, otherSegment.p0, p1) ||
                onSegment2(p0, otherSegment.p1, p1) ||
                onSegment2(otherSegment.p0, p0, otherSegment.p1) ||
                onSegment2(otherSegment.p0, p1, otherSegment.p1)
                ){

                return true;
               }
        }


        //general equation of line is y = mx + c where m is the slope
        //assume equation of line 1 as y1 = m1x1 + c1
        //=> -m1x1 + y1 = c1 ----(1)
        //assume equation of line 2 as y2 = m2x2 + c2
        //=> -m2x2 + y2 = c2 -----(2)
        //if line 1 and 2 intersect then x1=x2=x & y1=y2=y where (x,y) is the intersection point
        //so we will get below two equations
        //-m1x + y = c1 --------(3)
        //-m2x + y = c2 --------(4)

        double x, y;

        //lineA is vertical x1 = x2
        //slope will be infinity
        //so lets derive another solution
        if (abs(x1 - x2) < tolerance)
        {
            //compute slope of line 2 (m2) and c2
            double m2 = (double)(y4 - y3) / (double)(x4 - x3);
            double c2 = -m2 * x3 + y3;

            //equation of vertical line is x = c
            //if line 1 and 2 intersect then x1=c1=x
            //subsitute x=x1 in (4) => -m2x1 + y = c2
            // => y = c2 + m2x1
            x = (double)x1;
            y = c2 + m2 * x1;
        }
        //otherSegment is vertical x3 = x4
        //slope will be infinity
        //so lets derive another solution
        else if (abs(x3 - x4) < tolerance)
        {
            //compute slope of line 1 (m1) and c2
            double m1 = (double)(y2 - y1) / (double)(x2 - x1);
            double c1 = -m1 * x1 + y1;

            //equation of vertical line is x = c
            //if line 1 and 2 intersect then x3=c3=x
            //subsitute x=x3 in (3) => -m1x3 + y = c1
            // => y = c1 + m1x3
            x = (double)x3;
            y = c1 + m1 * x3;
        }
        //lineA & otherSegment are not vertical
        //(could be horizontal we can handle it with slope = 0)
        else
        {
            //compute slope of line 1 (m1) and c2
            double m1 = (double)(y2 - y1) / (double)(x2 - x1);
            double c1 = -m1 * x1 + y1;

            //compute slope of line 2 (m2) and c2
            double m2 = (double)(y4 - y3) / (double)(x4 - x3);
            double c2 = -m2 * x3 + y3;

            //solving equations (3) & (4) => x = (c1-c2)/(m2-m1)
            //plugging x value in equation (4) => y = c2 + m2 * x
            x = (c1 - c2) / (m2 - m1);
            y = c2 + m2 * x;

            //verify by plugging intersection point (x, y)
            //in orginal equations (1) & (2) to see if they intersect
            //otherwise x,y values will not be finite and will fail this check
            if (!(abs(-m1 * x + y - c1) < tolerance
                && abs(-m2 * x + y - c2) < tolerance))
            {
                //return default (no intersection)
                return false;
            }
        }

        //x,y can intersect outside the line segment since line is infinitely long
        //so finally check if x, y is within both the line segments
        IntersectionPoint p =  IntersectionPoint(x, y);
        if (onSegment(p0, p, p1) && onSegment(otherSegment.p0, p, otherSegment.p1))
        {
            Point pi = Point((int)round(x), (int)round(y));
            if (onEndpoint(pi) && otherSegment.onEndpoint(pi))
                return false;
            return true;

        }

        //return default (no intersection)
        return false;
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

    long long CountPoints();
    long long CountSegments();
	long long CountVertices();
	long long CountEdges();
	long long CountCells();

    std::vector<unsigned long long> GetIntersectingSegments();
    std::vector<unsigned long long> GetDegenerateSegments();
    std::vector<unsigned long long> GetPointsOnSegments();

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

    Point GetPoint(int index);
    Segment GetSegment(int index);
	c_Vertex GetVertex(long long index);
	c_Edge GetEdge(long long index);
	c_Cell GetCell(long long index);


private:
	std::vector<Point> points;
	std::vector<Segment> segments;
	voronoi_diagram<double> vd;

};
