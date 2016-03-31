// Voronoi.cpp : Defines the entry point for the console application.
//

#define BOOST_POLYGON_NO_DEPS
#define BOOST_NO_USER_CONFIG
#define BOOST_NO_COMPILER_CONFIG
#define BOOST_NO_STDLIB_CONFIG
#define BOOST_NO_PLATFORM_CONFIG
#define BOOST_HAS_STDINT_H

#define __GLIBC__ 0

#include "boost/polygon/voronoi.hpp"

struct Point {
	int X;
	int Y;
	Point(int x = 0, int y = 0) : X(x), Y(y) {}
};

struct Segment {
	Point p0;
	Point p1;
	Segment(Point a = Point(), Point b = Point()) : p0(a.X, a.Y), p1(b.X, b.Y) {}
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
};


struct c_Edge {
	long long start;
	long long end;

	bool isPrimary;

	size_t site1;
	size_t site2;
	
	bool isLinear;
	
    long cell;
    long twin;	
	
	c_Edge(long long start = -1, long long end = -1, bool isPrimary = false, size_t site1 = -1, size_t site2 = -1, bool isLinear = false, long cell = -1, long twin = -1) {
		this->start = start;
		this->end = end;
		this->isPrimary = isPrimary;
		this->site1 = site1;
		this->site2 = site2;
		this->isLinear = isLinear;
		this->cell = cell;
		this->twin = twin;	
	}
};

//A structure to identify a segment as part of one cell only
struct c_Cell{
    size_t cell_identifier;
	size_t site;
	bool contains_point;
	bool contains_segment;
	bool is_open;
	
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
	void GetEdges(std::vector<c_Vertex> &, std::vector<c_Edge> &);
  void GetCells(std::vector<c_Vertex> &, std::vector<c_Edge> &, std::vector<c_Cell> &);
	std::vector<Point> GetPoints();
	std::vector<Segment> GetSegments();
private:
	std::vector<Point> points;
	std::vector<Segment> segments;
	voronoi_diagram<double> vd;
};