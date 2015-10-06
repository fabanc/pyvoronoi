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
	size_t start;
	size_t end;

	bool isPrimary;

	size_t site1;
	size_t site2;

	c_Edge(size_t start = -1, size_t end = -1, bool isPrimary = false, size_t site1 = -1, size_t site2 = -1) {
		this->start = start;
		this->end = end;
		this->isPrimary = isPrimary;
		this->site1 = site1;
		this->site2 = site2;
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
	std::vector<Point> GetPoints();
	std::vector<Segment> GetSegments();
private:
	std::vector<Point> points;
	std::vector<Segment> segments;
	voronoi_diagram<double> vd;
};