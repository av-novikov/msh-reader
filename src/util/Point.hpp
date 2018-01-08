#ifndef POINT_HPP_
#define POINT_HPP_

#include <cmath>
#include <iostream>
#include <vector>
#include <new>

#define EQUALITY_TOLERANCE 1.E-9

namespace point
{
	enum HalfType { PLUS, MINUS };
	struct Id
	{
		int cell;	char nebr;
	};
	struct RegCellId
	{
		int cell;
		char nebr [2];
	};
	inline bool operator==(const RegCellId& id, const int _cell)
	{
		return (id.cell == _cell);
	};
	struct Interaction
	{
		std::vector<RegCellId> cells;
		std::vector<int> row_offset, colms;
		std::vector<double> trans;

		Interaction() {};
		Interaction(const std::vector<int>& _cells)
		{
			for (const auto& cell : _cells)
				cells.push_back(RegCellId{ cell, {-1, -1} });
		};
		~Interaction() {};
		Interaction& operator=(const Interaction& other)
		{
			cells = other.cells;
			row_offset = other.row_offset;
			colms = other.colms;
			trans = other.trans;
			return *this;
		};
	};
	enum PointType {INNER, BORDER};
	struct Point
	{
		PointType type;
		const int id;
		union
		{
			double coords[3];
			struct
			{
				double x;	double y;	double z;
			};
		};
		std::vector<int> cells;
		Interaction* int_reg;

		Point() : id(-1), int_reg(nullptr), type(INNER) {};
		Point(const double _x, const double _y, const double _z) : id(-1), x(_x), y(_y), z(_z), int_reg(nullptr), type(INNER) { };
		Point(const int _id, const double _x, const double _y, const double _z) : id(_id), x(_x), y(_y), z(_z), int_reg(nullptr), type(INNER) { };
		Point(const int _id, const double _x, const double _y, const double _z, const std::vector<int>& _cells) : Point(_id, _x, _y, _z)
		{
			this->cells = _cells;
		};
		~Point() 
		{ 
			if(int_reg != nullptr)
				delete int_reg; 
		};
		Point(const Point& a) : id(a.id)
		{
			(*this) = a;
		};
		Point& operator=(const Point& rhs)
		{
			x = rhs.x, y = rhs.y, z = rhs.z;
			cells = rhs.cells;
			int_reg = rhs.int_reg;
			type = rhs.type;
			return *this;
		};
		Point& operator/=(const double k)
		{
			x /= k;	y /= k;	z /= k;
			return *this;
		};
		Point& operator+=(const Point& rhs)
		{
			x += rhs.x;	y += rhs.y;	z += rhs.z;
			return *this;
		};
		Point& operator-=(const Point& rhs)
		{
			x -= rhs.x;	y -= rhs.y;	z -= rhs.z;
		};

		inline double norm() const { return sqrt(x * x + y * y + z * z); };
	};
	inline std::ostream& operator<<(std::ostream& os, const Point& a)
	{
		os << a.x << " " << a.y << " " << a.z << std::endl;
		return os;
	}
	inline bool operator==(const Point& a1, const Point& a2)
	{
		if ((fabs(a2.x - a1.x) > EQUALITY_TOLERANCE) ||
			(fabs(a2.y - a1.y) > EQUALITY_TOLERANCE) ||
			(fabs(a2.z - a1.z) > EQUALITY_TOLERANCE))
			return false;
		else
			return true;
	};
	inline Point operator-(const Point& rhs)
	{
		return Point(-rhs.x, -rhs.y, -rhs.z);
	};
	inline Point operator-(const Point& a1, const Point& a2)
	{
		return Point(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
	};
	inline Point operator+(const Point& rhs)
	{
		return Point(rhs.x, rhs.y, rhs.z);
	};
	inline Point operator+(const Point& a1, const Point& a2)
	{
		return Point(a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
	};
	inline Point operator*(const Point& a1, double k)
	{
		return Point(a1.x * k, a1.y * k, a1.z * k);
	};
	inline Point operator*(double k, const Point& a1)
	{
		return a1 * k;
	};
	inline Point operator/(const Point& a1, double k)
	{
		return Point(a1.x / k, a1.y / k, a1.z / k);
	};
	inline Point operator/(const Point& a1, const Point& a2)
	{
		return Point(a1.x / a2.x, a1.y / a2.y, a1.z / a2.z);
	};
	inline Point operator*(const Point& a1, const Point& a2)
	{
		return Point(a1.x * a2.x, a1.y * a2.y, a1.z * a2.z);
	};

	inline double dot_product(const Point& a1, const Point& a2)
	{
		return a1.x * a2.x + a1.y * a2.y + a1.z * a2.z;
	};
	inline Point vector_product(const Point& a1, const Point& a2)
	{
		return{ a1.y * a2.z - a1.z * a2.y,
				a1.z * a2.x - a1.x * a2.z,
				a1.x * a2.y - a1.y * a2.x };
	};
	inline double distance(const Point& a1, const Point& a2)
	{
		return sqrt(dot_product(a2 - a1, a2 - a1));
	};
	inline double square(const Point& a1, const Point& a2, const Point& a3)
	{
		const double a = distance(a1, a2);
		const double b = distance(a2, a3);
		const double c = distance(a3, a1);
		const double p = (a + b + c) / 2.0;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	};
	inline double square(const Point& a1, const Point& a2, const Point& a3, const Point& a4)
	{
		return square(a1, a2, a3) + square(a3, a4, a1);
	};
};

#endif /* POINT_HPP_ */