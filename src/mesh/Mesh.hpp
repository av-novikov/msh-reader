#ifndef MESH_HPP_
#define MESH_HPP_

#include "src/util/Point.hpp"
#include <vector>
#include <cassert>
#include <initializer_list>
#include <array>

namespace elem
{
	static const int MAX_ELEM_POINT_SIZE = 8;
	static const int MAX_ELEM_NEBR_SIZE = 6;
	static const int QUAD_VERT_SIZE = 4;
	static const int TRI_VERT_SIZE = 3;
	static const int MAX_INTERFACES_STENCIL = 6;

	enum EType { BORDER_TRI, BORDER_QUAD, FRAC_QUAD, PRISM, HEX, BORDER_HEX};
	inline const char num_of_verts(const EType type)
	{
		if (type == BORDER_TRI)
			return 3;
		else if (type == BORDER_QUAD || type == FRAC_QUAD)
			return 4;
		else if (type == PRISM)
			return 6;
		else if (type == HEX || type == BORDER_HEX)
			return 8;
	}
	inline const char num_of_nebrs(const EType type)
	{
		if (type == BORDER_TRI || type == BORDER_QUAD)
			return 1;
		else if (type == FRAC_QUAD)
			return 2;
		else if (type == PRISM)
			return 5;
		else if (type == HEX || type == BORDER_HEX)
			return 6;
	}
	struct Id { int cell; char nebr; };
	struct Nebr
	{
		Id nebr;
		std::array<Id, MAX_INTERFACES_STENCIL> stencilL, stencilR;
		double S;
		double L;
		point::Point cent;
		double omega1, omega2;
		point::Point n, nu;
	};

	class Element
	{
	public:
		int id;
		EType type;

		char verts_num;
		char nebrs_num;
		std::array<int, MAX_ELEM_POINT_SIZE> verts;
		std::array<Nebr, MAX_ELEM_NEBR_SIZE> nebrs;
		point::Point cent;
		double V;

		Element() {};
		Element(const EType _type, const int* _verts) : type(_type), verts_num(num_of_verts(_type)), nebrs_num(num_of_nebrs(_type))
		{
			for (int i = 0; i < verts_num; i++)
				verts[i] = _verts[i];
		};
		Element(const EType _type, const int* _verts, const int* _nebrs) : type(_type), verts_num(num_of_verts(_type)), nebrs_num(num_of_nebrs(_type))
		{
			for (int i = 0; i < verts_num; i++)
				verts[i] = _verts[i];
			for (int i = 0; i < nebrs_num; i++)
				nebrs[i].nebr.cell = _nebrs[i];
		};
		~Element() {};
	};
};
namespace mshreader
{
	class MshReader;
};
namespace snapshotter
{
	template<class modelType> class VTKSnapshotter;
};

namespace grid
{
	static const int stencil = 11;
	enum HalfType { LEFT, RIGHT };

	class Mesh
	{
		friend class mshreader::MshReader;
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractSolver;
	public:
		typedef elem::Element Cell;

		int inner_size, border_size, frac_size, pts_size;
		std::vector<point::Point> pts;
		std::vector<Cell> cells;
	protected:
		void set_geom_props();
		void count_types();
		void setNebrId();
		inline const char findNebrId(const Cell& cell_cur, const Cell& cell_nebr)
		{
			for (char i = 0; i < cell_nebr.nebrs_num; i++)
				if (cell_cur.id == cell_nebr.nebrs[i].nebr.cell)
					return i;
			exit(-1);
		};

	public:
		double Volume;

		Mesh();
		~Mesh();

		void process_geometry();

		size_t getCellsSize() const;
		inline Cell& getCell(const int i) {	return cells[i]; };
		inline const Cell& getCell(const int i) const { return cells[i]; };
	};
}

#endif /* MESH_HPP_ */
