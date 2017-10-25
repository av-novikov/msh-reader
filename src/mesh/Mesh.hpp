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

	enum EType { BORDER_TRI, BORDER_QUAD, FRAC_QUAD, PRISM, HEX };

	inline const int num_of_verts(const EType type)
	{
		if (type == EType::BORDER_TRI)
			return 3;
		else if (type == EType::BORDER_QUAD || type == EType::FRAC_QUAD)
			return 4;
		else if (type == EType::PRISM)
			return 6;
		else if (type == EType::HEX)
			return 8;
	}
	inline const int num_of_nebrs(const EType type)
	{
		if (type == EType::BORDER_TRI || type == EType::BORDER_QUAD)
			return 1;
		else if (type == EType::FRAC_QUAD)
			return 2;
		else if (type == EType::PRISM)
			return 5;
		else if (type == EType::HEX)
			return 6;
	}
	struct Nebr
	{
		int id;	
		double S;
		double L;
		point::Point cent;
	};

	class Element
	{
	public:
		int id;
		EType type;

		int verts_num;
		int nebrs_num;
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
				nebrs[i].id = _nebrs[i];
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
	class Mesh
	{
		friend class mshreader::MshReader;
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractSolver;
	public:
		typedef elem::Element Cell;
	protected:
		int inner_size, border_size, frac_size, pts_size;
		std::vector<point::Point> pts;
		std::vector<Cell> cells;

		int check_neighbors() const;
		bool are_adjanced(const Cell& el1, const Cell& el2);
		void set_geom_props();
		void count_types();

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
