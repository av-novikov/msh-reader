#ifndef MESH_HPP_
#define MESH_HPP_

#include "src/util/Point.hpp"
#include "src/mesh/Elem.hpp"
#include <vector>
#include <cassert>
#include <initializer_list>

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
	static const int stencil = -11;

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
		void set_interaction_regions();

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
