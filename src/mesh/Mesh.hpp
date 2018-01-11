#ifndef MESH_HPP_
#define MESH_HPP_

#include "src/util/Point.hpp"
#include "src/mesh/Elem.hpp"
#include <vector>
#include <cassert>
#include <initializer_list>
#include <functional>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

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
	static const int stencil = 12;
	const int MAX_REG_CELLS = 30;
	const double REL_TOL = 1.E-10;

	class Mesh
	{
		friend class mshreader::MshReader;
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractSolver;
	public:
		typedef elem::Element Cell;
		struct Perm_XY { double kx, ky; };
		typedef std::function<Perm_XY(const Cell&)> Perm_Getter;

		int inner_size, border_size, frac_size, pts_size;
		std::vector<point::Point> pts;
		std::vector<Cell> cells;
		Perm_Getter get_XY_perm;
	protected:
		double z_max;

		void setFrac();
		void set_geom_props();
		void count_types();
		void set_nearest();
		void setNebrId();
		inline const char findNebrId(const Cell& cell_cur, const Cell& cell_nebr)
		{
			for (char i = 0; i < cell_nebr.nebrs_num; i++)
				if (cell_cur.id == cell_nebr.nebrs[i].nebr.cell)
					return i;
			exit(-1);
		};
		void set_interaction_regions();

		adouble **buf_a, **buf_b, **buf_c, **buf_d, **inv_a, **t, **mult;
	public:
		double Volume;

		Mesh();
		~Mesh();

		void process_geometry();

		size_t getCalcCellsSize() const;
		inline Cell& getCell(const int i) {	return cells[i]; };
		inline const Cell& getCell(const int i) const { return cells[i]; };
		void calc_transmissibilities();
	};
}

#endif /* MESH_HPP_ */
