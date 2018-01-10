#ifndef OIL_HPP_
#define OIL_HPP_

#include "src/models/Variables.hpp"
#include "src/models/AbstractModel.hpp"
#include "src/models/Oil/Properties.hpp"

namespace oil
{
	typedef var::containers::TapeVar1Phase TapeVariable;
	class Oil : public AbstractModel<var::containers::Var1phase, Properties, var::BasicVariables, Oil>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		friend class OilSolver;
	protected:
		void setTrans();
		void setProps(const Properties& props);
		void makeDimLess();
		void setInitialState();

		TapeVariable* x;
		adouble* h;

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		grid::Mesh::Perm_XY getPerm_XY(const Cell& cell) 
		{
			return{ props_sk[0].kx, props_sk[0].ky };
		};
		inline const double getKz(const Cell& cell) const
		{
			return props_sk[0].kz;
		};

		inline adouble linearInterp(const Cell& cell, const int nebr_id, adouble f0, adouble f1, adouble f2)
		{
			const auto& nebr = cell.nebrs[nebr_id];
			return f0 * nebr.nearest_coeff[0] + f1 * nebr.nearest_coeff[1] + f2 * nebr.nearest_coeff[2];
		}

		/*double getTrans(const Cell& cell1, const int idx, const Cell& cell2) const
		{
			const double k1 = getPerm(cell1);
			const double k2 = getPerm(cell2);
			const double dist1 = getDistance(cell1, cell2);
			const double dist2 = getDistance(cell2, cell1);
			if(cell1.type == CellType::WELL)
				return props_sk[0].height * mesh->wellNebrs[idx].length * k1 * k2 / (k1 * dist2 + k2 * dist1);
			else
				return props_sk[0].height * cell1.length[idx] * k1 * k2 / (k1 * dist2 + k2 * dist1);
		};
		const double& getDistance(const Cell& cell, const Cell& beta) const
		{
			if (cell.type == CellType::WELL && beta.type != CellType::WELL)
			{
				for (const auto& nebr : mesh->wellNebrs)
					if (nebr.id == beta.id)
						return nebr.dist;
				return 0.0;
			}
			else
				return cell.getDistance(beta.id);
		};*/

		adouble solveInner(const Cell& cell);
		adouble solveBorder(const Cell& cell);
		adouble solveFrac(const Cell& cell);
	public:
		Oil();
		~Oil();

		void setPeriod(const int period);
		static const int var_size = VarContainer::size;

		double getRate(const size_t cell_idx);
	};
};

#endif /* OIL_HPP_ */
