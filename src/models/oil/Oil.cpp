#include "src/models/oil/Oil.hpp"
#include "src/util/utils.h"

#define ADOLC_ADVANCED_BRANCHING

using namespace oil;
using namespace std;
using namespace point;

Oil::Oil()
{
}
Oil::~Oil()
{
	delete[] x, h;
}
void Oil::setProps(const Properties& props)
{
	R_dim = props.R_dim;
	r_w = props.r_w;
	r_e = props.r_e;
	perfIntervals = props.perfIntervals;

	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	skeletonsNum = props.props_sk.size();
	props_sk = props.props_sk;
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].kx = MilliDarcyToM2(props_sk[j].kx);
		props_sk[j].ky = MilliDarcyToM2(props_sk[j].ky);
	}

	periodsNum = props.timePeriods.size();
	for (int i = 0; i < periodsNum; i++)
	{
		period.push_back(props.timePeriods[i]);

		if (leftBoundIsRate)
			rate.push_back(props.rates[i] / 86400.0);
		else
			pwf.push_back(props.pwf[i]);
	}

	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);

	alpha = props.alpha;

	makeDimLess();
}
void Oil::makeDimLess()
{
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;

	for (auto& props : props_sk)
	{
		props.kx /= (R_dim * R_dim);
		props.ky /= (R_dim * R_dim);
		props.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		props.beta /= (1.0 / P_dim);
		props.height /= R_dim;
		props.p_init /= P_dim;
		props.p_out /= P_dim;
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		if (leftBoundIsRate)
			rate[i] /= Q_dim;
		else
			pwf[i] /= P_dim;
	}

	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_ref /= P_dim;

	alpha /= t_dim;
}
void Oil::setInitialState()
{
	const auto& props = props_sk[0];
	for (size_t i = 0; i < cellsNum; i++)
	{
		const auto& cell = mesh->getCell(i);
		auto data = (*this)[i];
		data.u_prev.p = data.u_iter.p = data.u_next.p = props.p_init;
	}

	x = new TapeVariable[cellsNum];
	h = new adouble[var_size * cellsNum];

	mesh->get_XY_perm = bind(&Oil::getPerm_XY, this, placeholders::_1);
	setTrans();
}
void Oil::setPeriod(const int period)
{
	/*if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum * cells[it->first].hz / height_perf;
		}
		else {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = it->second * Q_sum / rate[period - 1];
		}
	}
	else*/
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}
}
double Oil::getRate(const size_t cell_idx)
{
	return 0.0;
}

void Oil::setTrans()
{
	mesh->calc_transmissibilities();
}

adouble Oil::solveInner(const Cell& cell)
{
	assert(cell.type == elem::HEX || cell.type == elem::PRISM);
	const auto& cur = x[cell.id];
	const auto prev = (*this)[cell.id].u_prev;

	adouble tmp;
	adouble H = props_sk[0].getPoro(cur.p) * props_oil.getDensity(cur.p) - props_sk[0].getPoro(prev.p) * props_oil.getDensity(prev.p);

	// Bottom
	const auto& nebr_bot = cell.nebrs[0];
	const auto& beta_bot = mesh->cells[nebr_bot.nebr.cell];
	adouble isInner_bot = (beta_bot.type != elem::BORDER_QUAD && beta_bot.type != elem::BORDER_TRI) ? true : false;
	const auto& anebr_bot = x[beta_bot.id];
	tmp = 0.0;
	condassign(tmp, isInner_bot, ht / cell.V * getVertTrans(cell, 0, beta_bot, nebr_bot.nebr.nebr) *
					linearInterp1d(props_oil.getDensity(cur.p) / props_oil.getViscosity(cur.p), nebr_bot.L,
						props_oil.getDensity(anebr_bot.p) / props_oil.getViscosity(anebr_bot.p), beta_bot.nebrs[nebr_bot.nebr.nebr].L) * (cur.p - anebr_bot.p));
	H += tmp;
	// Top
	const auto& nebr_top = cell.nebrs[1];
	const auto& beta_top = mesh->cells[nebr_top.nebr.cell];
	adouble isInner_top = (beta_top.type != elem::BORDER_QUAD && beta_top.type != elem::BORDER_TRI) ? true : false;
	const auto& anebr_top = x[beta_top.id];
	tmp = 0.0;
	condassign(tmp, isInner_top, ht / cell.V * getVertTrans(cell, 1, beta_top, nebr_top.nebr.nebr) *
					linearInterp1d(props_oil.getDensity(cur.p) / props_oil.getViscosity(cur.p), nebr_top.L,
						props_oil.getDensity(anebr_top.p) / props_oil.getViscosity(anebr_top.p), beta_top.nebrs[nebr_top.nebr.nebr].L) * (cur.p - anebr_top.p));
	H += tmp;
	// Laterals
	for (int i = 2; i < cell.nebrs_num; i++)
	{
		const auto& nebr = cell.nebrs[i];
		const Cell& beta = mesh->cells[nebr.nebr.cell];
		tmp = ht / cell.V * linearInterp(nebr, props_oil.getDensity(x[nebr.nearest_cells[0]].p) / props_oil.getViscosity(x[nebr.nearest_cells[0]].p),
			props_oil.getDensity(x[nebr.nearest_cells[1]].p) / props_oil.getViscosity(x[nebr.nearest_cells[1]].p),
			props_oil.getDensity(x[nebr.nearest_cells[2]].p) / props_oil.getViscosity(x[nebr.nearest_cells[2]].p));
		const auto& ireg1 = nebr.ireg[PLUS];
		for (int j = 0; j < ireg1->trans.size(); j++)
			H += tmp * ireg1->trans[nebr.ireg_id[PLUS]][j] * x[ireg1->cells[j].cell].p;
		const auto& ireg2 = nebr.ireg[MINUS];
		for (int j = 0; j < ireg2->trans.size(); j++)
			H += tmp * ireg2->trans[nebr.ireg_id[MINUS]][j] * x[ireg2->cells[j].cell].p;
	}
	return H;
}
adouble Oil::solveBorder(const Cell& cell)
{
	assert(cell.type == elem::BORDER_HEX);
	const auto& cur = x[cell.id];
	const auto& nebr = x[cell.nebrs[0].nebr.cell];
	//adouble H;
	//adouble rightIsPres = rightBoundIsPres;
	//condassign(H, rightIsPres, (cur.p - (adouble)(props_sk[0].p_out)) / P_dim, (cur.p - nebr.p) / P_dim);
	return cur.p - (adouble)(props_sk[0].p_out);
}
adouble Oil::solveFrac(const Cell& cell)
{
	assert(cell.type == elem::FRAC_HEX);
	return x[cell.id].p - Pwf;
}