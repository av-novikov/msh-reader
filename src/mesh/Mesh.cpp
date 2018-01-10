#include "src/mesh/Mesh.hpp"
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <utility>
#define _USE_MATH_DEFINES
#include <math.h>

#include "src/util/utils.h"

using namespace grid;
using namespace point;
using namespace std;

Mesh::Mesh()
{
	buf_a = new adouble*[MAX_REG_CELLS];
	buf_b = new adouble*[MAX_REG_CELLS];
	buf_c = new adouble*[MAX_REG_CELLS];
	buf_d = new adouble*[MAX_REG_CELLS];
	t =		new adouble*[MAX_REG_CELLS];
	mult =  new adouble*[MAX_REG_CELLS];
	for (int i = 0; i < MAX_REG_CELLS; i++)
	{
		buf_a[i] = new adouble[MAX_REG_CELLS];
		buf_b[i] = new adouble[MAX_REG_CELLS];
		buf_c[i] = new adouble[MAX_REG_CELLS];
		buf_d[i] = new adouble[MAX_REG_CELLS];
		t[i] = new adouble[MAX_REG_CELLS];
		mult[i] = new adouble[MAX_REG_CELLS];
	}
}
Mesh::~Mesh()
{
	for (int i = 0; i < MAX_REG_CELLS; i++)
		delete[] buf_a[i], buf_b[i], buf_c[i], buf_d[i], t[i], mult[i];
	delete[] buf_a, buf_b, buf_c, buf_d, t, mult;
}
void Mesh::process_geometry()
{
	setNebrId();
	set_geom_props();
	count_types();

	set_interaction_regions(); 
	set_nearest();
}
void Mesh::count_types()
{
	inner_size = border_size = frac_size = 0;
	Volume = 0.0;
	for (const auto& el : cells)
	{
		Volume += el.V;
		if (el.type == elem::HEX || 
			el.type == elem::BORDER_HEX || 
			el.type == elem::PRISM)							inner_size++;
		else if (el.type == elem::FRAC_QUAD)				frac_size++;
		else												border_size++;
	}
}
void Mesh::set_nearest()
{
	double dist = 0.0;

	// Find nearest cell centers for all inner faces
	for (int i = 0; i < inner_size; i++)
	{
		auto& cell = cells[i];
		
		if (cell.type == elem::HEX || cell.type == elem::PRISM)
			for (int j = 2; j < cell.nebrs_num; j++)
			{
				auto& nebr = cell.nebrs[j];
				nebr.nearest_cells[0] = cell.id;
				nebr.nearest_dists[0] = point::distance(nebr.cent, cells[nebr.nearest_cells[0]].cent);
				nebr.nearest_cells[1] = nebr.nebr.cell;
				nebr.nearest_dists[1] = point::distance(nebr.cent, cells[nebr.nearest_cells[1]].cent);

				const auto& cell1 = cells[min_element(nebr.ireg[MINUS]->cells.begin(), nebr.ireg[MINUS]->cells.end(), [&](const RegCellId& id1, const RegCellId& id2)
				{
					if (id1.cell == cell.id || id1.cell == nebr.nebr.cell)
						return false;
					if (id2.cell == cell.id || id2.cell == nebr.nebr.cell)
						return true;
					else
						return point::distance(nebr.cent, cells[id1.cell].cent) < point::distance(nebr.cent, cells[id2.cell].cent);
				})->cell];
				const auto& cell2 = cells[min_element(nebr.ireg[PLUS]->cells.begin(), nebr.ireg[PLUS]->cells.end(), [&](const RegCellId& id1, const RegCellId& id2)
				{
					if (id1.cell == cell.id || id1.cell == nebr.nebr.cell)
						return false;
					if (id2.cell == cell.id || id2.cell == nebr.nebr.cell)
						return true;
					else
						return point::distance(nebr.cent, cells[id1.cell].cent) < point::distance(nebr.cent, cells[id2.cell].cent); })->cell];

				nebr.nearest_cells[2] = point::distance(nebr.cent, cell1.cent) < point::distance(nebr.cent, cell2.cent) ? cell1.id : cell2.id;
				nebr.nearest_dists[2] = point::distance(nebr.cent, cells[nebr.nearest_cells[2]].cent);
			}
	}
}
void Mesh::set_geom_props()
{
	double tmp1, tmp2;

	for (auto& el : cells)
	{
		if (el.type == elem::HEX || el.type == elem::BORDER_HEX)
		{
			// element center
			tmp1 = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]]);
			tmp2 = square(pts[el.verts[2]], pts[el.verts[3]], pts[el.verts[0]]);
			el.nebrs[0].cent = (tmp1 * (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]]) + tmp2 * (pts[el.verts[2]] + pts[el.verts[3]] + pts[el.verts[0]])) / (tmp1 + tmp2) / 3.0;
			tmp1 = square(pts[el.verts[4]], pts[el.verts[5]], pts[el.verts[6]]);
			tmp2 = square(pts[el.verts[6]], pts[el.verts[7]], pts[el.verts[4]]);
			el.nebrs[1].cent = (tmp1 * (pts[el.verts[4]] + pts[el.verts[5]] + pts[el.verts[6]]) + tmp2 * (pts[el.verts[6]] + pts[el.verts[7]] + pts[el.verts[4]])) / (tmp1 + tmp2) / 3.0;
			el.cent = (el.nebrs[0].cent + el.nebrs[1].cent) / 2.0;

			// bottom
			el.nebrs[0].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[3]]);
			el.nebrs[0].L = distance(el.cent, el.nebrs[0].cent);
			// top
			el.nebrs[1].S = square(pts[el.verts[4]], pts[el.verts[5]], pts[el.verts[6]], pts[el.verts[7]]);
			el.nebrs[1].L = distance(el.cent, el.nebrs[1].cent);

			el.nebrs[2].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[5]], pts[el.verts[4]]);
			el.nebrs[2].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[5]] + pts[el.verts[4]]) / 4.0;
			el.nebrs[2].L = distance(el.cent, el.nebrs[2].cent);
			el.nebrs[2].n = vector_product(pts[el.verts[1]] - pts[el.verts[0]], pts[el.verts[4]] - pts[el.verts[0]]);
			el.nebrs[2].nu = { (el.nebrs[2].cent - el.cent).y, -(el.nebrs[2].cent - el.cent).x, 0.0 };

			el.nebrs[3].S = square(pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[6]], pts[el.verts[5]]);
			el.nebrs[3].cent = (pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[6]] + pts[el.verts[5]]) / 4.0;
			el.nebrs[3].L = distance(el.cent, el.nebrs[3].cent);
			el.nebrs[3].n = vector_product(pts[el.verts[2]] - pts[el.verts[1]], pts[el.verts[5]] - pts[el.verts[1]]);
			el.nebrs[3].nu = { (el.nebrs[3].cent - el.cent).y, -(el.nebrs[3].cent - el.cent).x, 0.0 };

			el.nebrs[4].S = square(pts[el.verts[2]], pts[el.verts[3]], pts[el.verts[7]], pts[el.verts[6]]);
			el.nebrs[4].cent = (pts[el.verts[2]] + pts[el.verts[3]] + pts[el.verts[7]] + pts[el.verts[6]]) / 4.0;
			el.nebrs[4].L = distance(el.cent, el.nebrs[4].cent);
			el.nebrs[4].n = vector_product(pts[el.verts[3]] - pts[el.verts[2]], pts[el.verts[6]] - pts[el.verts[2]]);
			el.nebrs[4].nu = { (el.nebrs[4].cent - el.cent).y, -(el.nebrs[4].cent - el.cent).x, 0.0 };

			el.nebrs[5].S = square(pts[el.verts[3]], pts[el.verts[0]], pts[el.verts[4]], pts[el.verts[7]]);
			el.nebrs[5].cent = (pts[el.verts[3]] + pts[el.verts[0]] + pts[el.verts[4]] + pts[el.verts[7]]) / 4.0;
			el.nebrs[5].L = distance(el.cent, el.nebrs[5].cent);
			el.nebrs[5].n = vector_product(pts[el.verts[0]] - pts[el.verts[3]], pts[el.verts[7]] - pts[el.verts[3]]);
			el.nebrs[5].nu = { (el.nebrs[5].cent - el.cent).y, -(el.nebrs[5].cent - el.cent).x, 0.0 };

			el.T[0] = 2.0 * square(el.cent, el.nebrs[2].cent, el.nebrs[3].cent);
			el.T[1] = 2.0 * square(el.cent, el.nebrs[3].cent, el.nebrs[4].cent);
			el.T[2] = 2.0 * square(el.cent, el.nebrs[4].cent, el.nebrs[5].cent);
			el.T[3] = 2.0 * square(el.cent, el.nebrs[5].cent, el.nebrs[2].cent);
			el.V = fabs(pts[el.verts[0]].z - pts[el.verts[4]].z) * (el.nebrs[0].S + el.nebrs[1].S) / 2.0;
		}
		else if (el.type == elem::PRISM)
		{
			el.cent = { 0, 0, 0 };
			for (int i = 0; i < el.verts_num; i++)
				el.cent += pts[el.verts[i]];
			el.cent /= el.verts_num;

			// top
			el.nebrs[0].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]]);
			el.nebrs[0].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]]) / 3.0;
			el.nebrs[0].L = distance(el.cent, el.nebrs[0].cent);
			// bottom
			el.nebrs[1].S = square(pts[el.verts[3]], pts[el.verts[4]], pts[el.verts[5]]);
			el.nebrs[1].cent = (pts[el.verts[3]] + pts[el.verts[4]] + pts[el.verts[5]]) / 3.0;
			el.nebrs[1].L = distance(el.cent, el.nebrs[1].cent);

			el.nebrs[2].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[4]], pts[el.verts[3]]);
			el.nebrs[2].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[4]] + pts[el.verts[3]]) / 4.0;
			el.nebrs[2].L = distance(el.cent, el.nebrs[2].cent);
			el.nebrs[2].n = vector_product(pts[el.verts[1]] - pts[el.verts[0]], pts[el.verts[3]] - pts[el.verts[0]]) / 2.0;
			el.nebrs[2].nu = { (el.nebrs[2].cent - el.cent).y, -(el.nebrs[2].cent - el.cent).x, 0.0 };


			el.nebrs[3].S = square(pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[5]], pts[el.verts[4]]);
			el.nebrs[3].cent = (pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[5]] + pts[el.verts[4]]) / 4.0;
			el.nebrs[3].L = distance(el.cent, el.nebrs[3].cent);
			el.nebrs[3].n = vector_product(pts[el.verts[2]] - pts[el.verts[1]], pts[el.verts[4]] - pts[el.verts[1]]) / 2.0;
			el.nebrs[3].nu = { (el.nebrs[3].cent - el.cent).y, -(el.nebrs[3].cent - el.cent).x, 0.0 };

			el.nebrs[4].S = square(pts[el.verts[2]], pts[el.verts[0]], pts[el.verts[3]], pts[el.verts[5]]);
			el.nebrs[4].cent = (pts[el.verts[2]] + pts[el.verts[0]] + pts[el.verts[3]] + pts[el.verts[5]]) / 4.0;
			el.nebrs[4].L = distance(el.cent, el.nebrs[4].cent);
			el.nebrs[4].n = vector_product(pts[el.verts[0]] - pts[el.verts[2]], pts[el.verts[5]] - pts[el.verts[2]]) / 2.0;
			el.nebrs[4].nu = { (el.nebrs[4].cent - el.cent).y, -(el.nebrs[4].cent - el.cent).x, 0.0 };

			el.T[0] = 2.0 * square(el.cent, el.nebrs[2].cent, el.nebrs[3].cent);
			el.T[1] = 2.0 * square(el.cent, el.nebrs[3].cent, el.nebrs[4].cent);
			el.T[2] = 2.0 * square(el.cent, el.nebrs[4].cent, el.nebrs[2].cent);
			el.V = fabs(pts[el.verts[0]].z - pts[el.verts[3]].z) * (el.nebrs[0].S + el.nebrs[1].S) / 2.0;
		}
		else if (el.type == elem::BORDER_TRI)
		{
			el.nebrs[0].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]]);
			el.cent = el.nebrs[0].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]]) / 3.0;
			el.nebrs[0].L = el.V = 0.0;
		}
		else if (el.type == elem::BORDER_QUAD)
		{
			el.nebrs[0].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[3]]);
			tmp1 = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]]);
			tmp2 = square(pts[el.verts[2]], pts[el.verts[3]], pts[el.verts[0]]);
			el.cent  = el.nebrs[0].cent = (tmp1 * (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]]) + tmp2 * (pts[el.verts[2]] + pts[el.verts[3]] + pts[el.verts[0]])) / (tmp1 + tmp2) / 3.0;
			el.nebrs[0].L = el.V = 0.0;
			el.nebrs[0].n = vector_product(pts[el.verts[1]] - pts[el.verts[0]], pts[el.verts[3]] - pts[el.verts[0]]);
		}
		else if (el.type == elem::FRAC_QUAD)
		{
			el.nebrs[0].S = el.nebrs[1].S = square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[3]]);
			el.cent = el.nebrs[0].cent = el.nebrs[1].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[3]]) / 4.0;
			el.nebrs[0].L = el.nebrs[1].L = el.V = 0.0;
			el.nebrs[0].n = vector_product(pts[el.verts[1]] - pts[el.verts[0]], pts[el.verts[3]] - pts[el.verts[0]]);
		}
	}

	// Cell Z-orientation transparency checking
	for (const auto& cell : cells)
	{
		if (cell.type == elem::HEX || cell.type == elem::BORDER_HEX || cell.type == elem::PRISM)
		{
			for (int i = 0; i < cell.verts_num / 2; i++)
			{
				assert(pts[cell.verts[i]].z == pts[cell.verts[0]].z);
				assert(pts[cell.verts[i]].z < pts[cell.verts[i + cell.verts_num / 2]].z);
			}
			for (int i = cell.verts_num / 2; i < cell.verts_num; i++)
				assert(pts[cell.verts[i]].z == pts[cell.verts[cell.verts_num / 2]].z);
		}
	}
};
void Mesh::set_interaction_regions()
{
	// Getting max z-axis coordinate
	z_max = max_element(pts.begin(), pts.end(), [&](const Point& pt1, const Point& pt2) {return pt1.z < pt2.z; })->z;
	
	bool isBorder = true;
	for (auto& pt : pts)
	{
		isBorder = true;
		// Check is not border
		for (const auto& cell : pt.cells)
			if (cells[cell].type != elem::BORDER_HEX)
			{
				isBorder = false;
				break;
			}
		if (isBorder)
			pt.type = point::BORDER;
		// All except top z-plane
		if (!isBorder && fabs(pt.z - z_max) > EQUALITY_TOLERANCE)
		{
			pt.int_reg = new Interaction(pt.cells);
			auto& cur_cells = pt.int_reg->cells;
			cur_cells.erase(remove_if(cur_cells.begin(), cur_cells.end(), [&](RegCellId cell_id) { return (cells[cell_id.cell].cent.z - pt.z < 0.0) ? true : false; }), cur_cells.end());
		}
	}
	
	// Order the cells in interaction regions
	for (auto& pt : pts)
	{
		if (pt.int_reg != nullptr)
		{
			auto& ireg_cells = pt.int_reg->cells;
			sort(ireg_cells.begin(), ireg_cells.end(), [&](RegCellId& id1, RegCellId& id2)
			{
				double phi1, phi2;
				const Cell& cell1 = cells[id1.cell];	const Cell& cell2 = cells[id2.cell];
				const double x1 = cell1.cent.x - pt.x;	const double y1 = cell1.cent.y - pt.y;
				const double x2 = cell2.cent.x - pt.x;	const double y2 = cell2.cent.y - pt.y;
				if (x1 < 0.0)
					phi1 = M_PI + atan(y1 / x1);
				else
					phi1 = (y1 < 0.0 ? 2.0 * M_PI + atan(y1 / x1) : atan(y1 / x1));
				if (x2 < 0.0)
					phi2 = M_PI + atan(y2 / x2);
				else
					phi2 = (y2 < 0.0 ? 2.0 * M_PI + atan(y2 / x2) : atan(y2 / x2));

				return phi1 < phi2;
			});
		}
	}

	// Fill pointers to interaction regions in cells neighbours
	using point::HalfType;
	vector<RegCellId>::iterator it;
	for (int i = 0; i < inner_size; i++)
	{
		Cell& cell = cells[i];
		if (cell.type == elem::HEX || cell.type == elem::BORDER_HEX)
		{
			// 1 or 5
			const auto& pt1 = pts[cell.verts[1]];
			if (pt1.int_reg != nullptr)
			{
				auto& cells1_plus = pt1.int_reg->cells;
				it = find(cells1_plus.begin(), cells1_plus.end(), cell.id);
				assert(it != cells1_plus.end());
				cell.nebrs[2].ireg[HalfType::PLUS] = cell.nebrs[3].ireg[HalfType::MINUS] = pt1.int_reg;
				if (vector_product(cell.nebrs[3].cent - pt1, cell.nebrs[2].cent - pt1).z > 0.0)
				{	it->nebr[0] = 3;	it->nebr[1] = 2; }
				else
				{	it->nebr[1] = 3;	it->nebr[0] = 2; }
			}
			// 0 or 4
			const auto& pt2 = pts[cell.verts[0]];
			if (pt2.int_reg != nullptr)
			{
				auto& cells1_minus = pt2.int_reg->cells;
				it = find(cells1_minus.begin(), cells1_minus.end(), cell.id);
				assert(it != cells1_minus.end());
				cell.nebrs[2].ireg[HalfType::MINUS] = cell.nebrs[5].ireg[HalfType::PLUS] = pt2.int_reg;
				if (vector_product(cell.nebrs[5].cent - pt2, cell.nebrs[2].cent - pt2).z > 0.0)
				{	it->nebr[0] = 5;	it->nebr[1] = 2; }
				else
				{	it->nebr[1] = 5;	it->nebr[0] = 2; }
			}
			// 2 or 6
			const auto& pt3 = pts[cell.verts[2]];
			if (pt3.int_reg != nullptr)
			{
				auto& cells2_plus = pt3.int_reg->cells;
				it = find(cells2_plus.begin(), cells2_plus.end(), cell.id);
				assert(it != cells2_plus.end());
				cell.nebrs[3].ireg[HalfType::PLUS] = cell.nebrs[4].ireg[HalfType::MINUS] = pt3.int_reg;
				if (vector_product(cell.nebrs[4].cent - pt3, cell.nebrs[3].cent - pt3).z > 0.0)
				{	it->nebr[0] = 4;	it->nebr[1] = 3; }
				else
				{	it->nebr[1] = 4;	it->nebr[0] = 3; }
			}
			// 3 or 7
			const auto& pt4 = pts[cell.verts[3]];
			if (pt4.int_reg != nullptr)
			{
				auto& cells3_plus = pt4.int_reg->cells;
				it = find(cells3_plus.begin(), cells3_plus.end(), cell.id);
				assert(it != cells3_plus.end());
				cell.nebrs[4].ireg[HalfType::PLUS] = cell.nebrs[5].ireg[HalfType::MINUS] = pt4.int_reg;
				if (vector_product(cell.nebrs[5].cent - pt4, cell.nebrs[4].cent - pt4).z > 0.0)
				{	it->nebr[0] = 5;	it->nebr[1] = 4; }
				else
				{	it->nebr[1] = 5;	it->nebr[0] = 4; }
			}
		}
		else if (cell.type == elem::PRISM)
		{
			// 1 or 4
			const auto& pt1 = pts[cell.verts[1]];
			if (pt1.int_reg != nullptr)
			{
				auto& cells1_plus = pt1.int_reg->cells;
				it = find(cells1_plus.begin(), cells1_plus.end(), cell.id);
				assert(it != cells1_plus.end());
				cell.nebrs[2].ireg[HalfType::PLUS] = cell.nebrs[3].ireg[HalfType::MINUS] = pt1.int_reg;
				if (vector_product(cell.nebrs[3].cent - pt1, cell.nebrs[2].cent - pt1).z > 0.0)
				{	it->nebr[0] = 3;	it->nebr[1] = 2; }
				else
				{	it->nebr[1] = 3;	it->nebr[0] = 2; }
			}
			// 0 or 3
			const auto& pt2 = pts[cell.verts[0]];
			if (pt2.int_reg != nullptr)
			{
				auto& cells1_minus = pt2.int_reg->cells;
				it = find(cells1_minus.begin(), cells1_minus.end(), cell.id);
				assert(it != cells1_minus.end());
				cell.nebrs[2].ireg[HalfType::MINUS] = cell.nebrs[4].ireg[HalfType::PLUS] = pt2.int_reg;
				if (vector_product(cell.nebrs[4].cent - pt2, cell.nebrs[2].cent - pt2).z > 0.0)
				{	it->nebr[0] = 4;	it->nebr[1] = 2; }
				else
				{	it->nebr[1] = 4;	it->nebr[0] = 2; }
			}
			// 2 or 5
			const auto& pt3 = pts[cell.verts[2]];
			if (pt3.int_reg != nullptr)
			{
				auto& cells2_plus = pt3.int_reg->cells;
				it = find(cells2_plus.begin(), cells2_plus.end(), cell.id);
				assert(it != cells2_plus.end());
				cell.nebrs[3].ireg[HalfType::PLUS] = cell.nebrs[4].ireg[HalfType::MINUS] = pt3.int_reg;
				if (vector_product(cell.nebrs[4].cent - pt3, cell.nebrs[3].cent - pt3).z > 0.0)
				{	it->nebr[0] = 4;	it->nebr[1] = 3; }
				else
				{	it->nebr[1] = 4;	it->nebr[0] = 3; }
			}
		}	
	}
}
void Mesh::calc_transmissibilities()
{
	typedef array<array<adouble, 2>, 2> Interface;
	
	double sum, abs_sum;
	adouble buf, cond;
	//LocalMatrix<double> a, b, c, d, tmp, t;
	/*double buf_a [MAX_REG_CELLS * MAX_REG_CELLS], buf_b[MAX_REG_CELLS * MAX_REG_CELLS], buf_c[MAX_REG_CELLS * MAX_REG_CELLS], buf_d[MAX_REG_CELLS * MAX_REG_CELLS];
	int ind_ai [3 * MAX_REG_CELLS], ind_aj [3 * MAX_REG_CELLS], ind_bi[3 * MAX_REG_CELLS], ind_bj[3 * MAX_REG_CELLS], 
			ind_ci[3 * MAX_REG_CELLS], ind_cj[3 * MAX_REG_CELLS], ind_di[3 * MAX_REG_CELLS], ind_dj[3 * MAX_REG_CELLS];*/
	int i_plus, i_minus;

	for (auto& pt : pts)
	{
		if (pt.int_reg != nullptr)
		{
			auto& reg = *pt.int_reg;

			if (reg.trans.size() == 0)
			{
				reg.trans.resize(reg.cells.size());
				for (auto& tr : reg.trans)
					tr.resize(reg.cells.size());

				reg.trans_doub.resize(reg.cells.size());
				for (auto& tr : reg.trans_doub)
					tr.resize(reg.cells.size());
			}

			vector<Interface> omega(reg.cells.size());
			for(int i = 0; i < reg.cells.size(); i++)
			{
				const auto& cell_id1 = reg.cells[i];
				const auto& cell1 = cells[cell_id1.cell];
				const auto& nebr11 = cell1.nebrs[cell_id1.nebr[0]];
				const auto& nebr12 = cell1.nebrs[cell_id1.nebr[1]];
				const auto perm1 = get_XY_perm(cell1);

				const auto& cell_id2 = (i < reg.cells.size() - 1 ? reg.cells[i + 1] : reg.cells[0]);
				const auto& cell2 = cells[cell_id2.cell];
				const auto& nebr21 = cell2.nebrs[cell_id2.nebr[0]];
				const auto& nebr22 = cell2.nebrs[cell_id2.nebr[1]];
				const auto perm2 = get_XY_perm(cell2);

				const auto& nu11 = nebr11.nu;	const auto& nu12 = nebr12.nu;
				const auto& nu21 = nebr21.nu;	const auto& nu22 = nebr22.nu;
				omega[i][0][0] = -(nebr12.n.x * perm1.kx * nu12.x + nebr12.n.y * perm1.ky * nu12.y) / cell1.T[(int)(cell_id1.nebr[1]) - 2];
				omega[i][0][1] = -(nebr12.n.x * perm1.kx * nu11.x + nebr12.n.y * perm1.ky * nu11.y) / cell1.T[(int)(cell_id1.nebr[1]) - 2];

				omega[i][1][0] = -(nebr21.n.x * perm2.kx * nu22.x + nebr21.n.y * perm2.ky * nu22.y) / cell2.T[(int)(cell_id2.nebr[0]) - 2];
				omega[i][1][1] = -(nebr21.n.x * perm2.kx * nu21.x + nebr21.n.y * perm2.ky * nu21.y) / cell2.T[(int)(cell_id2.nebr[0]) - 2];
			}

			for (int i = 0; i < reg.cells.size(); i++)
				for (int j = 0; j < reg.cells.size(); j++)
					buf_a[i][j] = buf_b[i][j] = buf_c[i][j] = buf_d[i][j] = 0.0;
			for (int i = 0; i < reg.cells.size(); i++)
			{
				i_plus = (i + 1 ? i < reg.cells.size() - 1 : 0);
				i_minus = (i - 1 ? i > 0 : reg.cells.size() - 1);
				// A
				buf_a[i][i_minus] =		omega[i][0][1];
				buf_a[i][i] =			omega[i][0][0] - omega[i][1][1];
				buf_a[i][i_plus] =		-omega[i][1][0];
				// B
				buf_b[i][i] =			omega[i][0][0] + omega[i][0][1];
				buf_b[i][i_plus] =		-omega[i][1][0] - omega[i][1][1];
				// C
				buf_c[i][i] =			omega[i][0][0];
				buf_c[i][i_minus] =		omega[i][0][1];
				// D
				buf_d[i][i] =			omega[i][0][0] + omega[i][0][1];
			}
			QR_Decomp<adouble> qr;
			inv_a = qr.Invert(buf_a, reg.cells.size());
			multiply<adouble>(buf_c, inv_a, mult, reg.cells.size());
			for (int i = 0; i < reg.cells.size(); i++)
				delete[] inv_a[i];
			delete[] inv_a;
			multiply<adouble>(mult, buf_b, t, reg.cells.size());
			matrix_add<adouble>(t, buf_d, -1.0, reg.cells.size());

			// Sum == 0
			for (int i = 0; i < reg.cells.size(); i++)
			{
				abs_sum = sum = 0.0;
				for (int j = 0; j < reg.cells.size(); j++)
				{
					buf = t[i][j];
					sum += buf.value();		abs_sum += fabs(buf.value());
				}
				/*for (int j = 0; j < reg.cells.size(); j++)
				{
					cond = t[i][j] / abs_sum < 10.0 * REL_TOL ? true : false;
					condassign(t[i][j], cond, (adouble)0.0);
				}*/
				//assert(fabs(sum) / abs_sum < REL_TOL);

				for (int j = 0; j < reg.cells.size(); j++)
				{
					reg.trans[i][j] = t[i][j];
					reg.trans_doub[i][j] = t[i][j].value();
				}
			}
		}
	}
}
size_t Mesh::getCellsSize() const
{
	return cells.size();
}
void Mesh::setNebrId()
{
	for (auto& el : cells)
		for (char i = 0; i < el.nebrs_num; i++)
			el.nebrs[i].nebr.nebr = findNebrId(el, cells[el.nebrs[i].nebr.cell]);
}