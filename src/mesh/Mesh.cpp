#include "src/mesh/Mesh.hpp"
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <utility>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace grid;
using namespace point;
using namespace std;

Mesh::Mesh()
{
}
Mesh::~Mesh()
{
	for (auto& pt : pts)
		delete pt.int_reg;
}
void Mesh::process_geometry()
{
	setNebrId();
	set_geom_props();
	count_types();

	set_interaction_regions();
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
			el.T[0] = 2.0 * square(el.cent, el.nebrs[2].cent, el.nebrs[3].cent);

			el.nebrs[3].S = square(pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[6]], pts[el.verts[5]]);
			el.nebrs[3].cent = (pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[6]] + pts[el.verts[5]]) / 4.0;
			el.nebrs[3].L = distance(el.cent, el.nebrs[3].cent);
			el.nebrs[3].n = vector_product(pts[el.verts[2]] - pts[el.verts[1]], pts[el.verts[5]] - pts[el.verts[1]]);
			el.nebrs[3].nu = { (el.nebrs[3].cent - el.cent).y, -(el.nebrs[3].cent - el.cent).x, 0.0 };
			el.T[1] = 2.0 * square(el.cent, el.nebrs[3].cent, el.nebrs[4].cent);

			el.nebrs[4].S = square(pts[el.verts[2]], pts[el.verts[3]], pts[el.verts[7]], pts[el.verts[6]]);
			el.nebrs[4].cent = (pts[el.verts[2]] + pts[el.verts[3]] + pts[el.verts[7]] + pts[el.verts[6]]) / 4.0;
			el.nebrs[4].L = distance(el.cent, el.nebrs[4].cent);
			el.nebrs[4].n = vector_product(pts[el.verts[3]] - pts[el.verts[2]], pts[el.verts[6]] - pts[el.verts[2]]);
			el.nebrs[4].nu = { (el.nebrs[4].cent - el.cent).y, -(el.nebrs[4].cent - el.cent).x, 0.0 };
			el.T[2] = 2.0 * square(el.cent, el.nebrs[4].cent, el.nebrs[5].cent);

			el.nebrs[5].S = square(pts[el.verts[3]], pts[el.verts[0]], pts[el.verts[4]], pts[el.verts[7]]);
			el.nebrs[5].cent = (pts[el.verts[3]] + pts[el.verts[0]] + pts[el.verts[4]] + pts[el.verts[7]]) / 4.0;
			el.nebrs[5].L = distance(el.cent, el.nebrs[5].cent);
			el.nebrs[5].n = vector_product(pts[el.verts[0]] - pts[el.verts[3]], pts[el.verts[7]] - pts[el.verts[3]]);
			el.nebrs[5].nu = { (el.nebrs[5].cent - el.cent).y, -(el.nebrs[5].cent - el.cent).x, 0.0 };
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
			el.T[0] = 2.0 * square(el.cent, el.nebrs[2].cent, el.nebrs[3].cent);

			el.nebrs[3].S = square(pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[5]], pts[el.verts[4]]);
			el.nebrs[3].cent = (pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[5]] + pts[el.verts[4]]) / 4.0;
			el.nebrs[3].L = distance(el.cent, el.nebrs[3].cent);
			el.nebrs[3].n = vector_product(pts[el.verts[2]] - pts[el.verts[1]], pts[el.verts[4]] - pts[el.verts[1]]) / 2.0;
			el.nebrs[3].nu = { (el.nebrs[3].cent - el.cent).y, -(el.nebrs[3].cent - el.cent).x, 0.0 };
			el.T[1] = 2.0 * square(el.cent, el.nebrs[3].cent, el.nebrs[4].cent);

			el.nebrs[4].S = square(pts[el.verts[2]], pts[el.verts[0]], pts[el.verts[3]], pts[el.verts[5]]);
			el.nebrs[4].cent = (pts[el.verts[2]] + pts[el.verts[0]] + pts[el.verts[3]] + pts[el.verts[5]]) / 4.0;
			el.nebrs[4].L = distance(el.cent, el.nebrs[4].cent);
			el.nebrs[4].n = vector_product(pts[el.verts[0]] - pts[el.verts[2]], pts[el.verts[5]] - pts[el.verts[2]]) / 2.0;
			el.nebrs[4].nu = { (el.nebrs[4].cent - el.cent).y, -(el.nebrs[4].cent - el.cent).x, 0.0 };
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
	Point* vert;
	auto get_phi = [](const Point& vert, const Point& nebr_cent) -> double
	{
		double phi;
		const double x = nebr_cent.x - vert.x;	const double y = nebr_cent.y - vert.y;
		if (x < 0.0)
			phi = M_PI + atan(y / x);
		else
			phi = (y < 0.0 ? 2.0 * M_PI + atan(y / x) : atan(y / x));
		return phi;
	};

	for (int i = 0; i < inner_size; i++)
	{
		Cell& cell = cells[i];
		if (cell.type == elem::HEX || cell.type == elem::BORDER_HEX)
		{
			// 1 or 5
			auto& cells1_plus = pts[cell.verts[1]].int_reg->cells;
			it = find(cells1_plus.begin(), cells1_plus.end(), cell.id);
			if (it != cells1_plus.end())
				vert = &pts[cell.verts[1]];
			else
				vert = &pts[cell.verts[5]];

			cell.nebrs[2].ireg[HalfType::PLUS] = cell.nebrs[3].ireg[HalfType::MINUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[2].cent) > get_phi(*vert, cell.nebrs[3].cent))
			{	it->nebr[0] = 3;	it->nebr[1] = 2;	}
			else
			{	it->nebr[1] = 3;	it->nebr[0] = 2;	}

			// 0 or 4
			auto& cells1_minus = pts[cell.verts[0]].int_reg->cells;
			it = find(cells1_minus.begin(), cells1_minus.end(), cell.id);
			if (it != cells1_minus.end())
				vert = &pts[cell.verts[0]];
			else
				vert = &pts[cell.verts[4]];

			cell.nebrs[2].ireg[HalfType::MINUS] = cell.nebrs[5].ireg[HalfType::PLUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[2].cent) > get_phi(*vert, cell.nebrs[5].cent))
			{	it->nebr[0] = 5;	it->nebr[1] = 2;	}
			else
			{	it->nebr[1] = 5;	it->nebr[0] = 2;	}

			// 2 or 6
			auto& cells2_plus = pts[cell.verts[2]].int_reg->cells;
			it = find(cells2_plus.begin(), cells2_plus.end(), cell.id);
			if (it != cells2_plus.end())
				vert = &pts[cell.verts[2]];
			else
				vert = &pts[cell.verts[6]];

			cell.nebrs[3].ireg[HalfType::PLUS] = cell.nebrs[4].ireg[HalfType::MINUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[3].cent) > get_phi(*vert, cell.nebrs[4].cent))
			{	it->nebr[0] = 4;	it->nebr[1] = 3;	}
			else
			{	it->nebr[1] = 4;	it->nebr[0] = 3;	}

			// 3 or 7
			auto& cells3_plus = pts[cell.verts[3]].int_reg->cells;
			it = find(cells3_plus.begin(), cells3_plus.end(), cell.id);
			if (it != cells3_plus.end())
				vert = &pts[cell.verts[3]];
			else
				vert = &pts[cell.verts[7]];

			cell.nebrs[4].ireg[HalfType::PLUS] = cell.nebrs[5].ireg[HalfType::MINUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[4].cent) > get_phi(*vert, cell.nebrs[5].cent))
			{	it->nebr[0] = 5;	it->nebr[1] = 4;	}
			else
			{	it->nebr[1] = 5;	it->nebr[0] = 4;	}
		}
		else if (cell.type == elem::PRISM)
		{
			// 1 or 4
			auto& cells1_plus = pts[cell.verts[1]].int_reg->cells;
			it = find(cells1_plus.begin(), cells1_plus.end(), cell.id);
			if (it != cells1_plus.end())
				vert = &pts[cell.verts[1]];
			else
				vert = &pts[cell.verts[4]];

			cell.nebrs[2].ireg[HalfType::PLUS] = cell.nebrs[3].ireg[HalfType::MINUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[2].cent) > get_phi(*vert, cell.nebrs[3].cent))
			{	it->nebr[0] = 3;	it->nebr[1] = 2;	}
			else
			{	it->nebr[1] = 3;	it->nebr[0] = 2;	}

			// 0 or 3
			auto& cells1_minus = pts[cell.verts[0]].int_reg->cells;
			it = find(cells1_minus.begin(), cells1_minus.end(), cell.id);
			if (it != cells1_minus.end())
				vert = &pts[cell.verts[0]];
			else
				vert = &pts[cell.verts[3]];

			cell.nebrs[2].ireg[HalfType::MINUS] = cell.nebrs[4].ireg[HalfType::PLUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[2].cent) > get_phi(*vert, cell.nebrs[4].cent))
			{	it->nebr[0] = 4;	it->nebr[1] = 2;	}
			else
			{	it->nebr[1] = 4;	it->nebr[0] = 2;	}
			
			// 2 or 5
			auto& cells2_plus = pts[cell.verts[2]].int_reg->cells;
			it = find(cells2_plus.begin(), cells2_plus.end(), cell.id);
			if (it != cells2_plus.end())
				vert = &pts[cell.verts[2]];
			else
				vert = &pts[cell.verts[5]];
			cell.nebrs[3].ireg[HalfType::PLUS] = cell.nebrs[4].ireg[HalfType::MINUS] = vert->int_reg;
			if (get_phi(*vert, cell.nebrs[3].cent) > get_phi(*vert, cell.nebrs[4].cent))
			{	it->nebr[0] = 4;	it->nebr[1] = 3;	}
			else
			{	it->nebr[1] = 4;	it->nebr[0] = 3;	}
		}	
	}
}
void Mesh::calc_transmissibilities()
{
	typedef array<array<double, 2>, 2> Interface;
	for (auto& pt : pts)
	{
		if (pt.type == point::INNER && fabs(pt.z - z_max) > EQUALITY_TOLERANCE)
		{
			auto& reg = *pt.int_reg;
			vector<Interface> omega(reg.cells.size());
			for(const auto& cell_id : reg.cells)
			{
				const auto& cell = cells[cell_id.cell];
				const auto& nebr1 = cells[cell_id.cell].nebrs[cell_id.nebr[0]];
				const auto& nebr2 = cells[cell_id.cell].nebrs[cell_id.nebr[1]];
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