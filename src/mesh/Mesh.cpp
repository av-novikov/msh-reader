#include "src/mesh/Mesh.hpp"
#include <algorithm>
#include <unordered_set>

using namespace grid;
using namespace point;
using namespace std;

Mesh::Mesh()
{
}
Mesh::~Mesh()
{
}
void Mesh::process_geometry()
{
	setNebrId();
	set_geom_props();
	count_types();
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
	int cur_id, prev_id, next_id, counter;

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
};
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