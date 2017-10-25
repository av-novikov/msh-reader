#include "src/mesh/Mesh.hpp"
#include <algorithm>

using namespace grid;
using namespace std;

Mesh::Mesh()
{
}
Mesh::~Mesh()
{
}
void Mesh::process_geometry()
{
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
		if (el.type == elem::HEX || el.type == elem::PRISM)	inner_size++;
		else if (el.type == elem::FRAC_QUAD)				frac_size++;
		else												border_size++;
	}
}
void Mesh::set_geom_props()
{
	for (auto& el : cells)
	{
		// element center
		el.cent = { 0.0, 0.0, 0.0 };
		for (int i = 0; i < el.verts_num; i++)
			el.cent += pts[el.verts[i]];
		el.cent /= el.verts_num;

		if (el.type == elem::EType::HEX)
		{
			el.nebrs[0].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[3]]);
			el.nebrs[0].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[3]]) / 4.0;
			el.nebrs[0].L = point::distance(el.cent, el.nebrs[0].cent);

			el.nebrs[1].S = point::square(pts[el.verts[4]], pts[el.verts[5]], pts[el.verts[6]], pts[el.verts[7]]);
			el.nebrs[1].cent = (pts[el.verts[4]] + pts[el.verts[5]] + pts[el.verts[6]] + pts[el.verts[7]]) / 4.0;
			el.nebrs[1].L = point::distance(el.cent, el.nebrs[1].cent);

			el.nebrs[2].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[5]], pts[el.verts[4]]);
			el.nebrs[2].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[5]] + pts[el.verts[4]]) / 4.0;
			el.nebrs[2].L = point::distance(el.cent, el.nebrs[2].cent);

			el.nebrs[3].S = point::square(pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[6]], pts[el.verts[5]]);
			el.nebrs[3].cent = (pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[6]] + pts[el.verts[5]]) / 4.0;
			el.nebrs[3].L = point::distance(el.cent, el.nebrs[3].cent);

			el.nebrs[4].S = point::square(pts[el.verts[2]], pts[el.verts[3]], pts[el.verts[7]], pts[el.verts[6]]);
			el.nebrs[4].cent = (pts[el.verts[2]] + pts[el.verts[3]] + pts[el.verts[7]] + pts[el.verts[6]]) / 4.0;
			el.nebrs[4].L = point::distance(el.cent, el.nebrs[4].cent);

			el.nebrs[5].S = point::square(pts[el.verts[3]], pts[el.verts[0]], pts[el.verts[4]], pts[el.verts[7]]);
			el.nebrs[5].cent = (pts[el.verts[3]] + pts[el.verts[0]] + pts[el.verts[4]] + pts[el.verts[7]]) / 4.0;
			el.nebrs[5].L = point::distance(el.cent, el.nebrs[5].cent);

			el.V = fabs(pts[el.verts[0]].z - pts[el.verts[4]].z) * (el.nebrs[0].S + el.nebrs[1].S) / 2.0;
		}
		else if (el.type == elem::EType::PRISM)
		{
			el.nebrs[0].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]]);
			el.nebrs[0].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]]) / 3.0;
			el.nebrs[0].L = point::distance(el.cent, el.nebrs[0].cent);

			el.nebrs[1].S = point::square(pts[el.verts[3]], pts[el.verts[4]], pts[el.verts[5]]);
			el.nebrs[1].cent = (pts[el.verts[3]] + pts[el.verts[4]] + pts[el.verts[5]]) / 3.0;
			el.nebrs[1].L = point::distance(el.cent, el.nebrs[1].cent);

			el.nebrs[2].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[4]], pts[el.verts[3]]);
			el.nebrs[2].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[4]] + pts[el.verts[3]]) / 4.0;
			el.nebrs[2].L = point::distance(el.cent, el.nebrs[2].cent);

			el.nebrs[3].S = point::square(pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[5]], pts[el.verts[4]]);
			el.nebrs[3].cent = (pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[5]] + pts[el.verts[4]]) / 4.0;
			el.nebrs[3].L = point::distance(el.cent, el.nebrs[3].cent);

			el.nebrs[4].S = point::square(pts[el.verts[2]], pts[el.verts[0]], pts[el.verts[3]], pts[el.verts[5]]);
			el.nebrs[4].cent = (pts[el.verts[2]] + pts[el.verts[0]] + pts[el.verts[3]] + pts[el.verts[5]]) / 4.0;
			el.nebrs[4].L = point::distance(el.cent, el.nebrs[4].cent);

			el.V = fabs(pts[el.verts[0]].z - pts[el.verts[3]].z) * (el.nebrs[0].S + el.nebrs[1].S) / 2.0;
		}
		else if (el.type == elem::EType::BORDER_TRI)
		{
			el.nebrs[0].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]]);
			el.nebrs[0].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]]) / 3.0;
			el.nebrs[0].L = el.V = 0.0;
		}
		else if (el.type == elem::EType::BORDER_QUAD)
		{
			el.nebrs[0].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[3]]);
			el.nebrs[0].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[3]]) / 4.0;
			el.nebrs[0].L = el.V = 0.0;
		}
		else if (el.type == elem::EType::FRAC_QUAD)
		{
			el.nebrs[0].S = el.nebrs[1].S = point::square(pts[el.verts[0]], pts[el.verts[1]], pts[el.verts[2]], pts[el.verts[3]]);
			el.nebrs[0].cent = el.nebrs[1].cent = (pts[el.verts[0]] + pts[el.verts[1]] + pts[el.verts[2]] + pts[el.verts[3]]) / 4.0;
			el.nebrs[0].L = el.nebrs[1].L = el.V = 0.0;
		}
	}
};
int Mesh::check_neighbors() const
{
	int sum;
	for (const auto& el : cells)
		for (int i = 0; i < el.nebrs_num; i++)
		{
			sum = 0;
			const auto& el_nebr = cells[el.nebrs[i].id];
			for (int j = 0; j < el_nebr.nebrs_num; j++)
				if (el_nebr.nebrs[j].id == el.id)
					if (el_nebr.nebrs[j].cent == el.nebrs[i].cent)
						break;
					else
					{
						return i;
						exit(-1);
					}
		}
}
size_t Mesh::getCellsSize() const
{
	return cells.size();
}