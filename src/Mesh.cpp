#include "Mesh.hpp"
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
	set_neighbors();
}
bool Mesh::are_adjanced(const elem::Element& el1, const elem::Element& el2)
{
	int sum = 0;

	for (int i = 0; i < el1.verts_num; i++)
		sum += (find(begin(el2.verts), end(el2.verts), el1.verts[i]) == end(el2.verts) ? 0 : 1);

	if (sum > 2) return true;
	else return false;
}
void Mesh::set_neighbors()
{
	/*for (int i = 0; i < inner_size; i++)
	{
		auto& el = elems[i];
		if (el.type == elem::EType::HEX)
			for (const auto& el_nebr : elems)
			{
				if (el_nebr.num != el.num)
				{
					bool b[6] = { false };
					auto find_in_verts4 = [&](int v1, int v2, int v3, int v4) -> bool
					{
						bool b1, b2, b3, b4;
						b1 = b2 = b3 = b4 = false;
						for (int j = 0; j < el_nebr.verts_num; j++)
						{
							if (!b1) if (v1 == el_nebr.verts[j]) { b1 = true; continue; }
							if (!b2) if (v2 == el_nebr.verts[j]) { b2 = true; continue; }
							if (!b3) if (v3 == el_nebr.verts[j]) { b3 = true; continue; }
							if (!b4) if (v4 == el_nebr.verts[j]) { b4 = true; continue; }
						}
						return b1 * b2 * b3 * b4;
					};
					if (!b[0]) if (find_in_verts4(el.verts[0], el.verts[1], el.verts[2], el.verts[3])) { el.nebrs[0].id = el_nebr.num; b[0] = true; };
					if (!b[1]) if (find_in_verts4(el.verts[4], el.verts[5], el.verts[6], el.verts[7])) { el.nebrs[1].id = el_nebr.num; b[1] = true; };
					if (!b[2]) if (find_in_verts4(el.verts[0], el.verts[1], el.verts[5], el.verts[4])) { el.nebrs[2].id = el_nebr.num; b[2] = true; };
					if (!b[3]) if (find_in_verts4(el.verts[1], el.verts[2], el.verts[6], el.verts[5])) { el.nebrs[3].id = el_nebr.num; b[3] = true; };
					if (!b[4]) if (find_in_verts4(el.verts[2], el.verts[3], el.verts[7], el.verts[6])) { el.nebrs[4].id = el_nebr.num; b[4] = true; };
					if (!b[5]) if (find_in_verts4(el.verts[3], el.verts[0], el.verts[4], el.verts[7])) { el.nebrs[5].id = el_nebr.num; b[5] = true; };
				}
			}
		else if (el.type == elem::EType::PRISM)
			for (const auto& el_nebr : elems)
			{
				if (el_nebr.num != el.num)
				{
					bool b[5] = { false };

					auto find_in_verts3 = [&](int v1, int v2, int v3) -> bool
					{
						bool b1, b2, b3;
						b1 = b2 = b3 = false;
						for (int j = 0; j < el_nebr.verts_num; j++)
						{
							if (!b1) if (v1 == el_nebr.verts[j]) { b1 = true; continue; }
							if (!b2) if (v2 == el_nebr.verts[j]) { b2 = true; continue; }
							if (!b3) if (v3 == el_nebr.verts[j]) { b3 = true; continue; }
						}
						return b1 * b2 * b3;
					};
					auto find_in_verts4 = [&](int v1, int v2, int v3, int v4) -> bool
					{
						bool b1, b2, b3, b4;
						b1 = b2 = b3 = b4 = false;
						for (int j = 0; j < el_nebr.verts_num; j++)
						{
							if (!b1) if (v1 == el_nebr.verts[j]) { b1 = true; continue; }
							if (!b2) if (v2 == el_nebr.verts[j]) { b2 = true; continue; }
							if (!b3) if (v3 == el_nebr.verts[j]) { b3 = true; continue; }
							if (!b4) if (v4 == el_nebr.verts[j]) { b4 = true; continue; }
						}
						return b1 * b2 * b3 * b4;
					};
					if (!b[0]) if (find_in_verts3(el.verts[0], el.verts[1], el.verts[2])) { el.nebrs[0].id = el_nebr.num; b[0] = true; };
					if (!b[1]) if (find_in_verts3(el.verts[3], el.verts[4], el.verts[5])) { el.nebrs[1].id = el_nebr.num; b[1] = true; };
					if (!b[2]) if (find_in_verts4(el.verts[0], el.verts[1], el.verts[4], el.verts[3])) { el.nebrs[2].id = el_nebr.num; b[2] = true; };
					if (!b[3]) if (find_in_verts4(el.verts[1], el.verts[2], el.verts[5], el.verts[4])) { el.nebrs[3].id = el_nebr.num; b[3] = true; };
					if (!b[4]) if (find_in_verts4(el.verts[2], el.verts[0], el.verts[3], el.verts[5])) { el.nebrs[4].id = el_nebr.num; b[4] = true; };
				}
			}

		for (int j = 0; i < el.nebrs_num; i++)
			assert(el.nebrs[j].id > 0 && el.nebrs[j].id < elems.size());

	}*/

	kdtree = kd_create(3);
	for (auto& el : elems)
		for(int i = 0; i < el.nebrs_num; i++)
			kd_insert3(kdtree, el.nebrs[i].cent.x, el.nebrs[i].cent.y, el.nebrs[i].cent.z, &el);

	double pos[3];
	struct kdres* res;
	const double range = 100 * EQUALITY_TOLERANCE;
	for (auto& el : elems)
	{
		for (int i = 0; i < el.nebrs_num; i++)
		{
			const point::Point& pt = el.nebrs[i].cent;
			res = kd_nearest_range3(kdtree, pt.x, pt.y, pt.z, range);
			while (!kd_res_end(res))
			{
				const auto& el_nebr = *(elem::Element*)kd_res_item(res, pos);
				if (el_nebr.num == el.num) continue;
				el.nebrs[0].id = el_nebr.num;
				break;
			}
			kd_res_free(res);
		}
	}

	int aaa = check_neighbors();
}
void Mesh::set_geom_props()
{
	for (auto& el : elems)
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
	}
};
int Mesh::check_neighbors() const
{
	int sum;
	for (const auto& el : elems)
		for (int i = 0; i < el.nebrs_num; i++)
		{
			sum = 0;
			const auto& el_nebr = elems[el.nebrs[i].id];
			for (int j = 0; j < el_nebr.nebrs_num; j++)
				if (el_nebr.nebrs[j].id == el.num)
					if (el_nebr.nebrs[j].cent == el.nebrs[i].cent)
						break;
					else
					{
						return i;
						exit(-1);
					}
		}
}