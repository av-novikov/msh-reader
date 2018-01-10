#ifndef ELEM_HPP_
#define ELEM_HPP_

#include <array>

namespace elem
{
	static const int MAX_ELEM_POINT_SIZE = 8;
	static const int MAX_ELEM_NEBR_SIZE = 6;
	static const int QUAD_VERT_SIZE = 4;
	static const int TRI_VERT_SIZE = 3;
	static const int MAX_INTERFACES_STENCIL = 6;
	static const int NEARESET_POINT_NUMBER = 3;

	enum EType { BORDER_TRI, BORDER_QUAD, FRAC_QUAD, PRISM, HEX, BORDER_HEX };
	inline const char num_of_verts(const EType type)
	{
		if (type == BORDER_TRI)
			return 3;
		else if (type == BORDER_QUAD || type == FRAC_QUAD)
			return 4;
		else if (type == PRISM)
			return 6;
		else if (type == HEX || type == BORDER_HEX)
			return 8;
	}
	inline const char num_of_nebrs(const EType type)
	{
		if (type == BORDER_TRI || type == BORDER_QUAD)
			return 1;
		else if (type == FRAC_QUAD)
			return 2;
		else if (type == PRISM)
			return 5;
		else if (type == HEX || type == BORDER_HEX)
			return 6;
	}

	struct Nebr
	{
		point::Id nebr;
		double S;
		double L;
		point::Point cent;
		point::Interaction* ireg [2];
		point::Point n, nu;
		int nearest_cells[NEARESET_POINT_NUMBER];
		double nearest_dists[NEARESET_POINT_NUMBER];
	};

	class Element
	{
	public:
		const int id;
		EType type;

		char verts_num;
		char nebrs_num;
		std::array<int, MAX_ELEM_POINT_SIZE> verts;
		std::array<Nebr, MAX_ELEM_NEBR_SIZE> nebrs;
		std::array<double, MAX_ELEM_NEBR_SIZE-2> T;
		point::Point cent;
		double V;

		Element() : id(-1) {};
		Element(const int _id, const EType _type, const int* _verts) : id(_id), type(_type), verts_num(num_of_verts(_type)), nebrs_num(num_of_nebrs(_type))
		{
			for (int i = 0; i < verts_num; i++)
				verts[i] = _verts[i];
		};
		Element(const int _id, const EType _type, const int* _verts, const int* _nebrs) : id(_id), type(_type), verts_num(num_of_verts(_type)), nebrs_num(num_of_nebrs(_type))
		{
			for (int i = 0; i < verts_num; i++)
				verts[i] = _verts[i];
			for (int i = 0; i < nebrs_num; i++)
				nebrs[i].nebr.cell = _nebrs[i];
		};
		~Element() {};
	};
};

#endif /* ELEM_HPP_ */
