#ifndef EDGE_HPP_
#define EDGE_HPP_

#include "src/util/Point.hpp"

namespace elem
{
	template <int N>
	class Edge
	{
	public:
		std::array<const point::Point*,N> verts;

	};
}

#endif /* EDGE_HPP_ */