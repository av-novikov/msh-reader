#ifndef MSHREADER_HPP_
#define MSHREADER_HPP_

#include <string>
#include <fstream>
#include <cstdlib>

#include "src/mesh/Mesh.hpp"

namespace mshreader
{
	class MshReader
	{
	private:
		std::string buf;
		std::string x, y, z;
		std::string::size_type sz;
		grid::Mesh* mesh;
		int vertInds[elem::MAX_ELEM_POINT_SIZE];
		int nebrInds[elem::MAX_ELEM_NEBR_SIZE];

		const std::string NODES_BEGIN = "$Nodes";
		const std::string NODES_END = "$EndNodes";
		const std::string ELEMS_BEGIN = "$Elements";
		const std::string ELEMS_END = "$EndElements";
	protected:
	public:
		MshReader();
		~MshReader();

		const grid::Mesh* read(const std::string filename, const double x_dim);
	};
};

#endif /* MSHREADER_HPP_ */