#ifndef MSHREADER_HPP_
#define MSHREADER_HPP_

#include <string>
#include <fstream>
#include <cstdlib>

#include "Mesh.hpp"

namespace mshreader
{
	class MshReader
	{
	private:
		std::string buf;
		std::string x, y, z;
		std::string::size_type sz;
		grid::Mesh* mesh;
		std::vector<elem::Element> border_elems;
		int vertInds[elem::MAX_ELEM_NEBR_SIZE];

		const std::string NODES_BEGIN = "$Nodes";
		const std::string NODES_END = "$EndNodes";
		const std::string ELEMS_BEGIN = "$Elements";
		const std::string ELEMS_END = "$EndElements";
		const std::string TRI_TYPE = "2";
		const std::string QUAD_TYPE = "3";
		const std::string HEX_TYPE = "5";
		const std::string PRISM_TYPE = "6";
	protected:
	public:
		MshReader();
		~MshReader();

		const grid::Mesh* read(const std::string filename);
	};
};

#endif /* MSHREADER_HPP_ */