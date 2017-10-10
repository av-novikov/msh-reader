#ifndef VTKSNAPSHOTTER_HPP_
#define VTKSNAPSHOTTER_HPP_

#include "Mesh.hpp"

namespace snapshotter
{
	class VTKSnapshotter
	{
	public:
		typedef grid::Mesh Mesh;
	protected:
		Mesh* mesh;
		const std::string prefix = "snaps/";
		std::string pattern;
		std::string replace(std::string filename, std::string from, std::string to);
		std::string getFileName(const int i);
	public:
		VTKSnapshotter();
		~VTKSnapshotter();

		void setGrid(Mesh* _mesh);

		void dump(const int i);
	};
};

#endif /* VTKSNAPSHOTTER_HPP_ */
