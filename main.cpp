#include <string>
#include <iostream>
#include <memory>

#include "src/mesh/MshReader.hpp"
#include "src/VTKSnapshotter.hpp"

using namespace std;

int main()
{
	mshreader::MshReader reader;
	std::shared_ptr<grid::Mesh> mesh = make_shared<grid::Mesh>(*reader.read("attempt2.nebr"));
	mesh->process_geometry();

	snapshotter::VTKSnapshotter snapshot;
	snapshot.setGrid(mesh.get());
	snapshot.dump(0);

	return 0;
}