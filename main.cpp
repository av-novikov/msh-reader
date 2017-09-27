#include <string>
#include <iostream>
#include <memory>

#include "src/MshReader.hpp"

using namespace std;

int main()
{
	mshreader::MshReader reader;
	std::shared_ptr<grid::Mesh> mesh = make_shared<grid::Mesh>(*reader.read("attempt2.msh"));
	mesh->process_geometry();

	return 0;
}