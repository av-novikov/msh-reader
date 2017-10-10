#include "MshReader.hpp"

using namespace mshreader;
using namespace std;

MshReader::MshReader() 
{
}
MshReader::~MshReader()
{
}
const grid::Mesh* MshReader::read(const string filename)
{
	mesh = new grid::Mesh;
	ifstream msh;
	msh.open(filename.c_str(), ifstream::in);

	// Trash
	do
		msh >> buf;
	while (buf != NODES_BEGIN);
	// Number of nodes
	msh >> mesh->pts_size;
	// Nodes
	msh >> buf;
	while (buf != NODES_END)
	{
		msh >> x; msh >> y;	msh >> z;
		mesh->pts.push_back(point::Point(stod(x, &sz), stod(y, &sz), stod(z, &sz)));
		msh >> buf;
	}
	assert(mesh->pts.size() == mesh->pts_size);

	// Trash
	do
	msh >> buf;
	while (buf != ELEMS_BEGIN);
	// Elements
	msh >> buf; msh >> buf;

	while (buf != ELEMS_END)
	{
		auto readElem = [&, this](const elem::EType type)
		{
			for (int i = 0; i < elem::num_of_verts(type); i++)
				msh >> vertInds[i];
			for (int i = 0; i < elem::num_of_nebrs(type); i++)
				msh >> nebrInds[i];

			mesh->elems.push_back(elem::Element(type, vertInds, nebrInds));
		};

		msh >> buf;
		if (buf == to_string(elem::BORDER_TRI))			readElem(elem::EType::BORDER_TRI);
		else if (buf == to_string(elem::BORDER_QUAD))	readElem(elem::EType::BORDER_QUAD);
		else if (buf == to_string(elem::FRAC_QUAD))		readElem(elem::EType::FRAC_QUAD);
		else if (buf == to_string(elem::PRISM))			readElem(elem::EType::PRISM);
		else if (buf == to_string(elem::HEX))			readElem(elem::EType::HEX);
		else											getline(msh, buf);
		msh >> buf;
	}

	msh.close();

	for (int i = 0; i < mesh->elems.size(); i++)
		mesh->elems[i].num = i;

	return mesh;
}