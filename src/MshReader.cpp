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
	msh >> buf;

	while (buf != ELEMS_END)
	{
		auto readElem = [&, this](const elem::EType type, bool isInner)
		{
			msh >> buf;	msh >> buf;	msh >> buf;
			for (int i = 0; i < elem::num_of_verts(type); i++)
			{
				msh >> vertInds[i];
				vertInds[i]--;
			}

			if (isInner)
				mesh->elems.push_back(elem::Element(type, vertInds));
			else
				border_elems.push_back(elem::Element(type, vertInds));
		};

		msh >> buf;
		if (buf == TRI_TYPE)		readElem(elem::EType::BORDER_TRI, false);
		else if (buf == QUAD_TYPE)	readElem(elem::EType::BORDER_QUAD, false);
		else if (buf == HEX_TYPE)	readElem(elem::EType::HEX, true);
		else if (buf == PRISM_TYPE)	readElem(elem::EType::PRISM, true);
		else						getline(msh, buf);

		msh >> buf;
	}

	msh.close();

	mesh->inner_size = mesh->elems.size();
	mesh->border_size = border_elems.size();
	mesh->elems.insert(end(mesh->elems), begin(border_elems), end(border_elems));

	for (int i = 0; i < mesh->elems.size(); i++)
		mesh->elems[i].num = i;

	return mesh;
}