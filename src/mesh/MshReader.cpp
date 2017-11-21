#include "src/mesh/MshReader.hpp"

using namespace mshreader;
using namespace std;

MshReader::MshReader() 
{
}
MshReader::~MshReader()
{
}
const grid::Mesh* MshReader::read(const string filename, const double x_dim)
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
	int counter = 0;
	vector<int> cells;
	while (buf != NODES_END)
	{
		msh >> x; msh >> y;	msh >> z;
		msh >> buf;
		cells.clear();
		cells.resize(stoi(buf, &sz));
		for (auto& cell_idx : cells)
		{
			msh >> buf;
			cell_idx = stoi(buf, &sz);
		}
		mesh->pts.push_back(point::Point(counter++, stod(x, &sz) / x_dim, stod(y, &sz) / x_dim, stod(z, &sz) / x_dim, cells));
		msh >> buf;
	}
	assert(mesh->pts.size() == mesh->pts_size);

	// Trash
	do
	msh >> buf;
	while (buf != ELEMS_BEGIN);
	// Elements
	msh >> buf; msh >> buf;

	counter = 0;
	while (buf != ELEMS_END)
	{
		auto readElem = [&, this](const elem::EType type)
		{
			for (int i = 0; i < elem::num_of_verts(type); i++)
				msh >> vertInds[i];
			for (int i = 0; i < elem::num_of_nebrs(type); i++)
				msh >> nebrInds[i];

			mesh->cells.push_back(elem::Element(counter++, type, vertInds, nebrInds));
		};

		msh >> buf;
		if (buf == to_string(elem::BORDER_TRI))			readElem(elem::BORDER_TRI);
		else if (buf == to_string(elem::BORDER_QUAD))	readElem(elem::BORDER_QUAD);
		else if (buf == to_string(elem::FRAC_QUAD))		readElem(elem::FRAC_QUAD);
		else if (buf == to_string(elem::PRISM))			readElem(elem::PRISM);
		else if (buf == to_string(elem::HEX))			readElem(elem::HEX);
		else if (buf == to_string(elem::BORDER_HEX))	readElem(elem::BORDER_HEX);
		else											getline(msh, buf);
		msh >> buf;
	}

	msh.close();

	return mesh;
}