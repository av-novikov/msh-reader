#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkHexahedron.h>

#include "src/util/VTKSnapshotter.hpp"

#include "src/models/oil/Oil.hpp"

using namespace std;
using namespace snapshotter;

template<class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter(const Model* _model) : model(_model), mesh(_model->getMesh())
{
	R_dim = model->R_dim;
	pattern = prefix + "Mesh_%{STEP}.vtu";

	// Write cell types
	types = new int[mesh->cells.size()];
	for (int i = 0; i < mesh->inner_size; i++)
	{
		const auto& el = mesh->cells[i];
		if (el.type == elem::HEX)
			types[i] = VTK_HEXAHEDRON;
		else if (el.type == elem::PRISM)
			types[i] = VTK_WEDGE;
	}
}
template<class modelType>
VTKSnapshotter<modelType>::~VTKSnapshotter()
{
	delete types;
}
template<class modelType>
string VTKSnapshotter<modelType>::replace(string filename, string from, string to)
{
	size_t start_pos = 0;
	while ((start_pos = filename.find(from, start_pos)) != string::npos)
	{
		filename.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return filename;
}
template<class modelType>
std::string VTKSnapshotter<modelType>::getFileName(const int i)
{
	string filename = pattern;
	return replace(filename, "%{STEP}", to_string(i));
}

template<class modelType>
void VTKSnapshotter<modelType>::dump(const int i)
{
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto cells = vtkSmartPointer<vtkCellArray>::New();
	auto type = vtkSmartPointer<vtkIntArray>::New();
	type->SetName("type");

	points->Allocate(mesh->pts.size());
	cells->Allocate(mesh->cells.size());

	for (const auto& pt : mesh->pts)
		points->InsertNextPoint(pt.x * R_dim, pt.y * R_dim, pt.z * R_dim);

	for (int i = 0; i < mesh->inner_size; i++)
	{
		const auto& el = mesh->cells[i];
		if (el.type == elem::HEX)
		{
			auto vtkCell = vtkSmartPointer<vtkHexahedron>::New();
			for (int j = 0; j < el.verts_num; j++)
				vtkCell->GetPointIds()->SetId(j, el.verts[j]);
			cells->InsertNextCell(vtkCell);
		}
		else if (el.type == elem::PRISM)
		{
			auto vtkCell = vtkSmartPointer<vtkWedge>::New();
			for (int j = 0; j < el.verts_num; j++)
				vtkCell->GetPointIds()->SetId(j, el.verts[j]);
			cells->InsertNextCell(vtkCell);
		}

		type->InsertNextValue(el.type);
	}

	grid->SetPoints(points);
	grid->SetCells(types, cells);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(type);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<oil::Oil>;