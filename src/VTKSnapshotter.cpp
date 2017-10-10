#include "VTKSnapshotter.hpp"
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

using namespace std;
using namespace snapshotter;

VTKSnapshotter::VTKSnapshotter()
{
	pattern = prefix + "Mesh_%{STEP}.vtu";
}
VTKSnapshotter::~VTKSnapshotter()
{
}
void VTKSnapshotter::setGrid(Mesh* _mesh)
{
	mesh = _mesh;
}
string VTKSnapshotter::replace(string filename, string from, string to)
{
	size_t start_pos = 0;
	while ((start_pos = filename.find(from, start_pos)) != string::npos)
	{
		filename.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return filename;
}
std::string VTKSnapshotter::getFileName(const int i)
{
	string filename = pattern;
	return replace(filename, "%{STEP}", to_string(i));
}

void VTKSnapshotter::dump(const int i)
{
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto cells = vtkSmartPointer<vtkCellArray>::New();
	auto type = vtkSmartPointer<vtkIntArray>::New();
	type->SetName("type");

	points->Allocate(mesh->pts.size());
	cells->Allocate(mesh->elems.size());

	for (const auto& pt : mesh->pts)
		points->InsertNextPoint(pt.x, pt.y, pt.z);

	for (int i = 0; i < mesh->inner_size; i++)
	{
		const auto& el = mesh->elems[i];
		if (el.type == elem::HEX)
		{
			auto vtkCell = vtkSmartPointer<vtkHexahedron>::New();
			for (int j = 0; j < el.verts_num; j++)
				vtkCell->GetPointIds()->SetId(j, el.verts[j]);

			cells->InsertNextCell(vtkCell);
			type->InsertNextValue(el.type);
		}
		/*else if (el.type == elem::PRISM)
		{
			auto vtkCell = vtkSmartPointer<vtkWedge>::New();
			for (int j = 0; j < el.verts_num; j++)
				vtkCell->GetPointIds()->SetId(j, el.verts[j]);

			cells->InsertNextCell(vtkCell);
		}

		type->InsertNextValue(el.type);*/
	}

	grid->SetPoints(points);
	grid->SetCells(VTK_HEXAHEDRON, cells);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(type);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}