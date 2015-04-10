// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_transformation.h"

#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

namespace {

const char *kScalarFieldFile = "data/convective_cell_ftle.vtk";

const char *kTransformOutputFile = "transform.vtk";

const int kConnectivity = 6;
const int kNeighborLimit = 10;
const double kNeighborThreshold = 0.5;  // In percentage.

}

void transform_test() {
  printf("transform_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  vtkSmartPointer<vtkStructuredPoints> label_field =
      vtkSmartPointer<vtkStructuredPoints>(
          WatershedTransformation::transform(
              field, kConnectivity, kNeighborLimit, kNeighborThreshold));

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kTransformOutputFile);
  writer->SetInputData(label_field);
  writer->Write();

  printf("} transform_test\n");
}

int main() {
  transform_test();

  return 0;
}

