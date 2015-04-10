// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_transformation.h"

#include "util.h"

#include <algorithm>
#include <limits>
#include <queue>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

const int kDirections[26][3] = {
    {0, 0, 1}, {0, 0, -1},
    {0, 1, 0}, {0, -1, 0},
    {1, 0, 0}, {-1, 0, 0},

    {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1},
    {1, 0, 1}, {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
    {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},

    {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
    {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};

int encode(int *index, int *dimensions) {
  return (index[2] * dimensions[1] + index[1]) * dimensions[0] + index[0];
}

int decode(int code, int *dimensions, int *index) {
  index[0] = code % dimensions[0];
  code /= dimensions[0];
  index[1] = code % dimensions[1];
  index[2] = code / dimensions[1];
}

bool outside(int *index, int *dimensions) {
  for (int c = 0; c < 3; c++) {
    if (index[c] < 0 || index[c] >= dimensions[c]) {
      return true;
    }
  }

  return false;
}

}

vtkStructuredPoints *WatershedTransformation::transform(
    vtkStructuredPoints *scalar_field,
    bool six_connectivity,
    int neighbor_limit,
    double neighbor_threshold) {
  int dimensions[3];
  double spacing[3];

  scalar_field->GetDimensions(dimensions);
  scalar_field->GetSpacing(spacing);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int *dist_2_lower = new int[nx * ny * nz];
  std::fill(dist_2_lower, dist_2_lower + nx * ny * nz,
            std::numeric_limits<int>::max());

  int connectivity = six_connectivity ? 6 : 26;

  std::queue<int> bfs_queue;

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = encode(curr_index, dimensions);
        double curr_scalar =
            scalar_field->GetPointData()->GetScalars()->GetTuple1(curr_code);

        for (int d = 0; d < connectivity; d++) {
          int next_x = x + kDirections[d][0];
          int next_y = y + kDirections[d][1];
          int next_z = z + kDirections[d][2];

          int next_index[] = {next_x, next_y, next_z};
          if (outside(next_index, dimensions)) {
            continue;
          }

          int next_code = encode(next_index, dimensions);
          double next_scalar =
              scalar_field->GetPointData()->GetScalars()->GetTuple1(next_code);

          if (next_scalar < curr_scalar) {
            dist_2_lower[curr_code] = 1;
            break;
          }
        }

        if (dist_2_lower[curr_code] == 1) {
          bfs_queue.push(curr_code);
        }
      }
    }
  }

  while (!bfs_queue.empty()) {
    int curr_code = bfs_queue.front();
    bfs_queue.pop();

    int curr_index[3];
    decode(curr_code, dimensions, curr_index);

    double curr_scalar =
        scalar_field->GetPointData()->GetScalars()->GetTuple1(curr_code);

    for (int d = 0; d < connectivity; d++) {
      int next_x = curr_index[0] + kDirections[d][0];
      int next_y = curr_index[1] + kDirections[d][1];
      int next_z = curr_index[2] + kDirections[d][2];

      int next_index[] = {next_x, next_y, next_z};
      if (outside(next_index, dimensions)) {
        continue;
      }

      int next_code = encode(next_index, dimensions);
      if (dist_2_lower[next_code] < std::numeric_limits<int>::max()) {
        continue;
      }

      double next_scalar =
          scalar_field->GetPointData()->GetScalars()->GetTuple1(next_code);

      if (next_scalar == curr_scalar) {
        dist_2_lower[next_code] = dist_2_lower[curr_code] + 1;
        queue.push(next_code);
      }
    }
  }

  for (int p = 0; p < nx * ny * nz; p++) {
    if (dist_2_lower[p] == std::numeric_limits<int>::max()) {
      int curr_index[3];
      decode(p, dimensions, curr_index);

      bool has_lower = false;
      for (int d = 0; d < connectivity; d++) {
        int next_x = curr_index[0] + kDirections[d][0];
        int next_y = curr_index[1] + kDirections[d][1];
        
      }
    }
  }

  delete [] dist_2_lower;
}

