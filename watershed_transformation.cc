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

struct PointComparator {
  vtkStructuredPoints *scalar_field_;
  int *dist_2_lower_;

  PointComparator(vtkStructuredPoints *scalar_field, int *dist_2_lower)
      : scalar_field_(scalar_field), dist_2_lower_(dist_2_lower) {
  }

  bool operator () (int p_1, int p_2) const {
    if (p_1 == p_2) {
      return false;
    }

    if (p_1 == -1) {
      return false;
    }

    if (p_2 == -1) {
      return true;
    }

    double s_1 = scalar_field_->GetPointData()->GetScalars()->GetTuple1(p_1);
    double s_2 = scalar_field_->GetPointData()->GetScalars()->GetTuple1(p_2);

    if (s_1 != s_2) {
      return s_1 < s_2;
    }

    return dist_2_lower_[p_1] < dist_2_lower_[p_2];
  }
};

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

  /// DEBUG ///
  printf("Dimensions: %d, %d, %d.\n", nx, ny, nz);

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
        bfs_queue.push(next_code);
      }
    }
  }

  // Check every local minima.
  for (int p = 0; p < nx * ny * nz; p++) {
    if (dist_2_lower[p] == std::numeric_limits<int>::max()) {
      int curr_index[3];
      decode(p, dimensions, curr_index);

      double curr_scalar =
          scalar_field->GetPointData()->GetScalars()->GetTuple1(p);

      bool has_lower = false;
      for (int d = 0; d < connectivity; d++) {
        int next_x = curr_index[0] + kDirections[d][0];
        int next_y = curr_index[1] + kDirections[d][1];
        int next_z = curr_index[2] + kDirections[d][2];

        int next_index[] = {next_x, next_y, next_z};
        if (outside(next_index, dimensions)) {
          continue;
        }

        int next_code = encode(next_index, dimensions);
        double next_scalar =
            scalar_field->GetPointData()->GetScalars()->GetTuple1(next_code);

        if (next_scalar < curr_scalar) {
          has_lower = true;
          break;
        }
      }

      if (has_lower) {
        report_error("False local minima found.\n");
      }
    }
  }

  vtkStructuredPoints *label_field = vtkStructuredPoints::New();
  label_field->CopyStructure(scalar_field);

  vtkSmartPointer<vtkDoubleArray> label_array =
      vtkSmartPointer<vtkDoubleArray>::New();
  label_array->SetNumberOfComponents(1);
  label_array->SetNumberOfTuples(nx * ny * nz);

  int *point_order = new int[nx * ny * nz];

  for (int p = 0; p < nx * ny * nz; p++) {
    point_order[p] = p;
  }

  PointComparator point_comparator(scalar_field, dist_2_lower);
  std::sort(point_order, point_order + nx * ny * nz, point_comparator);

  int *labels = new int[nx * ny * nz];
  std::fill(labels, labels + nx * ny * nz, -1);

  int num_labels = 0;

  for (int order = 0; order < nx * ny * nz; order++) {
    int start_code = point_order[order];
    if (labels[start_code] != -1) {
      continue;
    }

    double start_scalar =
        scalar_field->GetPointData()->GetScalars()->GetTuple1(start_code);

    if (dist_2_lower[start_code] == std::numeric_limits<int>::max()) {
      labels[start_code] = num_labels++;

      std::queue<int> bfs_queue;
      bfs_queue.push(start_code);
      while (!bfs_queue.empty()) {
        int curr_code = bfs_queue.front();
        bfs_queue.pop();

        int curr_index[3];
        decode(curr_code, dimensions, curr_index);

        for (int d = 0; d < connectivity; d++) {
          int next_x = curr_index[0] + kDirections[d][0];
          int next_y = curr_index[1] + kDirections[d][1];
          int next_z = curr_index[2] + kDirections[d][2];

          int next_index[] = {next_x, next_y, next_z};
          if (outside(next_index, dimensions)) {
            continue;
          }

          int next_code = encode(next_index, dimensions);

          if (labels[next_code] != -1) {
            if (labels[next_code] != labels[start_code]) {
              report_error("Conflict in local minima labelling.\n");
            }

            continue;
          }

          double next_scalar =
              scalar_field->GetPointData()->GetScalars()->GetTuple1(next_code);
          if (next_scalar != start_scalar) {
            continue;
          }

          labels[next_code] = labels[start_code];
          bfs_queue.push(next_code);
        }
      }
    } else {
      int candidate = -1;
      int curr_index[3];
      decode(start_code, dimensions, curr_index);
      for (int d = 0; d < connectivity; d++) {
        int next_x = curr_index[0] + kDirections[d][0];
        int next_y = curr_index[1] + kDirections[d][1];
        int next_z = curr_index[2] + kDirections[d][2];

        int next_index[] = {next_x, next_y, next_z};
        if (outside(next_index, dimensions)) {
          continue;
        }

        int next_code = encode(next_index, dimensions);
        if (point_comparator(next_code, candidate)) {
          candidate = next_code;
        }
      }

      if (scalar_field->GetPointData()->GetScalars()->GetTuple1(candidate)
          > start_scalar) {
        report_error("Conflict in non-minima labelling.\n");
      }

      if (labels[candidate] == -1) {
        report_error("Points are misordered.\n");
      }

      labels[start_code] = labels[candidate];
    }
  }

  /// DEBUG ///
  printf("num_labels = %d\n", num_labels);

  for (int p = 0; p < nx * ny * nz; p++) {
    label_array->SetTuple1(p, static_cast<double>(labels[p]));
  }

  label_field->GetPointData()->SetScalars(label_array);

  delete [] dist_2_lower;
  delete [] point_order;
  delete [] labels;

  return label_field;
}

