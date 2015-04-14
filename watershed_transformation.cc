// Author: Mingcheng Chen (linyufly@gmail.com)

#include "watershed_transformation.h"

#include "util.h"

#include <algorithm>
#include <limits>
#include <map>
#include <queue>
#include <vector>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

struct PersistencePair {
  int label_1_, label_2_;
  double persistence_;

  PersistencePair() {
  }

  PersistencePair(int label_1, int label_2, double persistence)
      : label_1_(label_1), label_2_(label_2), persistence_(persistence) {}

  bool operator < (const PersistencePair &other) const {
    return persistence_ < other.persistence_;
  }
};

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

int find_root(int node, int *father) {
  if (father[node] == node) {
    return node;
  }

  return father[node] = find_root(father[node], father);
}

// Merge n_1 to n_2. The basin of n_2 should have smaller scalar.
void merge(int n_1, int n_2, int *father) {
  int r_1 = find_root(n_1, father);
  int r_2 = find_root(n_2, father);

  if (r_1 != r_2) {
    father[r_1] = r_2;
  }
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

  double max_scalar = -std::numeric_limits<double>::max();
  double min_scalar = std::numeric_limits<double>::max();

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = encode(curr_index, dimensions);
        double curr_scalar =
            scalar_field->GetPointData()->GetScalars()->GetTuple1(curr_code);

        if (curr_scalar < min_scalar) {
          min_scalar = curr_scalar;
        }

        if (curr_scalar > max_scalar) {
          max_scalar = curr_scalar;
        }

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

  std::vector<double> label_scalars;

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

      label_scalars.push_back(start_scalar);

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

  double **min_neighbor = create_matrix<double>(nx * ny * nz, neighbor_limit + 1);
  for (int p = 0; p < nx * ny * nz; p++) {
    min_neighbor[p][0] = scalar_field->GetPointData()->GetScalars()->GetTuple1(p);
  }

  for (int dist = 1; dist <= neighbor_limit; dist++) {
    for (int p = 0; p < nx * ny * nz; p++) {
      int curr_index[3];
      decode(p, dimensions, curr_index);

      min_neighbor[p][dist] = min_neighbor[p][dist - 1];
      for (int d = 0; d < connectivity; d++) {
        int next_x = curr_index[0] + kDirections[d][0];
        int next_y = curr_index[1] + kDirections[d][1];
        int next_z = curr_index[2] + kDirections[d][2];

        int next_index[] = {next_x, next_y, next_z};
        if (outside(next_index, dimensions)) {
          continue;
        }

        int next_code = encode(next_index, dimensions);
        if (labels[next_code] != labels[p]) {
          continue;
        }

        min_neighbor[p][dist] = std::min(min_neighbor[p][dist],
                                         min_neighbor[next_code][dist - 1]);
      }
    }
  }

  int *father = new int[num_labels];
  for (int l = 0; l < num_labels; l++) {
    father[l] = l;
  }

  double absolute_threshold = (max_scalar - min_scalar) * neighbor_threshold;

  for (int p = 0; p < nx * ny * nz; p++) {
    int curr_index[3];
    decode(p, dimensions, curr_index);

    for (int d = 0; d < connectivity; d++) {
      int next_x = curr_index[0] + kDirections[d][0];
      int next_y = curr_index[1] + kDirections[d][1];
      int next_z = curr_index[2] + kDirections[d][2];

      int next_index[] = {next_x, next_y, next_z};
      if (outside(next_index, dimensions)) {
        continue;
      }

      int next_code = encode(next_index, dimensions);

      int root_1 = find_root(labels[p], father);
      int root_2 = find_root(labels[next_code], father);
      if (root_1 == root_2) {
        continue;
      }

      double scalar_1 =
          scalar_field->GetPointData()->GetScalars()->GetTuple1(p);
      double scalar_2 =
          scalar_field->GetPointData()->GetScalars()->GetTuple1(next_code);

      if (scalar_1 < label_scalars[root_1]
          || scalar_2 < label_scalars[root_2]) {
        printf("scalar_1 = %lf, minima_1 = %lf\n",
               scalar_1, label_scalars[root_1]);
        printf("scalar_2 = %lf, minima_2 = %lf\n",
               scalar_2, label_scalars[root_2]);
        report_error("Minima is larger than ridge.\n");
      }

      if (scalar_1 - min_neighbor[p][neighbor_limit] < absolute_threshold
          && scalar_2 - min_neighbor[next_code][neighbor_limit]
          < absolute_threshold) {
        if (label_scalars[root_1] > label_scalars[root_2]) {
          merge(root_1, root_2, father);
        } else {
          merge(root_2, root_1, father);
        }
      }
    }
  }

  for (int p = 0; p < nx * ny * nz; p++) {
    label_array->SetTuple1(
        p, static_cast<double>(label_scalars[find_root(labels[p], father)]));
  }

  label_field->GetPointData()->SetScalars(label_array);

  delete [] dist_2_lower;
  delete [] point_order;
  delete [] labels;
  delete [] father;

  delete_matrix(min_neighbor);

  return label_field;
}

vtkStructuredPoints *WatershedTransformation::transform(
    vtkStructuredPoints *scalar_field,
    bool six_connectivity,
    double persistence_threshold) {
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

  double max_scalar = -std::numeric_limits<double>::max();
  double min_scalar = std::numeric_limits<double>::max();

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = encode(curr_index, dimensions);
        double curr_scalar =
            scalar_field->GetPointData()->GetScalars()->GetTuple1(curr_code);

        if (curr_scalar < min_scalar) {
          min_scalar = curr_scalar;
        }

        if (curr_scalar > max_scalar) {
          max_scalar = curr_scalar;
        }

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

  std::vector<double> label_scalars;

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

      label_scalars.push_back(start_scalar);

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

  std::vector<std::map<int, double> > persistence_links(num_labels);

  for (int p = 0; p < nx * ny * nz; p++) {
    int curr_index[3];
    decode(p, dimensions, curr_index);

    for (int d = 0; d < connectivity; d++) {
      int next_x = curr_index[0] + kDirections[d][0];
      int next_y = curr_index[1] + kDirections[d][1];
      int next_z = curr_index[2] + kDirections[d][2];

      int next_index[] = {next_x, next_y, next_z};
      if (outside(next_index, dimensions)) {
        continue;
      }

      int next_code = encode(next_index, dimensions);

      int label_1 = labels[p];
      int label_2 = labels[next_code];

      if (label_1 == label_2) {
        continue;
      }

      double scalar_1 =
          scalar_field->GetPointData()->GetScalars()->GetTuple1(p);
      double scalar_2 =
          scalar_field->GetPointData()->GetScalars()->GetTuple1(next_code);

      if (scalar_1 < label_scalars[label_1]
          || scalar_2 < label_scalars[label_2]) {
        printf("scalar_1 = %lf, minima_1 = %lf\n",
               scalar_1, label_scalars[label_1]);
        printf("scalar_2 = %lf, minima_2 = %lf\n",
               scalar_2, label_scalars[label_2]);
        report_error("Minima is larger than ridge.\n");
      }

      double persistence = std::min(scalar_1 - label_scalars[label_1],
                                    scalar_2 - label_scalars[label_2]);

      bool need_update = false;
      if (persistence_links[label_1].find(label_2) == persistence_links[label_1].end()) {
        need_update = true;
      } else if (persistence_links[label_1][label_2] > persistence) {
        need_update = true;
      }

      if (need_update) {
        persistence_links[label_1][label_2] = persistence_links[label_2][label_1] = persistence;
      }
    }
  }
 
  std::vector<PersistencePair> pairs;
  for (int l = 0; l < num_labels; l++) {
    for (std::map<int, double>::iterator itr = persistence_links[l].begin();
         itr != persistence_links[l].end(); ++itr) {
      if (itr->first == l) {
        report_error("Persistence pair contains two same labels.\n");
      }
      if (itr->first > l) {
        continue;
      }
      pairs.push_back(PersistencePair(l, itr->first, itr->second));
    }
  }

  std::sort(pairs.begin(), pairs.end());

  int *father = new int[num_labels];
  for (int l = 0; l < num_labels; l++) {
    father[l] = l;
  }

  double absolute_threshold = (max_scalar - min_scalar) * persistence_threshold;

  for (int p = 0; p < pairs.size(); p++) {
    int l_1 = pairs[p].label_1_;
    int l_2 = pairs[p].label_2_;
    double persistence = pairs[p].persistence_;
    if (persistence > absolute_threshold) {
      break;
    }

    int r_1 = find_root(l_1, father);
    int r_2 = find_root(l_2, father);

    if (r_1 == r_2) {
      continue;
    }

    if (label_scalars[r_1] > label_scalars[r_2]) {
      
    }
  }


  for (int p = 0; p < nx * ny * nz; p++) {
    label_array->SetTuple1(
        p, static_cast<double>(label_scalars[find_root(labels[p], father)]));
  }

  label_field->GetPointData()->SetScalars(label_array);

  delete [] dist_2_lower;
  delete [] point_order;
  delete [] labels;
  delete [] father;

  delete_matrix(min_neighbor);

  return label_field;

}

