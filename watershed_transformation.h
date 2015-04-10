// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef WATERSHED_TRANSFORMATION_H_
#define WATERSHED_TRANSFORMATION_H_

class vtkStructuredPoints;

class WatershedTransformation {
 public:
  static vtkStructuredPoints *transform(
      vtkStructuredPoints *scalar_field,
      bool six_connectivity,
      int neighbor_limit,
      double neighbor_threshold);
};

#endif  // WATERSHED_TRANSFORMATION_H_

