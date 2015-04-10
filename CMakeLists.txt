cmake_minimum_required(VERSION 2.8)

project(WatershedTransformation)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

add_executable(WatershedTransformationTest watershed_transformation_test.cc
                                           watershed_transformation.cc
                                           util.cc)
target_link_libraries(WatershedTransformationTest ${VTK_LIBRARIES})

