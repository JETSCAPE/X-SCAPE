#include "Cornelius.h"

namespace JetscapeCornelius {

Cornelius::Cornelius()
    : number_elements(0),
      cube_dimension(0),
      initialized(false),
      print_initialized(false) {}

Cornelius::~Cornelius() {
  // Close the file if it is open
  if (output_file.is_open()) {
    output_file.close();
  }
}

void Cornelius::init_print_cornelius(std::string filename) {
  output_file.open(filename);
  print_initialized = true;
}

void Cornelius::find_surface_2d(
    std::array<std::array<double, STEPS>, STEPS>& cu) {
  if (!initialized || cube_dimension != 2) {
    std::cerr << "Cornelius not initialized for 2D case." << std::endl;
    exit(1);
  }
  std::array<int, 2> c_i = {0, 1};
  std::array<double, 2> c_v = {0, 0};
  cube_2d.init_square(cu, c_i, c_v, dx);
  cube_2d.construct_lines(value);
  number_elements = cube_2d.get_number_lines();
  for (int i = 0; i < number_elements; i++) {
    for (int j = 0; j < DIM; j++) {
      normals[i][j] = cube_2d.get_lines()[i].get_normal()[j];
      centroids[i][j] = cube_2d.get_lines()[i].get_centroid()[j];
    }
  }
}

void Cornelius::find_surface_3d(
    std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu) {
  std::array<double, DIM> position = {0};
  surface_3d(cu, position, false);
}

void Cornelius::find_surface_3d_print(
    std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu,
    std::array<double, DIM>& position) {
  surface_3d(cu, position, true);
}

void Cornelius::surface_3d(
    std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu,
    std::array<double, DIM>& position, bool do_print) {
  if (!initialized || cube_dimension != 3) {
    std::cerr << "Cornelius not initialized for 3D case." << std::endl;
    exit(1);
  }
  // Check if the cube actually contains surface elements.
  // If all or none of the elements are below the criterion, no surface
  // elements exist.
  int value_greater = 0;
  for (const auto& array2d : cu) {
    for (const auto& array1d : array2d) {
      for (double element : array1d) {
        value_greater += (element >= this->value);
      }
    }
  }
  if (value_greater == 0 || value_greater == 8) {
    number_elements = 0;
    return;
  }
  const int c_i = 0;
  const double c_v = 0.0;
  cube_3d.init_cube(cu, c_i, c_v, dx);
  cube_3d.construct_polygons(value);
  number_elements = cube_3d.get_number_polygons();
  auto& polys = cube_3d.get_polygons();
  for (int i = 0; i < number_elements; ++i) {
    auto& poly = polys[i];
    auto& normal = poly.get_normal();
    auto& centroid = poly.get_centroid();
    normals[i] = normal;
    centroids[i] = centroid;
    if (print_initialized && do_print) {
      poly.print(output_file, position);
    }
  }
}

void Cornelius::find_surface_4d(
    std::array<std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>,
               STEPS>& cu) {
  if (!initialized || cube_dimension != 4) {
    std::cerr << "Cornelius not initialized for 4D case." << std::endl;
    exit(1);
  }
  // Check if the cube actually contains surface elements.
  // If all or none of the elements are below the criterion, no surface
  // elements exist.
  int value_greater = 0;
  for (const auto& array3d : cu) {
    for (const auto& array2d : array3d) {
      for (const auto& array1d : array2d) {
        for (double element : array1d) {
          value_greater += (element >= this->value);
        }
      }
    }
  }
  if (value_greater == 0 || value_greater == 16) {
    number_elements = 0;
    return;
  }
  cube_4d.init_hypercube(cu, dx);
  cube_4d.construct_polyhedra(value);
  number_elements = cube_4d.get_number_polyhedra();
  auto& polys = cube_4d.get_polyhedra();
  for (int i = 0; i < number_elements; ++i) {
    auto& poly = polys[i];
    auto& normal = poly.get_normal();
    auto& centroid = poly.get_centroid();
    normals[i] = normal;
    centroids[i] = centroid;
  }
}

std::vector<std::vector<double>> Cornelius::get_normals() {
  std::vector<std::vector<double>> normals_vector(
      number_elements, std::vector<double>(cube_dimension));
  for (int i = 0; i < number_elements; i++) {
    std::copy(normals[i].begin() + (DIM - cube_dimension),
              normals[i].begin() + DIM, normals_vector[i].begin());
  }
  return normals_vector;
}

std::vector<std::vector<double>> Cornelius::get_centroids() {
  std::vector<std::vector<double>> centroids_vector(
      number_elements, std::vector<double>(cube_dimension));
  for (int i = 0; i < number_elements; i++) {
    std::copy(centroids[i].begin() + (DIM - cube_dimension),
              centroids[i].begin() + DIM, centroids_vector[i].begin());
  }
  return centroids_vector;
}

double Cornelius::get_centroid_element(int index_surface_element,
                                       int element_centroid) {
  if (index_surface_element >= number_elements ||
      element_centroid >= cube_dimension) {
    throw std::out_of_range(
        "Cornelius error: asking for an element which does not exist.");
  }
  return centroids[index_surface_element]
                  [element_centroid + (DIM - cube_dimension)];
}

double Cornelius::get_normal_element(int index_surface_element,
                                     int element_normal) {
  if (index_surface_element >= number_elements ||
      element_normal >= cube_dimension) {
    throw std::out_of_range(
        "Cornelius error: asking for an element which does not exist.");
  }
  return normals[index_surface_element]
                [element_normal + (DIM - cube_dimension)];
}

}  // namespace JetscapeCornelius