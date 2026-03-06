#ifndef CORNELIUS_H
#define CORNELIUS_H

/**
 * This code is based on cornelius++ version 1.3:
 * Copyright 2012
 * Pasi Huovinen and Hannu Holopainen
 *
 * This subroutine is aimed to be used as a part of the fluid dynamical models
 * of the heavy-ion physics community. Permission to use it for any purpose
 * except for any commercial purpose is granted, provided that any publication
 * cites the paper describing the algorithm:
 *   P. Huovinen and H. Petersen, arXiv:1206.3371 [nucl-th]
 *   Eur.Phys.J.A 48 (2012) 171
 *
 * Permission to distribute this subroutine is granted, provided that no fee is
 * charged, and that this copyright and permission notice appear in all the
 * copies. Permission to modify this subroutine is granted provided that the
 * modified subroutine is made publicly available latest when any results
 * obtained using the modified subroutine are published, the modified subroutine
 * is distributed under terms similar to this notice, and the modified code
 * carries both the original copyright notice and notices stating that you
 * modified it, and a relevant date/year.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of FITNESS FOR A PARTICULAR
 * PURPOSE.
 *
 * Last update 03.08.2012 Hannu Holopainen
 *
 * Last update 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa: Modernized
 * the code to C++17 features, more memory efficient and faster. Added unit
 * tests and also tests to ensure that the new code is compatible with the
 * old one.
 *
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <ostream>
#include <vector>

#include "Cube.h"
#include "GeneralGeometryElement.h"
#include "Hypercube.h"
#include "Square.h"

namespace JetscapeCornelius {

/**
 *
 * A class for finding a constant value surface from a 2-4 dimensional
 * cube. The values at corners and length of the sides are needed as
 * an input.
 *
 * Algorithm by Pasi Huovinen. This code is based on the original FORTRAN
 * code by Pasi Huovinen.
 *
 * 23.04.2012 Hannu Holopainen
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class Cornelius : public GeneralGeometryElement {
 private:
  static constexpr int STEPS = 2; /**< Number of steps for the discretization */
  static constexpr int DIM = 4;   /**< Dimension of the space (default is 4D) */
  static constexpr int MAX_ELEMENTS = 10; /**< Maximum number of elements */

  int number_elements; /**< Number of surface elements found */
  std::array<std::array<double, DIM>, MAX_ELEMENTS>
      normals; /**< Array to store the normals of the surface elements */
  std::array<std::array<double, DIM>, MAX_ELEMENTS>
      centroids; /**< Array to store the centroids of the surface elements */
  int cube_dimension; /**< Dimension of the cube (2, 3, or 4) */
  bool initialized;   /**< Flag to indicate if Cornelius has been initialized */
  bool print_initialized; /**< Flag to indicate if printing is initialized */
  double value;           /**< Threshold value for surface detection */
  std::array<double, DIM> dx; /**< Array of step sizes in each dimension */
  std::ofstream
      output_file; /**< Output file stream for printing surface elements */

  Square cube_2d;    /**< 2D cube for surface detection */
  Cube cube_3d;      /**< 3D cube for surface detection */
  Hypercube cube_4d; /**< 4D cube for surface detection */

  /**
   * @brief Processes and finds the surface elements for a 3D cube.
   *
   * @param cu Values at the corners of the cube as a 3d table so that value
   *                      [0][0][0] is at (0,0,0) and [1][1][1] is at
   * (dx1,dx2,dx3).
   * @param position Absolute position at the point [0][0][0] in form
   * (0,x1,x2,x3).
   * @param do_print Boolean flag indicating if the surface elements should be
   * printed.
   */
  void surface_3d(
      std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu,
      std::array<double, DIM>& position, bool do_print);

 public:
  /**
   * @brief Default constructor for the Cornelius class.
   */
  Cornelius();

  /**
   * @brief Destructor for the Cornelius class.
   */
  ~Cornelius();

  /**
   * @brief Initializes the Cornelius object.
   *
   * @param dimension The dimension of the cube (2, 3, or 4).
   * @param new_value The value for surface.
   * @param new_dx Length of the sides of the cube. Must contain as many
   * elements as the dimension of the problem (dx1,dx2,...).
   */
  inline void init_cornelius(int dimension, double new_value,
                             const std::array<double, DIM>& new_dx) {
    cube_dimension = dimension;
    value = new_value;
    std::fill_n(dx.begin(), DIM - cube_dimension, 1.0);
    std::copy(new_dx.begin(), new_dx.begin() + cube_dimension,
              dx.begin() + (DIM - cube_dimension));
    initialized = true;
  }

  /**
   * @brief Initializes the output file for printing surface elements.
   *
   * @param filename The name of the output file.
   */
  void init_print_cornelius(std::string filename);

  /**
   * @brief Finds surface elements in a 2D cube.
   *
   * @param cu Values at the corners of the cube as a 2d table so that value
   *                  [0][0] is at (0,0,0) and [1][1] is at (dx1,dx2).
   */
  void find_surface_2d(std::array<std::array<double, STEPS>, STEPS>& cu);

  /**
   * @brief Finds surface elements in a 3D cube.
   *
   * @param cu Values at the corners of the cube as a 3d table so that value
   *                  [0][0][0] is at (0,0,0) and [1][1][1] is at (dx1,dx2,dx3).
   */
  void find_surface_3d(
      std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu);

  /**
   * @brief Finds the surface elements in 3-dimensional case and prints the
   * actual triangles which are found by the algorithm.
   *
   * @param cu Values at the corners of the cube as a 3d table so that value
   *                  [0][0][0] is at (0,0,0) and [1][1][1] is at (dx1,dx2,dx3).
   * @param position Absolute position at the point [0][0][0] in form
   * (0,x1,x2,x3).
   */
  void find_surface_3d_print(
      std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu,
      std::array<double, DIM>& position);

  /**
   * @brief Finds surface elements in a 4D hypercube.
   *
   * @param cu Values at the corners of the cube as a 4d table so that value
   *                  [0][0][0][0] is at (0,0,0,0) and [1][1][1][1] is at
   *                  (dx1,dx2,dx3,dx4).
   */
  void find_surface_4d(
      std::array<
          std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>,
          STEPS>& cu);

  /**
   * @brief Gets the number of surface elements found.
   *
   * @return The number of surface elements.
   */
  inline int get_number_elements() { return number_elements; }

  /**
   * @brief Normal vectors as a 2d table with the following number of indices
   * [number of elements][dimension of the problem]. This gives \sigma_\mu
   * without factors(sqrt(-g)) from the metric.
   *
   * @return A vector of vectors with dimensions [number of
   * elements][dimension of the problem] containing the normal vectors of the
   * surface elements.
   */
  std::vector<std::vector<double>> get_normals();

  /**
   * @brief Centroid vectors as a 2d table with the following number of indices
   * [number of elements][dimension of the problem].
   *
   * @return A vector of vectors representing the centroids.
   */
  std::vector<std::vector<double>> get_centroids();

  /**
   * @brief Gets a specific centroid element.
   *
   * @param index_surface_element The index of the surface element. Valid
   *               values are [0,number of elements in this cube].
   * @param element_centroid The index of the centroid element. Valid values are
   *               [0,dimension of the problem].
   * @return The value of the specified centroid element.
   */
  double get_centroid_element(int index_surface_element, int element_centroid);

  /**
   * @brief Gets a specific normal element.
   *
   * @param index_surface_element The index of the surface element. Valid
   *               values are [0,number of elements in this cube].
   * @param element_normal The index of the normal element. Valid values are
   *               [0,dimension of the problem].
   * @return The value of the specified normal element.
   */
  double get_normal_element(int index_surface_element, int element_normal);
};

}  // namespace JetscapeCornelius

#endif  // CORNELIUS_H