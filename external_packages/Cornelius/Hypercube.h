#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include "Cube.h"
#include "GeneralGeometryElement.h"
#include "Polyhedron.h"

namespace JetscapeCornelius {

/**
 * @class Hypercube
 * @brief A class representing a 4-dimensional hypercube for geometric
 * operations.
 *
 * The Hypercube class provides methods to initialize a hypercube, split into
 * cubes, check for ambiguities, construct polyhedra, and access geometric
 * elements.
 *
 * 13.10.2011 Hannu Holopainen
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 */
class Hypercube : public GeneralGeometryElement {
 private:
  static constexpr int DIM = 4;     ///< Dimension of the space.
  static constexpr int STEPS = 2;   ///< Number of steps for discretization.
  static constexpr int NCUBES = 8;  ///< Number of cubes in the hypercube.
  static constexpr int MAX_POLYHEDRONS =
      10;  ///< Maximum number of polyhedrons.

  std::array<std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>,
             STEPS>
      hypercube;                      ///< 4D array representing the hypercube.
  std::vector<Polyhedron> polyhedra;  ///< Vector to store the polyhedra.
  std::array<Cube, NCUBES>
      cubes;  ///< Array to store the cubes in the hypercube.
  std::array<const Polygon*, NCUBES * 10>
    polygon_refs;  ///< Temporary references to polygons collected from cubes.

  int number_polyhedra;        ///< Number of polyhedra in the hypercube.
  bool ambiguous;              ///< Indicates if the hypercube is ambiguous.
  std::array<double, DIM> dx;  ///< Delta values for discretization.
  std::array<char, NCUBES * 10> not_used;  ///< Array to track unused polygons.

  // Temporary array for storing a cube
  std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>
      cube;  ///< 3D array representing a cube

 public:
  /**
   * @brief Default constructor for the Hypercube class.
   */
  Hypercube();

  /**
   * @brief Destructor for the Hypercube class.
   */
  ~Hypercube();

  /**
   * @brief Initializes the hypercube with given parameters.
   * @param hc 4D array representing the hypercube.
   * @param new_dx Delta values for discretization.
   */
  inline void init_hypercube(
      std::array<
          std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>,
          STEPS>& hc,
      std::array<double, DIM>& new_dx) {
    hypercube = hc;
    dx = new_dx;
    number_polyhedra = 0;
    ambiguous = false;
    polyhedra.clear();
    polyhedra.reserve(MAX_POLYHEDRONS);
  }

  /**
   * @brief Splits the hypercube into cubes.
   * @param value The value used to decide if the point is inside the surface.
   * @return The number of points below the value.
   */
  int split_to_cubes(double value);

  /**
   * @brief Checks if the hypercube is ambiguous based on a given value.
   * @param number_points_below_value The number of points below the value.
   */
  inline void check_ambiguity(int number_points_below_value) {
    ambiguous = std::any_of(cubes.begin(), cubes.end(),
                            [](Cube& cube) { return cube.is_ambiguous(); });

    if (!ambiguous) {
      int number_lines = 0;
      for (auto& cube : cubes) {
        number_lines += cube.get_number_lines();
      }

      if (number_points_below_value > 8) {
        number_points_below_value = 16 - number_points_below_value;
      }
      if (number_lines == 24 && number_points_below_value == 2) {
        ambiguous = true;
      }
    }
  }

  /**
   * @brief Constructs polyhedra within the hypercube based on a given value.
   * @param value The value used to construct polyhedra.
   */
  void construct_polyhedra(double value);

  /**
   * @brief Gets the number of polyhedra in the hypercube.
   * @return The number of polyhedra.
   */
  inline int get_number_polyhedra() { return number_polyhedra; }

  /**
   * @brief Gets the polyhedra in the hypercube.
   * @return A reference to the array of polyhedra.
   */
  inline std::vector<Polyhedron>& get_polyhedra() { return polyhedra; }

  /**
   * @brief Checks if the hypercube is ambiguous.
   * @return True if the hypercube is ambiguous, false otherwise.
   */
  inline bool is_ambiguous() { return ambiguous; }
};

}  // namespace JetscapeCornelius

#endif  // HYPERCUBE_H
