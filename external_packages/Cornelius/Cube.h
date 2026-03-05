#ifndef CUBE_H
#define CUBE_H

#include <array>
#include <cstring>
#include <iostream>
#include <vector>

#include "GeneralGeometryElement.h"
#include "Line.h"
#include "Polygon.h"
#include "Square.h"

namespace JetscapeCornelius {

/**
 * @class Cube
 * @brief Represents a cube in a geometric space, derived from
 * GeneralGeometryElement.
 *
 * The Cube class encapsulates the properties and operations related to a cube
 * in a geometric space. It provides methods to initialize the cube, construct
 * polygons within the cube, split the cube into squares, and check for
 * ambiguity.
 *
 * 13.10.2011 Hannu Holopainen
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class Cube : public GeneralGeometryElement {
 private:
  static constexpr int DIM = 4;           ///< Dimension of the space.
  static constexpr int CUBE_DIM = 4;      ///< Dimension of the cube.
  static constexpr int NSQUARES = 6;      ///< Number of squares in the cube.
  static constexpr int STEPS = 2;         ///< Number of steps.
  static constexpr int MAX_POLYGONS = 8;  ///< Maximum number of polygons.

  std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>
      cube;                       ///< 3D array representing the cube.
  std::vector<Polygon> polygons;  ///< Vector to store the polygons in the cube.
  std::array<Square, NSQUARES> squares;  ///< Array of squares in the cube.
  std::array<Line, NSQUARES * 2> lines;  ///< Array of lines in the squares.
  std::array<const Line*, NSQUARES * 2>
      line_refs;  ///< Array to hold the line refs.

  int number_lines;            ///< Number of lines in the cube.
  int number_polygons;         ///< Number of polygons in the cube.
  bool ambiguous;              ///< Indicates if the cube is ambiguous.
  int const_i;                 ///< Index for the constant dimension.
  double const_value;          ///< Value for the constant dimension.
  int x1, x2, x3;              ///< Indices for dimensions.
  std::array<double, DIM> dx;  ///< Delta values for discretization.
  std::array<char, NSQUARES * 2> not_used;  ///< Array to track unused lines.

  // Temporary array for storing a square
  std::array<std::array<double, STEPS>, STEPS>
      square;  ///< 2D array representing a square

 public:
  /**
   * @brief Default constructor for the Cube class.
   */
  Cube();

  /**
   * @brief Destructor for the Cube class.
   */
  ~Cube();

  /**
   * @brief Initializes the cube with given parameters.
   * @param cu 3D array representing the cube.
   * @param new_const_i Index for the constant dimension.
   * @param new_const_value Value for the constant dimension.
   * @param new_dx Delta values for discretization.
   */
  inline void init_cube(
      std::array<std::array<std::array<double, STEPS>, STEPS>, STEPS>& cu,
      int new_const_i, double new_const_value,
      std::array<double, DIM>& new_dx) {
    cube = cu;
    const_i = new_const_i;
    const_value = new_const_value;
    dx = new_dx;
    // Use a lookup table for x1/x2/x3 assignment
    constexpr int x_lookup[4][3] = {
        {1, 2, 3},  // new_const_i == 0
        {0, 2, 3},  // new_const_i == 1
        {0, 1, 3},  // new_const_i == 2
        {0, 1, 2}   // new_const_i == 3
    };
    x1 = x_lookup[new_const_i][0];
    x2 = x_lookup[new_const_i][1];
    x3 = x_lookup[new_const_i][2];
    number_lines = number_polygons = 0;
    ambiguous = false;
    polygons.clear();
    polygons.reserve(MAX_POLYGONS);
    polygons.emplace_back();  // Default to construct 1 Polygon
  }

  /**
   * @brief Constructs polygons within the cube based on a given value.
   * @param value The value used to construct polygons.
   */
  void construct_polygons(double value);

  /**
   * @brief Splits the cube into squares.
   */
  void split_to_squares();

  /**
   * @brief Checks if the cube is ambiguous based on the number of lines.
   * @param number_lines The number of lines in the cube.
   */
  inline void check_ambiguity(int number_lines) {
    // Check if any squares may have ambiguous elements
    for (int i = 0; i < NSQUARES; i++) {
      if (squares[i].is_ambiguous()) {
        ambiguous = true;
        return;
      }
    }
    // If the surface is not ambiguous already, it is still possible to
    // have a ambiguous case if we have exactly 6 lines, i.e. the surface
    // elements are at the opposite corners
    if (number_lines == 6) {
      ambiguous = true;
    }
  }

  /**
   * @brief Checks if the cube is ambiguous.
   * @return True if the cube is ambiguous, false otherwise.
   */
  inline bool is_ambiguous() { return ambiguous; }

  /**
   * @brief Gets the number of polygons in the cube.
   * @return The number of polygons.
   */
  inline int get_number_polygons() { return number_polygons; }

  /**
   * @brief Gets the number of lines in the cube.
   * @return The number of lines.
   */
  inline int get_number_lines() { return number_lines; }

  /**
   * @brief Gets the polygons in the cube.
   * @return A reference to the array of polygons.
   */
  inline std::vector<Polygon>& get_polygons() { return polygons; }
};

}  // namespace JetscapeCornelius

#endif  // CUBE_H