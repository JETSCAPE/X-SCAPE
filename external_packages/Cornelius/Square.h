#ifndef SQUARE_H
#define SQUARE_H

#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>

#include "GeneralGeometryElement.h"
#include "Line.h"

namespace JetscapeCornelius {

/**
 * @class Square
 * @brief A class representing a square element.
 *
 * This class handles squares. Finds the edges of the surface and also
 * the points which are always outside the surface so that we can determine
 * correct direction for normal vector.
 *
 * 13.10.2011 Hannu Holopainen
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class Square : public GeneralGeometryElement {
 private:
  static constexpr int DIM = 4;         ///< Dimension of the space.
  static constexpr int SQUARE_DIM = 2;  ///< Dimension of the square.
  static constexpr int MAX_POINTS = 4;  ///< Maximum number of points.
  static constexpr int MAX_LINES = 2;   ///< Maximum number of lines.

  static constexpr double ALMOST_ONE = 1.0 - 1e-9;  ///< Almost one value.
  static constexpr double ALMOST_ZERO = 1e-9;       ///< Almost zero value.
  std::array<std::array<double, SQUARE_DIM>, SQUARE_DIM>
      points;  ///< Points of the square.
  std::array<std::array<double, SQUARE_DIM>, MAX_POINTS>
      cuts;  ///< Points of the cuts.
  std::array<std::array<double, SQUARE_DIM>, MAX_POINTS>
      out;                                    ///< Points outside the square.
  std::array<int, DIM - SQUARE_DIM> const_i;  ///< Indices for constraints.
  std::array<double, DIM - SQUARE_DIM>
      const_value;                    ///< Values for constraints.
  int x1, x2;                         ///< Indices for dimensions of the square.
  std::array<double, DIM> dx;         ///< Delta values for lines.
  int number_cuts;                    ///< Number of cuts.
  int number_lines;                   ///< Number of lines.
  std::array<Line, MAX_LINES> lines;  ///< Lines in the square.
  bool ambiguous;                     ///< Indicates if the square is ambiguous.

  std::array<std::array<double, DIM>, SQUARE_DIM>
      points_temp;                   ///< Temporary points.
  std::array<double, DIM> out_temp;  ///< Temporary outside point.

 public:
  /**
   * @brief Default constructor for the Square class.
   */
  Square();

  /**
   * @brief Destructor for the Square class.
   */
  ~Square();

  /**
   * @brief Initializes the square with given parameters.
   * @param sq Points defining the square.
   * @param c_i Constant indices for the square.
   * @param c_v Values for the constant indices.
   * @param dex Delta values for lines.
   */
  inline void init_square(
      std::array<std::array<double, SQUARE_DIM>, SQUARE_DIM>& sq,
      std::array<int, DIM - SQUARE_DIM>& c_i,
      std::array<double, DIM - SQUARE_DIM>& c_v, std::array<double, DIM>& dex) {
    points = sq;
    const_i = c_i;
    const_value = c_v;
    dx = dex;
    x1 = x2 = -1;
    for (int i = 0; i < DIM; i++) {
      if (i != const_i[0] && i != const_i[1]) {
        (x1 < 0 ? x1 : x2) = i;
      }
    }
    number_cuts = number_lines = 0;
    ambiguous = false;
  }

  /**
   * @brief Constructs lines within the square based on a given value.
   * @param value The value used to construct lines.
   */
  void construct_lines(double value);

  /**
   * @brief Determines the ends of edges based on a given value.
   *
   * This finds the points from the edges where the surface crosses
   * the given square.
   *
   * The edges are gone through in the following order
   *      4
   *    -----
   *  1 |   | 3
   *    -----
   *      2
   * Since in the case where the ends are on each edge the ends are
   * assumed to be connected like \\
   *
   * @param value The value used to find ends of edges.
   */
  void ends_of_edge(double value);

  /**
   * @brief Finds the points outside the square based on a given value.
   * @param value The value used to find outside points.
   */
  void find_outside(double value);

  /**
   * @brief Checks if the square is ambiguous.
   * @return True if the square is ambiguous, false otherwise.
   */
  inline bool is_ambiguous() { return ambiguous; }

  /**
   * @brief Gets the number of lines in the square.
   * @return The number of lines.
   */
  inline int get_number_lines() { return number_lines; }

  /**
   * @brief Gets the lines in the square.
   * @return A reference to the array of lines.
   */
  inline std::array<Line, MAX_LINES>& get_lines() { return lines; }
};

}  // namespace JetscapeCornelius

#endif  // SQUARE_H