#include "Square.h"

namespace JetscapeCornelius {

Square::Square() : ambiguous(false), number_cuts(0), number_lines(0) {}

Square::~Square() = default;

void Square::construct_lines(double value) {
  // Check the corner points to see if there are lines
  // Each comparison gives 0 or 1 → pack into a mask
  unsigned mask = (points[0][0] >= value) | ((points[0][1] >= value) << 1) |
                  ((points[1][0] >= value) << 2) |
                  ((points[1][1] >= value) << 3);
  // If all corners are above or below this value, there are no lines in this
  // square
  // If all 0 (0000) or all 1s (1111 = 15), no cuts
  if (mask == 0 || mask == 15) {
    number_lines = 0;
    return;
  }
  // Find the cut points and the points which are always outside of the
  // surface. Also find_outside() arranges cuts so that first two cuts form
  // a line as defined in the algorithm (in case there are 4 cuts)
  ends_of_edge(value);
  if (number_cuts > 0)
    find_outside(value);
  number_lines = 0;
  // Unroll for the most common case: number_cuts == 2
  if (number_cuts == 2) {
    // First endpoint
    points_temp[0][x1] = cuts[0][0];
    points_temp[0][x2] = cuts[0][1];
    points_temp[0][const_i[0]] = const_value[0];
    points_temp[0][const_i[1]] = const_value[1];
    // Second endpoint
    points_temp[1][x1] = cuts[1][0];
    points_temp[1][x2] = cuts[1][1];
    points_temp[1][const_i[0]] = const_value[0];
    points_temp[1][const_i[1]] = const_value[1];
    // Outside point
    out_temp[x1] = out[0][0];
    out_temp[x2] = out[0][1];
    out_temp[const_i[0]] = const_value[0];
    out_temp[const_i[1]] = const_value[1];
    lines[0].init_line(points_temp, out_temp, const_i);
    number_lines = 1;
    return;
  }
  // General case (number_cuts == 4)
  for (int i = 0; i < number_cuts; i += 2) {
    // First endpoint
    points_temp[0][x1] = cuts[i][0];
    points_temp[0][x2] = cuts[i][1];
    points_temp[0][const_i[0]] = const_value[0];
    points_temp[0][const_i[1]] = const_value[1];

    // Second endpoint
    points_temp[1][x1] = cuts[i + 1][0];
    points_temp[1][x2] = cuts[i + 1][1];
    points_temp[1][const_i[0]] = const_value[0];
    points_temp[1][const_i[1]] = const_value[1];

    // Outside point (paired with i/2)
    out_temp[x1] = out[i / 2][0];
    out_temp[x2] = out[i / 2][1];
    out_temp[const_i[0]] = const_value[0];
    out_temp[const_i[1]] = const_value[1];

    // Create the line
    lines[number_lines++].init_line(points_temp, out_temp, const_i);
  }
}

void Square::ends_of_edge(double value) {
  // Precompute values
  const double top_left = points[0][0] - value;
  const double top_right = points[0][1] - value;
  const double bottom_left = points[1][0] - value;
  const double bottom_right = points[1][1] - value;

  // Fast add_cut lambda
  auto add_cut = [this](double x, double y) {
    cuts[number_cuts][0] = x;
    cuts[number_cuts][1] = y;
    ++number_cuts;
  };

  // Edge 1
  if (top_left * bottom_left < 0) {
    add_cut(top_left / (points[0][0] - points[1][0]) * dx[x1], 0.0);
  } else if (points[0][0] == value && points[1][0] < value) {
    add_cut(ALMOST_ZERO * dx[x1], 0.0);
  } else if (points[1][0] == value && points[0][0] < value) {
    add_cut(ALMOST_ONE * dx[x1], 0.0);
  }

  // Edge 2
  if (top_left * top_right < 0) {
    add_cut(0.0, top_left / (points[0][0] - points[0][1]) * dx[x2]);
  } else if (points[0][0] == value && points[0][1] < value) {
    add_cut(0.0, ALMOST_ZERO * dx[x2]);
  } else if (points[0][1] == value && points[0][0] < value) {
    add_cut(0.0, ALMOST_ONE * dx[x2]);
  }

  // Edge 3
  if (bottom_left * bottom_right < 0) {
    add_cut(dx[x1], bottom_left / (points[1][0] - points[1][1]) * dx[x2]);
  } else if (points[1][0] == value && points[1][1] < value) {
    add_cut(dx[x1], ALMOST_ZERO * dx[x2]);
  } else if (points[1][1] == value && points[1][0] < value) {
    add_cut(dx[x1], ALMOST_ONE * dx[x2]);
  }

  // Edge 4
  if (top_right * bottom_right < 0) {
    add_cut(top_right / (points[0][1] - points[1][1]) * dx[x1], dx[x2]);
  } else if (points[0][1] == value && points[1][1] < value) {
    add_cut(ALMOST_ZERO * dx[x1], dx[x2]);
  } else if (points[1][1] == value && points[0][1] < value) {
    add_cut(ALMOST_ONE * dx[x1], dx[x2]);
  }

  if (number_cuts != 0 && number_cuts != 2 && number_cuts != 4) {
    std::cerr << "Error in ends_of_edge: number_cuts " << number_cuts
              << std::endl;
    exit(1);
  }
}

void Square::find_outside(double value) {
  if (number_cuts == 4) {
    // If there are 4 cuts, the surface is ambiguous
    ambiguous = true;

    // Compute the value in the middle of the square
    double value_middle =
        0.25 * (points[0][0] + points[0][1] + points[1][0] + points[1][1]);
    // The default is that cuts are connected as \\ here.
    // If both value_middle and (0,0) are above or below the criterion
    // the cuts should be like // and we have to switch order in cuts
    // Determine if cuts need to be swapped
    if ((points[0][0] < value && value_middle < value) ||
        (points[0][0] > value && value_middle > value)) {
      std::swap(cuts[1], cuts[2]);
    }

    // The center is below, so the middle point is always outside the surface
    if (value_middle < value) {
      out[0][0] = 0.5 * dx[x1];
      out[0][1] = 0.5 * dx[x2];
      out[1][0] = 0.5 * dx[x1];
      out[1][1] = 0.5 * dx[x2];
    } else {  // The center is above
      // Cuts are \\ here so bottom left and top right corners are outside
      if (points[0][0] < value) {
        out[0][0] = 0;
        out[1][0] = dx[x1];
        out[0][1] = 0;
        out[1][1] = dx[x2];
        // Cuts are // here so bottom right and top left corners are outside
      } else {
        out[0][0] = dx[x1];
        out[0][1] = 0;
        out[1][0] = 0;
        out[1][1] = dx[x2];
      }
    }
  } else {
    // This is the normal case (not ambiguous)
    out[0][0] = out[0][1] = out[1][0] = out[1][1] = 0.0;
    int number_out = 0;

    if (points[0][0] < value) {
      out[0][0] += 0 * dx[x1];
      out[0][1] += 0 * dx[x2];
      number_out++;
    }
    if (points[0][1] < value) {
      out[0][0] += 0 * dx[x1];
      out[0][1] += 1 * dx[x2];
      number_out++;
    }
    if (points[1][0] < value) {
      out[0][0] += 1 * dx[x1];
      out[0][1] += 0 * dx[x2];
      number_out++;
    }
    if (points[1][1] < value) {
      out[0][0] += 1 * dx[x1];
      out[0][1] += 1 * dx[x2];
      number_out++;
    }

    if (number_out > 0) {
      for (int i = 0; i < SQUARE_DIM; i++) {
        for (int j = 0; j < SQUARE_DIM; j++) {
          out[i][j] /= number_out;
        }
      }
    }
  }
}

}  // namespace JetscapeCornelius
