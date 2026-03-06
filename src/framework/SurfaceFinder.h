/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// This is a general basic class for a hyper-surface finder

#ifndef SURFACEFINDER_H_
#define SURFACEFINDER_H_

#include <array>
#include <fstream>
#include <omp.h>
#include <vector>

#include "FluidEvolutionHistory.h"
#include "RealType.h"
#include "SurfaceCellInfo.h"

namespace Jetscape {

/**
 * @brief The `SurfaceFinder` class is responsible for identifying the
 * freeze-out hypersurface in a hydrodynamic evolution history.
 *
 * This class utilizes the Cornelius algorithm to find the hypersurface where
 * the temperature equals a specified cutoff value. It processes the bulk
 * evolution history, checks for intersections with the hypersurface, and
 * constructs a list of `SurfaceCellInfo` objects that represent the properties
 * of each surface cell.
 */
class SurfaceFinder {
 private:
  /// @brief The cutoff temperature for identifying the freeze-out hypersurface.
  Jetscape::real T_cut;

  /// @brief Reference to the bulk evolution history containing fluid
  /// properties.
  const EvolutionHistory &bulk_info;

  /// @brief Flag indicating whether the hydrodynamic medium is boost invariant.
  bool boost_invariant;

  /// @brief List of surface cells that form the identified hypersurface.
  std::vector<SurfaceCellInfo> surface_cell_list;

 public:
  /**
   * @brief Prepares a `SurfaceCellInfo` object with the given parameters.
   *
   * This function constructs a `SurfaceCellInfo` object by populating its
   * fields with the provided space-time coordinates, normal vector components,
   * and fluid cell properties. The resulting `SurfaceCellInfo` encapsulates all
   * necessary information about a surface element for further processing or
   * output.
   *
   * @param tau Proper time coordinate of the surface cell.
   * @param x X coordinate of the surface cell.
   * @param y Y coordinate of the surface cell.
   * @param eta Pseudorapidity coordinate of the surface cell.
   * @param da0 Normal vector component in the tau direction.
   * @param da1 Normal vector component in the x direction.
   * @param da2 Normal vector component in the y direction.
   * @param da3 Normal vector component in the eta direction.
   * @param fluid_cell The `FluidCellInfo` object containing fluid properties at
   * the surface cell location.
   * @return A fully populated `SurfaceCellInfo` object representing the surface
   * element.
   */
  SurfaceFinder(const Jetscape::real T_in, const EvolutionHistory &bulk_data);

  /**
   * @brief Destructor for the `SurfaceFinder` class.
   *
   * This destructor clears the `surface_cell_list` to free up memory resources.
   * It ensures that all dynamically allocated memory for surface cells is
   * released when a `SurfaceFinder` object goes out of scope or is explicitly
   * deleted.
   */
  ~SurfaceFinder();

  /**
   * @brief Finds and constructs the full hypersurface in a 3D space-time grid.
   *
   * This function iterates through a predefined space-time grid to identify
   * the freeze-out hypersurface using the Cornelius algorithm. It initializes
   * a 3D grid, iterates over time and spatial coordinates, checks for
   * intersections, and extracts surface elements to store them in
   * `surface_cell_list`.
   */
  void Find_full_hypersurface();

  /**
   * @brief gets the number of surface cells in the identified hypersurface.
   * @return The total number of surface cells in the `surface_cell_list`.
   */
  int get_number_of_surface_cells() const { return (surface_cell_list.size()); }

  /**
   * @brief Retrieves a `SurfaceCellInfo` object from the `surface_cell_list` by
   * index.
   * @return A `SurfaceCellInfo` object corresponding to the specified index.
   */
  SurfaceCellInfo get_surface_cell_with_idx(int idx) const {
    return (surface_cell_list[idx]);
  }

  /**
   * @brief Retrieves the entire list of surface cells as a vector.
   * @return A vector containing all `SurfaceCellInfo` objects in the
   * `surface_cell_list`.
   */
  std::vector<SurfaceCellInfo> get_surface_cells() const {
    return (surface_cell_list);
  }

  /**
   * @brief Retrieves the entire list of surface cells as a vector.
   * @return A vector containing all `SurfaceCellInfo` objects in the
   * `surface_cell_list`.
   */
  std::vector<SurfaceCellInfo> get_surface_cells_vector() const {
    return (surface_cell_list);
  }

  /**
   * @brief Checks if the temperature values in the cube intersect the cutoff
   * temperature.
   *
   * @param tau Central value of tau.
   * @param x Central value of x.
   * @param y Central value of y.
   * @param dt Time step size.
   * @param dx X step size.
   * @param dy Y step size.*
   * @param cube 3D array to store temperature values of the grid cell.
   * @note explicit multidimensional std::array types for 2x2x2 and
   * 2x2x2x2 cubes 3D cube: [2][2][2]
   * @return True if the temperature values intersect the cutoff temperature,
   * false otherwise.
   */
  bool check_intersect_3D(
      Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real dt,
      Jetscape::real dx, Jetscape::real dy,
      std::array<std::array<std::array<double, 2>, 2>, 2> &cube);

  /**
   * @brief Finds and constructs the full hypersurface in a 3D space-time grid.
   *
   * This function iterates through a predefined space-time grid to identify
   * the freeze-out hypersurface using the Cornelius algorithm. It initializes
   * a 3D grid, iterates over time and spatial coordinates, checks for
   * intersections, and extracts surface elements to store them in
   * `surface_cell_list`.
   *
   * @note This function dynamically allocates memory for a 3D cube and ensures
   * proper cleanup.
   */
  void Find_full_hypersurface_3D();

  /**
   * @brief Checks for an intersection between the hypersurface and a 4D
   * space-time grid cell.
   *
   * This function determines whether the freeze-out hypersurface intersects
   * a given space-time grid cell by evaluating the temperature values at its
   * corners. If an intersection is detected, it updates the provided `cube`
   * with temperature values.
   *
   * @param tau   Proper time coordinate at the center of the cell.
   * @param x     X coordinate at the center of the cell.
   * @param y     Y coordinate at the center of the cell.
   * @param eta   Pseudorapidity coordinate at the center of the cell.
   * @param dt    Proper time step size.
   * @param dx    X step size.
   * @param dy    Y step size.
   * @param deta  Eta step size.
   * @param cube  4D array to store temperature values at grid points.
   *
   * @note 4D cube: [2][2][2][2]
   *
   * @return True if an intersection occurs, false otherwise.
   */
  bool check_intersect_4D(
      Jetscape::real tau, Jetscape::real x, Jetscape::real y,
      Jetscape::real eta, Jetscape::real dt, Jetscape::real dx,
      Jetscape::real dy, Jetscape::real deta,
      std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube);

  /**
   * @brief Finds the full hypersurface in 4D space-time by identifying
   * isothermal surfaces.
   *
   * This function iterates over a 4D grid in space-time (τ, x, y, η) and
   * extracts hypersurface elements where the temperature crosses the critical
   * value
   * (`T_cut`). It employs the Cornelius isosurface finder to locate and store
   * surface elements.
   *
   * @details
   * - Initializes the grid and retrieves limits from `bulk_info`.
   * - Allocates a 4D array (`cube`) to store temperature values at neighboring
   * grid points.
   * - Loops over time (`τ`), space-time rapidity (`η`), and transverse plane
   * (`x, y`).
   * - Calls `check_intersect_4D()` to determine intersections.
   * - If an intersection is detected, `Cornelius` finds the isothermal
   * hypersurface.
   * - Surface elements are extracted, including centroids and normal vectors.
   * - Fluid properties are retrieved and stored as `SurfaceCell` elements.
   * - Cleans up dynamically allocated memory at the end.
   *
   * @note This function is computationally expensive due to its iteration over
   * a 4D space-time grid. Parallelization with OpenMP or Kokkos can improve
   * performance.
   *
   * @see check_intersect_4D(), Cornelius::find_surface_4d(),
   * PrepareASurfaceCell()
   */
  void Find_full_hypersurface_4D();

  /**
   * @brief Prepares a SurfaceCellInfo object containing surface and fluid
   * properties.
   *
   * This function constructs a `SurfaceCellInfo` structure by populating it
   * with spatial coordinates, normal vectors, hydrodynamic properties, and flow
   * velocities from a given `FluidCellInfo` object.
   *
   * @param tau Proper time coordinate of the surface cell.
   * @param x X-coordinate of the surface cell.
   * @param y Y-coordinate of the surface cell.
   * @param eta Space-time rapidity of the surface cell.
   * @param da0 Normal vector component in the τ-direction.
   * @param da1 Normal vector component in the x-direction.
   * @param da2 Normal vector component in the y-direction.
   * @param da3 Normal vector component in the η-direction.
   * @param fluid_cell The hydrodynamic fluid cell containing thermodynamic
   * properties.
   *
   * @return A `SurfaceCellInfo` object populated with spatial, normal vector,
   *         thermodynamic, and flow velocity data.
   *
   *
   * @note The four-velocity transformation ensures compatibility with the
   * space-time rapidity coordinate.
   * @see SurfaceCellInfo, FluidCellInfo
   */
  SurfaceCellInfo PrepareASurfaceCell(Jetscape::real tau, Jetscape::real x,
                                      Jetscape::real y, Jetscape::real eta,
                                      Jetscape::real da0, Jetscape::real da1,
                                      Jetscape::real da2, Jetscape::real da3,
                                      const FluidCellInfo fluid_cell);

  /**
   * @brief Writes the surface cell information to a specified file.
   *
   * This function appends the string representation of each `SurfaceCellInfo`
   * object in the provided vector to a file with the given filename. Each
   * surface cell's data is written in a formatted manner for easy readability.
   *
   * @param surface_cells A vector containing `SurfaceCellInfo` objects to be
   * written to the file.
   * @param filename The name of the file where the surface cell information
   * will be appended.
   */
  void WriteSurfaceToFile(const std::vector<SurfaceCellInfo> &surface_cell_list,
                          std::string filename);
};

}  // namespace Jetscape

#endif  // SURFACEFINDER_H_
