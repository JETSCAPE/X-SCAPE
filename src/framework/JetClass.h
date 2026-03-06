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

#ifndef JetClass_h
#define JetClass_h

#include <stdio.h>
#include <math.h>
#include "JetScapeParticles.h"
#include "JetScapeConstants.h"
#include "FourVector.h"
#include "fjcore.hh"

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

using std::ostream;

namespace Jetscape {

// class Parton;
class Vertex;
class FourVector;

/**************************************************************************************************/

//  JET CLASS

/*************************************************************************************************/

/**
 * @class Jet
 * @brief Placeholder class representing a jet object in the JETSCAPE framework.
 *
 * The Jet class is currently a dummy implementation, intended to be expanded
 * once the jet graph structure is fully implemented.
 */
class Jet {
  /**
   * @brief Default constructor for Jet.
   */
  Jet(){};

  /**
   * @brief Default destructor for Jet.
   */
  ~Jet(){};
};

/**************************************************************************************************/

//  VERTEX CLASS

/*************************************************************************************************/

/**
 * @class Vertex
 * @brief Represents a space-time vertex in the JETSCAPE framework.
 *
 * A Vertex stores its four-dimensional position using the FourVector class.
 * In future extensions, connections to parents, children, and siblings will be
 * handled through a graph structure.
 */
class Vertex {
 public:
  /**
   * @brief Default constructor. Initializes the vertex at the origin.
   */
  Vertex() { x_in_.Set(0, 0, 0, 0); }

  /**
   * @brief Constructor with explicit coordinates.
   * @param x Spatial x-coordinate.
   * @param y Spatial y-coordinate.
   * @param z Spatial z-coordinate.
   * @param t Time coordinate.
   */
  Vertex(double x, double y, double z, double t) { x_in_.Set(x, y, z, t); }

  /**
   * @brief Constructor from an existing FourVector.
   * @param x FourVector specifying the vertex location.
   */
  Vertex(FourVector &x) { set_location(x); }

  /**
   * @brief Destructor for Vertex.
   */
  virtual ~Vertex();
  /**
   * @brief Set the location of the vertex.
   * @param x FourVector specifying the new location.
   */
  void set_location(FourVector &x) { x_in_ = x; }

  /**
   * @brief Get the current location of the vertex.
   * @return Reference to the FourVector representing the location.
   */
  FourVector &x_in() { return (x_in_); }

  /**
   * @brief Stream output operator for Vertex.
   *
   * Prints the space-time coordinates of the vertex in the order:
   * `x y z t`.
   *
   * @param output Output stream.
   * @param vertex Vertex to print.
   * @return Reference to the modified output stream.
   */
  friend ostream &operator<<(ostream &output, Vertex &vertex) {
    output << vertex.x_in().x() << " " << vertex.x_in().y() << " "
           << vertex.x_in().z() << " " << vertex.x_in().t();
    return output;
  }

 protected:
  FourVector x_in_;  // location of the vertex
  // parents and siblings from Graph structure later ...
};

};  // namespace Jetscape

#endif /* JetClass_h */
