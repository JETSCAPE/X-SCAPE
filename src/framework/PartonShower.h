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

// PartonShower with graph from GTL

#ifndef PARTONSHOWER_H
#define PARTONSHOWER_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
//#include "JetScapeParticles.hpp"
#include "JetScapeLogger.h"
//#include "helper.h"

#include <map>

using std::shared_ptr;

namespace Jetscape {

class Vertex;
class Parton;

/**
 * @class PartonShower
 * @brief Represents a parton shower as a directed graph.
 *
 * This class provides an interface for constructing and analyzing parton
 * showers in a graph-theoretic framework. Each **node** corresponds to a
 * space-time splitting vertex, and each **edge** corresponds to a parton
 * propagating between two vertices.
 *
 * - Internally inherits from GTL's `graph` class.
 * - Provides mappings from GTL nodes/edges to JETSCAPE `Vertex` and `Parton`
 * objects.
 * - Supports export to multiple graph formats (GML, Graphviz `.gv`, GraphML).
 * - Provides utilities for analyzing shower structure, such as retrieving
 *   final-state partons and vertex/parton counts.
 */
class PartonShower : public graph {
 public:
  /**
   * @brief Construct an empty parton shower.
   */
  PartonShower();

  /**
   * @brief Destructor. Clears final parton storage.
   */
  virtual ~PartonShower();

  /**
   * @brief Create a new vertex in the shower graph.
   * @param v Pointer to a `Vertex` object to associate with the new node.
   * @return The GTL node created in the graph.
   */
  node new_vertex(std::shared_ptr<Vertex> v);

  /**
   * @brief Create a new parton (edge) in the shower graph.
   * @param s Source vertex (parent).
   * @param t Target vertex (child).
   * @param p Pointer to a `Parton` object to associate with the edge.
   * @return The unique integer ID of the edge.
   */
  int new_parton(node s, node t, std::shared_ptr<Parton> p);

  /**
   * @brief Clone the parton shower graph.
   * @note not yet working for 2--2 !!! Fix it !!!
   * @return A unique pointer to a new `PartonShower` object.
   */
  virtual std::unique_ptr<PartonShower> Clone();

  /**
   * @brief Insert a new parton into the shower graph.
   * @param e Existing edge to split.
   * @param v Vertex to insert at the split point.
   * @param p Parton to insert on the new edge.
   */
  void InsertParton(edge e, std::shared_ptr<Vertex> v,
                    std::shared_ptr<Parton> p);

  /**
   * @brief Insert a new parton after an existing edge in the shower graph.
   * @param e Existing edge to split.
   * @param v Vertex to insert at the split point (will copy position of
   * target).
   * @param p Parton to insert on the new edge.
   */
  void InsertPartonAfter(edge e, std::shared_ptr<Vertex> v,
                         std::shared_ptr<Parton> p);

  /**
   * @brief Insert a new edge into the shower graph.
   * @param e Existing edge to split.
   * @param eIns Edge to insert at the split point.
   */
  void InsertEdge(edge e, edge eIns){};

  /**
   * @brief Retrieve the vertex object associated with a given graph node.
   * @param n Node in the graph.
   * @return Shared pointer to the corresponding `Vertex`.
   */
  std::shared_ptr<Vertex> GetVertex(node n) { return vMap[n]; }

  /**
   * @brief Retrieve the parton object associated with a given graph edge.
   * @param e Edge in the graph.
   * @return Shared pointer to the corresponding `Parton`.
   */
  std::shared_ptr<Parton> GetParton(edge e) { return pMap[e]; }

  /**
   * @brief Retrieve the parton associated with the nth edge.
   * @param n Index of the edge (0-based).
   * @return Shared pointer to the `Parton`.
   */
  std::shared_ptr<Parton> GetPartonAt(int n);

  /**
   * @brief Retrieve the vertex associated with the nth node.
   * @param n Index of the node (0-based).
   * @return Shared pointer to the `Vertex`.
   */
  std::shared_ptr<Vertex> GetVertexAt(int n);

  /**
   * @brief Retrieve the graph edge associated with a given parton.
   * @param p Shared pointer to a `Parton`.
   * @return GTL edge handle corresponding to the parton.
   */
  edge GetEdge(std::shared_ptr<Parton> p) { return eMap.at(p); }

  /**
   * @brief Retrieve the graph node associated with a given vertex.
   * @param v Shared pointer to a `Vertex`.
   * @return GTL node handle corresponding to the vertex.
   */
  node GetNode(std::shared_ptr<Vertex> v) { return nMap.at(v); }

  /**
   * @brief Retrieve the nth node in the graph.
   * @param n Index of the node (0-based).
   * @return GTL node handle.
   */
  node GetNodeAt(int n);

  /**
   * @brief Retrieve the nth edge in the graph.
   * @param n Index of the edge (0-based).
   * @return GTL edge handle.
   */
  edge GetEdgeAt(int n);

  /**
   * @brief Check if a node corresponds to an N-to-1 splitting vertex.
   * @param n Node to check.
   */
  bool IsNtoOne(node n);

  /**
   * @brief Check if a node corresponds to a 1-to-1 splitting vertex.
   * @param n Node to check.
   */
  bool IsOneToOne(node n);

  /**
   * @brief Check if a node corresponds to a 1-to-2 splitting vertex.
   * @param n Node to check.
   */
  bool IsOneToTwo(node n);

  /**
   * @brief Check if a node corresponds to an end vertex (no outgoing edges).
   * @param n Node to check.
   */
  bool IsEndNode(node n);

  bool IsHighSplitEdge(edge e);
  bool IsHighSplitNode(node n);

  double GetSplitDeltaR(node n);
  double GetSplitZ(node n);
  double GetSplitKt(node n);

  double GetNextNodeTime(node n);
  double GetNodeTime(node n);  /// {return GetVertex(n)->x_in().t();}

  edge GetHighSplitEdge(node n);
  edge GetLowSplitEdge(node n);

  /**
   * @note to be implemented (see Py8ShowerPSG.cc for example)
   */
  void ReCalculateSplit(node n){};

  void GetBfsSortedListOfOneToTwoNodes(vector<node> &nl);
  void GetBfsSortedListOfNodes(vector<node> &nl);
  void GetBfsSortedListOfNodesAndEdges(vector<node> &nl, vector<edge> &el);
  void GetDfsSortedListOfNodes(vector<node> &nl);

  /**
   * @brief Get the number of parent partons of the nth edge.
   * @param n Edge index.
   * @return Number of incoming edges to the source node of the nth edge.
   */
  int GetNumberOfParents(int n);

  /**
   * @brief Get the number of child partons of the nth edge.
   * @param n Edge index.
   * @return Number of outgoing edges from the target node of the nth edge.
   */
  int GetNumberOfChilds(int n);

  /**
   * @brief Retrieve all final-state partons in the shower.
   *
   * Final partons are defined as those with no children and with
   * `pstat() > -10`.
   *
   * @return Vector of shared pointers to final-state `Parton` objects.
   */
  vector<std::shared_ptr<Parton>> GetFinalPartons();

  /**
   * @brief Retrieve all final-state partons in FastJet format.
   * @return Vector of `fjcore::PseudoJet` corresponding to final-state partons.
   */
  vector<fjcore::PseudoJet> GetFinalPartonsForFastJet();

  /**
   * @brief Retrieve all partons in the shower.
   * @return Vector of shared pointers to all `Parton` objects in the shower.
   */
  vector<std::shared_ptr<Parton>> GetPartons();

  /**
   * @brief Retrieve all final edges in the shower.
   * @return Vector of GTL edge handles corresponding to final-state partons.
   */
  vector<edge> GetFinalEdges();

  /**
   * @brief Clear the list of final partons.
   */
  void ClearFinalPartonList() { pFinal.clear(); }

  /**
   * @brief Clear the list of all partons.
   */
  void ClearPartonList() { pAll.clear(); }

  /**
   * @brief Get the total number of partons (edges).
   * @return Number of edges in the graph.
   */
  int GetNumberOfPartons() const { return number_of_edges(); }

  /**
   * @brief Get the total number of vertices (nodes).
   * @return Number of nodes in the graph.
   */
  int GetNumberOfVertices() const { return number_of_nodes(); }

  /**
   * @brief Handler for saving node information in GML format.
   * @param o Output stream.
   * @param n Node to save.
   */
  void save_node_info_handler(ostream *o, node n) const;

  /**
   * @brief Handler for saving edge information in GML format.
   * @param o Output stream.
   * @param n Edge to save.
   */
  void save_edge_info_handler(ostream *o, edge n) const;

  /**
   * @brief Handler for loading edge information from GML.
   * @param e Edge to load.
   * @param read Pointer to GML key-value pairs.
   */
  void load_edge_info_handler(edge e, GML_pair *read);

  /**
   * @brief Handler for loading node information from GML.
   * @param n Node to load.
   * @param read Pointer to GML key-value pairs.
   */
  void load_node_info_handler(node n, GML_pair *read);

  /**
   * @brief Handler called before clearing the graph.
   *
   * Resets stored vertex and parton pointers to avoid dangling references.
   */
  void pre_clear_handler();

  /**
   * @brief Print all vertices in the shower.
   * @note Non-verbose mode prints to `std::cout`, verbose mode uses logger.
   */
  void PrintVertices() { PrintNodes(false); }

  /**
   * @brief Print all partons in the shower.
   * @note Non-verbose mode prints to `std::cout`, verbose mode uses logger.
   */
  void PrintPartons() { PrintEdges(false); }

  /**
   * @brief Print all nodes in the graph.
   * @param verbose If true, use logger output; otherwise print to `std::cout`.
   */
  void PrintNodes(bool verbose = true);

  /**
   * @brief Print all edges in the graph.
   * @param verbose If true, use logger output; otherwise print to `std::cout`.
   */
  void PrintEdges(bool verbose = true);

  /**
   * @brief Save the shower graph to a GML file.
   * @param fName Output filename.
   */
  void SaveAsGML(string fName) { save(fName.c_str()); }

  /**
   * @brief Save the shower graph in Graphviz DOT (`.gv`) format.
   * @param fName Output filename.
   */
  void SaveAsGV(string fName);

  /**
   * @brief Save the shower graph in GraphML XML format.
   * @param fName Output filename.
   */
  void SaveAsGraphML(string fName);

 private:
  /// Mapping from nodes to vertices.
  node_map<std::shared_ptr<Vertex>> vMap;

  /// Mapping from edges to partons.
  edge_map<std::shared_ptr<Parton>> pMap;

  std::map<std::shared_ptr<Parton>, edge> eMap;
  std::map<std::shared_ptr<Vertex>, node> nMap;

  /// Cached list of final-state partons.
  vector<std::shared_ptr<Parton>> pFinal;

  /// Cached list of all partons.
  vector<std::shared_ptr<Parton>> pAll;

  /// Cached list of final edges.
  vector<edge> eFinal;
};

}  // end namespace Jetscape
#endif
