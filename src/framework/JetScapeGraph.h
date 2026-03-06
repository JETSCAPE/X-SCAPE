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

#ifndef JETSCAPEGRAPH_H
#define JETSCAPEGRAPH_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
//#include "JetScapeParticles.h"
#include "JetScapeLogger.h"
//#include "helper.h"

// REMARK: Still keeping the orginal PartonShower class for now. Even though
//    in the future it should be derived from the templated based class!!!

namespace Jetscape {

class Vertex;
class Parton;
class JetScapeParticleBase;
class Hadron;

template <class T>
/**
 * @brief Graph wrapper used by the JETSCAPE framework to store
 *        vertices and particles (templated by particle type).
 *
 * This class adapts the underlying GTL `graph` to hold
 * `Vertex` objects as nodes and objects of type `T` as edge-associated
 * particles. It provides convenience accessors, traversal helpers,
 * and serialization hooks used throughout the framework.
 *
 * @tparam T Type used to represent particles stored on edges (e.g.
 *           `JetScapeParticleBase` or `Hadron`).
 */
class JetScapeGraph : public graph {
 public:
  /**
   * @brief Construct an empty `JetScapeGraph`.
   */
  JetScapeGraph<T>();

  /**
   * @brief Virtual destructor.
   */
  virtual ~JetScapeGraph<T>();

  /**
   * @brief Add a new vertex to the graph.
   *
   * @param v Shared pointer to the `Vertex` object to insert.
   * @return The graph `node` identifier for the newly inserted vertex.
   */
  node new_vertex(std::shared_ptr<Vertex> v);

  /**
   * @brief Create a new particle (edge) between two nodes.
   *
   * @param s Source node identifier.
   * @param t Target node identifier.
   * @param p Shared pointer to the particle object to attach to the edge.
   * @return Integer status or the created edge id (implementation-specific).
   */
  int new_particle(node s, node t, std::shared_ptr<T> p);

  // virtual std::unique_ptr<JetScapeGraph> Clone(); //not yet working for 2--2
  // !!! Fix it !!!

  /**
   * @brief Insert a particle on an existing edge and associate a vertex.
   *
   * @param e Edge where the particle should be inserted.
   * @param v Vertex to associate with the insertion.
   * @param p Particle to insert.
   */
  void InsertParticle(edge e, std::shared_ptr<Vertex> v, std::shared_ptr<T> p);

  /**
   * @brief Insert a particle after a given edge, associating a vertex.
   *
   * @param e Edge after which to insert the new particle.
   * @param v Vertex to associate with the new particle.
   * @param p Particle to insert.
   */
  void InsertParticleAfter(edge e, std::shared_ptr<Vertex> v,
                           std::shared_ptr<T> p);

  /**
   * @brief Placeholder for edge insertion logic.
   *
   * Current implementation is empty; provided for API compatibility.
   */
  void InsertEdge(edge e, edge eIns){};

  /**
   * @brief Get the `Vertex` associated with a node.
   *
   * @param n Node identifier.
   * @return Shared pointer to the `Vertex` stored at node `n`.
   */
  std::shared_ptr<Vertex> GetVertex(node n) { return vMap[n]; }

  /**
   * @brief Get the particle stored on an edge.
   *
   * @param e Edge identifier.
   * @return Shared pointer to the particle of type `T` stored on edge `e`.
   */
  std::shared_ptr<T> GetParticle(edge e) { return pMap[e]; }

  /**
   * @brief Get the n-th particle in the internal ordering.
   *
   * @param n Index of the particle to retrieve.
   * @return Shared pointer to the particle at index `n`.
   */
  std::shared_ptr<T> GetParticleAt(int n);

  /**
   * @brief Get the n-th vertex in the internal ordering.
   *
   * @param n Index of the vertex to retrieve.
   * @return Shared pointer to the vertex at index `n`.
   */
  std::shared_ptr<Vertex> GetVertexAt(int n);

  /**
   * @brief Lookup the edge corresponding to a given particle object.
   *
   * @param p Shared pointer to a particle.
   * @return Edge identifier associated with `p`.
   */
  edge GetEdge(std::shared_ptr<T> p) { return eMap.at(p); }

  /**
   * @brief Lookup the node corresponding to a given vertex object.
   *
   * @param v Shared pointer to a vertex.
   * @return Node identifier associated with `v`.
   */
  node GetNode(std::shared_ptr<Vertex> v) { return nMap.at(v); }

  /**
   * @brief Return the node at a given index.
   *
   * @param n Index of the node.
   * @return Node identifier at index `n`.
   */
  node GetNodeAt(int n);

  /**
   * @brief Return the edge at a given index.
   *
   * @param n Index of the edge.
   * @return Edge identifier at index `n`.
   */
  edge GetEdgeAt(int n);

  /**
   * @brief Return true if the node has N incoming edges and one outgoing.
   *
   * @param n Node identifier.
   * @return True if node is N->1 topology.
   */
  bool IsNtoOne(node n);

  /**
   * @brief Return true if the node has exactly one incoming and one outgoing.
   *
   * @param n Node identifier.
   * @return True if node is 1->1 topology.
   */
  bool IsOneToOne(node n);

  /**
   * @brief Return true if the node represents a 1->2 splitting.
   *
   * @param n Node identifier.
   * @return True if node is a one-to-two split.
   */
  bool IsOneToTwo(node n);

  /**
   * @brief Return true if the node is an end node (no outgoing edges).
   *
   * @param n Node identifier.
   * @return True if node has no children.
   */
  bool IsEndNode(node n);

  /**
   * @brief Return whether the given edge corresponds to a "high" split.
   *
   * Implementation-specific: used to decide which branch corresponds to the
   * higher-energy or higher-kt daughter.
   *
   * @param e Edge identifier.
   * @return True if edge is marked as the high-split.
   */
  bool IsHighSplitEdge(edge e);

  /**
   * @brief Return whether the node contains a high-split edge.
   *
   * @param n Node identifier.
   * @return True if node contains a high-split.
   */
  bool IsHighSplitNode(node n);

  /**
   * @brief Return delta-R for the splitting at node `n`.
   *
   * @param n Node identifier.
   * @return Delta-R between split daughters.
   */
  double GetSplitDeltaR(node n);

  /**
   * @brief Return the momentum fraction z for the splitting at node `n`.
   *
   * @param n Node identifier.
   * @return The splitting momentum fraction z.
   */
  double GetSplitZ(node n);

  /**
   * @brief Return transverse momentum (kT) for splitting at node `n`.
   *
   * @param n Node identifier.
   * @return kT of the splitting.
   */
  double GetSplitKt(node n);

  /**
   * @brief Get the time of the next node in chronological order from `n`.
   *
   * @param n Node identifier.
   * @return Time value (implementation-specific units).
   */
  double GetNextNodeTime(node n);

  /**
   * @brief Get the time associated with node `n`.
   *
   * @param n Node identifier.
   * @return Time stored in the vertex at node `n`.
   */
  double GetNodeTime(node n);  // {return GetVertex(n)->x_in().t();}

  /**
   * @brief Return the edge corresponding to the high split for node `n`.
   *
   * @param n Node identifier.
   * @return Edge identifier for the high-split branch.
   */
  edge GetHighSplitEdge(node n);

  /**
   * @brief Return the edge corresponding to the low split for node `n`.
   *
   * @param n Node identifier.
   * @return Edge identifier for the low-split branch.
   */
  edge GetLowSplitEdge(node n);

  /**
   * @brief Recalculate splitting information for node `n`.
   *
   * Currently a placeholder; platform-specific implementations such as
   * Py8ShowerPSG.cc provide examples.
   *
   * @param n Node identifier.
   */
  void ReCalculateSplit(
      node n){};  // to be implemented (see Py8ShowerPSG.cc for example)

  /**
   * @brief Populate `nl` with nodes that are 1->2 splits in BFS order.
   *
   * @param nl Output vector to receive node identifiers.
   */
  void GetBfsSortedListOfOneToTwoNodes(vector<node> &nl);

  /**
   * @brief Populate `nl` with nodes in BFS (breadth-first) order.
   *
   * @param nl Output vector to receive node identifiers.
   */
  void GetBfsSortedListOfNodes(vector<node> &nl);

  /**
   * @brief Populate `nl` and `el` with nodes and edges in BFS order.
   *
   * @param nl Output vector for nodes.
   * @param el Output vector for edges.
   */
  void GetBfsSortedListOfNodesAndEdges(vector<node> &nl, vector<edge> &el);

  /**
   * @brief Populate `nl` with nodes in DFS (depth-first) order.
   *
   * @param nl Output vector to receive node identifiers.
   */
  void GetDfsSortedListOfNodes(vector<node> &nl);

  /**
   * @brief Return number of parents (incoming edges) for node `n`.
   *
   * @param n Node identifier.
   * @return Number of parent nodes.
   */
  int GetNumberOfParents(int n);

  /**
   * @brief Return number of children (outgoing edges) for node `n`.
   *
   * @param n Node identifier.
   * @return Number of child nodes.
   */
  int GetNumberOfChilds(int n);

  /**
   * @brief Return a vector of final particles (no further daughters).
   *
   * @return Vector of shared pointers to final particles of type `T`.
   */
  vector<std::shared_ptr<T>> GetFinalParticles();

  /**
   * @brief Convert final particles into FastJet `PseudoJet` objects.
   *
   * @return Vector of `fjcore::PseudoJet` representing final particles.
   */
  vector<fjcore::PseudoJet> GetFinalParticlesForFastJet();

  /**
   * @brief Clear the internal list of final particles.
   */
  void ClearFinalParticleList() { pFinal.clear(); }

  /**
   * @brief Return number of parton edges stored in the graph.
   *
   * @return Number of edges (partons).
   */
  int GetNumberOfPartons() const { return number_of_edges(); }

  /**
   * @brief Return number of vertices stored in the graph.
   *
   * @return Number of nodes (vertices).
   */
  int GetNumberOfVertices() const { return number_of_nodes(); }

  /**
   * @brief Serialization hook used by the GTL save routine: write node info.
   *
   * @param o Output stream to write to.
   * @param n Node identifier whose info should be written.
   */
  void save_node_info_handler(ostream *o, node n) const;

  /**
   * @brief Serialization hook used by the GTL save routine: write edge info.
   *
   * @param o Output stream to write to.
   * @param n Edge identifier whose info should be written.
   */
  void save_edge_info_handler(ostream *o, edge n) const;

  /**
   * @brief Deserialization hook: load edge info from a GML pair.
   *
   * @param e Edge identifier to populate.
   * @param read Pointer to GML_pair containing data read from file.
   */
  void load_edge_info_handler(edge e, GML_pair *read);

  /**
   * @brief Deserialization hook: load node info from a GML pair.
   *
   * @param n Node identifier to populate.
   * @param read Pointer to GML_pair containing data read from file.
   */
  void load_node_info_handler(node n, GML_pair *read);

  /**
   * @brief Hook called before clearing the graph; allows cleanup.
   */
  void pre_clear_handler();

  /**
   * @brief Convenience wrapper to print vertices (nodes).
   */
  void PrintVertices() { PrintNodes(false); }

  /**
   * @brief Convenience wrapper to print partons (edges).
   */
  void PrintPartons() { PrintEdges(false); }

  /**
   * @brief Print node information to stdout.
   *
   * @param verbose When true, print extended information.
   */
  void PrintNodes(bool verbose = true);

  /**
   * @brief Print edge information to stdout.
   *
   * @param verbose When true, print extended information.
   */
  void PrintEdges(bool verbose = true);

  /**
   * @brief Save graph in GML format to file `fName`.
   *
   * @param fName Output filename.
   */
  void SaveAsGML(string fName) { save(fName.c_str()); }

  /**
   * @brief Save graph in GraphViz (GV) format to file `fName`.
   *
   * @param fName Output filename.
   */
  void SaveAsGV(string fName);

  /**
   * @brief Save graph in GraphML format to file `fName`.
   *
   * @param fName Output filename.
   */
  void SaveAsGraphML(string fName);

 private:
  /**
   * @brief Map from node id to stored `Vertex` pointer.
   */
  node_map<std::shared_ptr<Vertex>> vMap;

  /**
   * @brief Map from edge id to stored particle pointer of type `T`.
   */
  edge_map<std::shared_ptr<T>> pMap;

  // REMARK: check if needed clearing in destructor ...
  /**
   * @brief Reverse lookup map from particle pointer to edge id.
   */
  std::map<std::shared_ptr<T>, edge> eMap;

  /**
   * @brief Reverse lookup map from vertex pointer to node id.
   */
  std::map<std::shared_ptr<Vertex>, node> nMap;

  /**
   * @brief Cached list of final particles (particles with no daughters).
   */
  vector<std::shared_ptr<T>> pFinal;
};

// typedef JetScapeGraph<JetScapeParticleBase> EventGraph;
// typedef JetScapeGraph<Hadron> HadronGraph;

}  // end namespace Jetscape
#endif
