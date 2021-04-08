/**
 * @file constraint_graph.h
 * @brief The constraint graph implementation with boost graph and simple related data structures
 * @author Keren Zhu
 * @date 11/25/2019
 */

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include "global/global.h"

PROJECT_NAMESPACE_BEGIN

class ConstraintGraph {
public:
  /// @brief default constructor
  explicit ConstraintGraph() { _cg.clear(); }
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS,
                                boost::no_property,
                                boost::property<boost::edge_weight_t, IntType>>
      graph_t;
  typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<graph_t>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<graph_t>::edge_iterator edge_iterator;
  typedef boost::property_map<graph_t, boost::vertex_index_t>::type IndexMap;
  typedef boost::graph_traits<graph_t>::adjacency_iterator adjacency_iterator;
  /// @brief return the boost graph
  graph_t &boostGraph() { return _cg; }
  /// @brief construct the graph with number of vertices
  /// @param the number of vertices
  void allocateVertices(IndexType numVex) {
    _cg.clear();
    _cg = graph_t(numVex);
  }
  /// @brief get the number of nodes
  /// @return the number of nodes
  IndexType numNodes() const { return boost::num_vertices(_cg); }
  /// @brief get the source node index
  /// @return the index of the source node
  IndexType sourceNodeIdx() const { return numNodes() - 2; }
  /// @brief get the target node index
  /// @return the index of the target node
  IndexType targetNodeIdx() const { return numNodes() - 1; }
  /// @brief get the number of cells
  /// @return the number of cell nodes in the graph
  IndexType numCellNodes() const { return numNodes() - 2; }
  /// @brief add edge to the graph
  /// @param the source node index
  /// @param the target node index
  void addEdge(IndexType sourceIdx, IndexType targetIdx, IntType weight = 1) {
    boost::add_edge(boost::vertex(sourceIdx, _cg),
                    boost::vertex(targetIdx, _cg), weight, _cg);
  }
  /// @brief remove a edge from the graph
  /// @param the source index
  /// @param the target index
  void removeEdge(IndexType sourceIdx, IndexType targetIdx) {
    boost::remove_edge(boost::vertex(sourceIdx, _cg),
                       boost::vertex(targetIdx, _cg), _cg);
  }
  /// @brief determine whether the graph has one specific edge
  /// @param the source index of the edge
  /// @param the target index of the edge
  /// @return true if has edge. false if not
  bool hasEdge(IndexType sourceIdx, IndexType targetIdx) {
    auto edge = boost::edge(boost::vertex(sourceIdx, _cg),
                            boost::vertex(targetIdx, _cg), _cg);
    return edge.second;
  }

  void clear() { _cg.clear(); }

private:
  graph_t _cg; ///< The boost graph
};


/// @brief a directed edge representing a constraint
class ConstraintEdge {
public:
  explicit ConstraintEdge(IndexType source, IndexType target)
      : _source(source), _target(target) {}
  /// @brief Get the index of source vertex
  /// @return the index of source vertex
  IndexType source() const { return _source; }
  /// @brief Get the index of target vertex
  /// @return the index of target vertex
  IndexType target() const { return _target; }
  /// @brief Get the weight of this edge
  /// @return the weight of this edge
  IntType weight() const { return 1; }
  /// @brief to string for debuging
  std::string toStr() const {
    std::stringstream ss;
    ss << "source " << _source << " target " << _target;
    return ss.str();
  }
  bool operator<(const ConstraintEdge &rhs) const {
    if (_source == rhs.source()) {
      return _target < rhs.target();
    }
    return _source < rhs.source();
  }
  bool operator==(const ConstraintEdge &rhs) const {
    return _source == rhs.source() && _target == rhs.target();
  }

private:
  IndexType _source; ///< The index of source vertex
  IndexType _target; ///< The index of target vertex
                     // IntType _weight;  ///< The weight of this edge
};

/// @brief Just a wrappers for ConstraintEdge
class Constraints {

public:
  /// @brief default constructor
  explicit Constraints() = default;
  /// @breif clear the constraint edges
  void clear() { _edges.clear(); }
  /// @brief get the constraint edges
  /// @return the constraint edges
  const std::set<ConstraintEdge> &edges() const { return _edges; }
  /// @brief get the constraint edges
  /// @return the constraint edges
  std::set<ConstraintEdge> &edges() { return _edges; }
  /// @brief add a constraint edge
  /// @param the source cell index
  /// @param the target cell index
  /// @param the weight of the edge
  void addConstraintEdge(IndexType sourceIdx, IndexType targetIdx,
                         IntType weight) {
    _edges.insert(ConstraintEdge(sourceIdx, targetIdx));
  }
  bool hasEdgeNoDirection(IndexType sourceIdx, IndexType targetIdx) {
    auto it = _edges.find(ConstraintEdge(sourceIdx, targetIdx));
    if (it != _edges.end()) {
      return true;
    }
    it = _edges.find(ConstraintEdge(targetIdx, sourceIdx));
    if (it != _edges.end()) {
      return true;
    }
    return false;
  }
  void removeConstraintEdge(IndexType sourceIdx, IndexType targetIdx) {
    auto it = _edges.find(ConstraintEdge(sourceIdx, targetIdx));
    if (it != _edges.end()) {
      _edges.erase(it);
    }
    it = _edges.find(ConstraintEdge(targetIdx, sourceIdx));
    if (it != _edges.end()) {
      _edges.erase(it);
    }
  }

private:
  std::set<ConstraintEdge> _edges; ///< The constraint edges
};

PROJECT_NAMESPACE_END

