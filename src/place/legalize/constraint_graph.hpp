/**
 * @file constraint_graph.h
 * @brief The constraint graph implementation with boost graph and simple related data structures
 * @author Keren Zhu
 * @date 11/25/2019
 */

#pragma once

#include "global/global.h"

PROJECT_NAMESPACE_BEGIN


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

