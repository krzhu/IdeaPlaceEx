/**
 * @file CGLegalizer.h
 * @brief The "legalization" solver with constraint graph + LP
 * @author Keren Zhu
 * @date 11/25/2019
 */

#ifndef IDEAPLACE_CG_LEGALIZER_H_
#define IDEAPLACE_CG_LEGALIZER_H_

#include "ConstraintGraph.h"
#include "db/Database.h"

PROJECT_NAMESPACE_BEGIN

/// @brief a directed edge representing a constraint
class ConstraintEdge
{
    public:
        explicit ConstraintEdge(IndexType source, IndexType target, IntType weight):
            _source(source), _target(target), _weight(weight)
    {}
        /// @brief Get the index of source vertex
        /// @return the index of source vertex
        IndexType source() const { return _source; }
        /// @brief Get the index of target vertex
        /// @return the index of target vertex
        IndexType target() const { return _target; }
        /// @brief Get the weight of this edge
        /// @return the weight of this edge
        IntType weight() const { return _weight; }
        /// @brief to string for debuging
        std::string toStr() const 
        {
            std::stringstream ss;
            ss << "source "<< _source <<" target "<< _target << " weight "<< _weight;
            return ss.str();
        }
    private:
        IndexType _source;  ///< The index of source vertex
        IndexType _target; ///< The index of target vertex
        IntType _weight;  ///< The weight of this edge
};


/// @brief Just a wrappers for ConstraintEdge
class Constraints
{

    public:
        /// @brief default constructor
        explicit Constraints() = default;
        /// @breif clear the constraint edges
        void clear() { _edges.clear(); }
        /// @brief get the constraint edges
        /// @return the constraint edges
        const std::vector<ConstraintEdge> & edges() const { return _edges; }
        /// @brief get the constraint edges
        /// @return the constraint edges
        std::vector<ConstraintEdge> & edges() { return _edges; }
        /// @brief add a constraint edge
        /// @param the source cell index
        /// @param the target cell index
        /// @param the weight of the edge
        void addConstraintEdge(IndexType sourceIdx, IndexType targetIdx, IntType weight)
        {
            _edges.emplace_back(ConstraintEdge(sourceIdx, targetIdx, weight));
        }
        
    private:
        std::vector<ConstraintEdge> _edges; ///< The constraint edges
};

class CGLegalizer
{
    private:
        class BoxEdge
        {
            public:
                explicit BoxEdge(LocType coord_, IndexType cellIdx_, bool isTop_)
                    : coord(coord_), cellIdx(cellIdx_), isTop(isTop_) {}
                bool operator<(const BoxEdge &rhs) const
                {
                    if (this->coord == rhs.coord)
                    {
                        if (this->isTop == rhs.isTop)
                        {
                            return this->cellIdx < rhs.cellIdx;
                        }
                        else
                        {
                            if (!this->isTop)
                            {
                                return false;
                            }
                            else
                            {
                                return true;
                            }
                        }
                    }
                    else
                    {
                        return this->coord < rhs.coord;
                    }
                }
            public:
                LocType coord; ///< Coordinate of the edge
                IndexType cellIdx; ///< The index of the cell 
                bool isTop; ///< True: top/right edge. False: bottom/left edge

        };
    public:
        /// @brief Constructor
        /// @param The database of IdeaPlaceEx
        explicit CGLegalizer(Database &db) : _db(db) {}
        /// @brief legalize the design
        void legalize()
        {
            this->generateConstraints();
        }
    private:
        /// @brief Generate the constraints (not optimal in number of constraints). Based on sweeping algorithm
        void generateConstraints();
        /// @brief construct constraint graph from two constraints
        void constructConstraintGraphs();
        /// @brief perform DFS-based transitive reduction
        /// @param the constraint graph
        void dagTransitiveReduction(ConstraintGraph & cg);
        /// @brief remove transitive edge from the graph
        /// @param constraint graph
        /// @param Vector2D of edge matrix
        /// @param number of cells/vertices?
        /// @param vector of whether visited
        /// @param Vector2D of reachable
        /// @return has transitive edge
        bool dfsRemoveTransitiveEdge(ConstraintGraph &cg, Vector2D<IntType> &edgeMat, IndexType node, 
                std::vector<bool> &visited, Vector2D<IntType> &reachable, ConstraintGraph::IndexMap &idxMap);
        /// @brief delete edges in initializing irredundant edges
        /// @param is horizontal? false-> vertical
        /// @param orders
        /// @param cands
        /// @param index of cell
        /// @param the overlap any vector
        void initIrredundantEdgesDelete(bool isHor, std::vector<IndexType> &orders,
                std::vector<IntType> &cand, IndexType cellIdx, std::vector<bool> &overlapAny);
        /// @brief insert edges in initializing irredundant edges
        /// @param is horizontal? false -> vertical
        /// @param orders
        /// @param cands
        /// @param index of cell
        void initIrredundantEdgesInsert(bool isHor , std::vector<IndexType> &orders, 
                std::vector<IntType> &cand, IndexType cellIdx);
        /// @brief get necessary edges
        void getNecessaryEdges();
        /// @brief perform DFS on the graph
        /// @param first: dp table
        /// @param second: the visited vector
        /// @param third: the node index
        /// @param fourth: constraint graph
        /// @param fifth: IndexMap for the boost graph
        void dfsGraph(Vector2D<IntType>& dpTab, std::vector<IntType> &visited, IndexType nodeIdx, ConstraintGraph &cg,
                ConstraintGraph::IndexMap &idxMap);
        /// @brief add edge greedily
        /// @param first: node index i
        /// @param second: node index j
        void addEdgeGreedy(IndexType i, IndexType j);
        /// @brief reload the constraints from the boost-based constraint graph
        void readloadConstraints();
    private:
        Database &_db; ///< The database of IdeaPlaceEx
        ConstraintGraph _hCG; ///< The horizontal constraint graph
        ConstraintGraph _vCG; ///< The vertical constraint graph
        Constraints _hConstraints; ///< The horizontal constraint edges
        Constraints _vConstraints; ///< The vertical constraint edges
};


PROJECT_NAMESPACE_END

#endif //IDEAPLACE_CG_LEGALIZER_H_
