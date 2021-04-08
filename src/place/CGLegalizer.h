/**
 * @file CGLegalizer.h
 * @brief The "legalization" solver with constraint graph + LP
 * @author Keren Zhu
 * @date 11/25/2019
 */

#ifndef IDEAPLACE_CG_LEGALIZER_H_
#define IDEAPLACE_CG_LEGALIZER_H_

#include "place/legalize/constraint_graph.hpp"
#include "place/legalize/lp_legalize.hpp"

PROJECT_NAMESPACE_BEGIN

/// @brief The LP solver for legalization
class LpLegalizeSolver {
  typedef ::klib::lp::LpModel lp_solver_type;
  typedef ::klib::lp::LpTrait lp_trait;
  typedef lp_trait::variable_type lp_variable_type;
  typedef lp_trait::expr_type lp_expr_type;

public:
  explicit LpLegalizeSolver(Database &db, Constraints &constraints,
                            bool isHor = true, IntType optHpwl = 0,
                            IntType optArea = 1)
      : _db(db), _constrains(constraints), _isHor(isHor), _optHpwl(optHpwl),
        _optArea(optArea) {} //_solver = SolverType(&_ilpModel); }
  /// @brief solve the problem
  bool solve();
  // @brief dump out the solutions to the database
  void exportSolution();
  /// @brief evaluate the objective function and return the value
  RealType evaluateObj();
  /// @brief set the maximum width or height (_wStar)
  /// @param the maximum width or height in the hpwl optimization problem
  void setWStar(RealType wStar) { _wStar = wStar; }

private:
  /// @brief solve the LP
  bool solveLp();
  /* Varibles functions */
  /// @brief add ILP variables
  void addIlpVars();
  /// @brief calculate the number of variables32
  IndexType numVars() const;
  /// @brief add location variables
  void addLocVars();
  /// @brief add wirelegth variables
  void addWirelengthVars();
  /// @brief add area variables
  void addAreaVars();
  /// @brief add sym group varibales
  void addSymVars();
  /* Obj functions */
  /// @brief set the objective function
  void configureObjFunc();
  /// @brief add wire length objective
  void addWirelengthObj();
  /// @brief add area objective
  void addAreaObj();
  /// @brief add relaxed symmetric penalty
  void addSymObj();
  /* Constraint functions */
  /// @brief add constraints
  void addIlpConstraints();
  /// @brief add area constraints
  void addBoundaryConstraints();
  /// @brief add topology constraints
  void addTopologyConstraints();
  /// @brief add symmetry constraints
  void addSymmetryConstraints();
  void addSymmetryConstraintsWithEqu();
  void addSymmetryConstraintsRex();
  /// @brief add hpwl constraints
  void addHpwlConstraints();
  /// @brief add current flow constraint
  void addCurrentFlowConstraints();

private:
  /* Configurations - Inputs */
  Database &_db;            ///< The database for the Ideaplace
  Constraints &_constrains; ///< The constraints edges to be honored
  bool _isHor = true;       ///< Whether solving horizontal or vertical
  IntType _optHpwl = 0;     ///< Whether optimizing HPWL in ILP problems
  IntType _optArea = 1;     ///< Whether optimizing area in ILP problems
  /* Optimization supporting variables */
  lp_solver_type _solver; ///<  LP sovler
  lp_expr_type _obj;      ///< The objective function of the ILP model
  std::vector<lp_variable_type>
      _locs; ///< The location variables of the ILP model
  std::vector<lp_variable_type>
      _wlL; ///< The left wirelength variables of the ILP model
  std::vector<lp_variable_type>
      _wlR;              ///< The right wirelength variables of the ILP model
  lp_variable_type _dim; ///< The variable for area optimization
  RealType _wStar = 0;   ///< The optimal W found in legalization step
  std::vector<lp_variable_type>
      _symLocs; ///< The variable for symmetric group axises
  std::vector<lp_variable_type>
      _symRexLeft; ///< The variables representing the left extrems of sym axis
                   ///< of each group
  std::vector<lp_variable_type>
      _symRexRight; ///< The variables representing the left extrems of sym axis
                    ///< of each group
#ifdef MULTI_SYM_GROUP
  bool _isMultipleSymGrp = true;
#else
  bool _isMultipleSymGrp = false;
#endif
  bool _relaxEqualityConstraint = false;
  bool _useCurrentFlowConstraint = false;
  // SolverType _solver; ///< Solver
  /*  Optimization Results */
  RealType _largeNum = 900000.0; ///< A large number
};

class CGLegalizer {
private:
  class BoxEdge {
  public:
    explicit BoxEdge(LocType coord_, IndexType cellIdx_, bool isTop_)
        : coord(coord_), cellIdx(cellIdx_), isTop(isTop_) {}
    bool operator<(const BoxEdge &rhs) const {
      if (this->coord == rhs.coord) {
        if (this->isTop == rhs.isTop) {
          return this->cellIdx < rhs.cellIdx;
        } else {
          if (!this->isTop) {
            return false;
          } else {
            return true;
          }
        }
      } else {
        return this->coord < rhs.coord;
      }
    }
    std::string toStr() const {
      std::stringstream ss;
      ss << "LocType " << coord << " cellIdx " << cellIdx << " isTop " << isTop;
      return ss.str();
    }

  public:
    LocType coord;     ///< Coordinate of the edge
    IndexType cellIdx; ///< The index of the cell
    bool isTop;        ///< True: top/right edge. False: bottom/left edge
  };

public:
  /// @brief Constructor
  /// @param The database of IdeaPlaceEx
  explicit CGLegalizer(Database &db) : _db(db) {}
  /// @brief legalize the design
  bool legalize();

private:
  /// @brief Generate the constraints (not optimal in number of constraints).
  /// Based on sweeping algorithm
  void generateHorConstraints();
  void generateVerConstraints();
  /// @brief linear programming-based legalization
  /// @param if solving horizontal or vertical
  /// @return the resulting objective function. if negative, failed
  RealType lpLegalization(bool isHor);
  /// @brief LP-based detailed placement. For optimizing wire length
  bool lpDetailedPlacement();

private:
  Database &_db;             ///< The database of IdeaPlaceEx
  Constraints _hConstraints; ///< The horizontal constraint edges
  Constraints _vConstraints; ///< The vertical constraint edges
  RealType _wStar; ///< The width from the objective function of the first LP
  RealType _hStar; ///< The width from the objective function of the first LP
};

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_CG_LEGALIZER_H_
