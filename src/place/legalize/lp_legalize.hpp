/**
 * @file lp_legalize.h
 * @brief Legalize the placement using linear programming
 * @author Keren Zhu
 * @date 04/07/2021
 */

#pragma once


#include "db/Database.h"
#include "constraint_graph.hpp"
#include "util/linear_programming.h"

PROJECT_NAMESPACE_BEGIN

namespace lp_legalize {

  struct LEGALIZE_HORIZONTAL_DIRECTION {};
  struct LEGALIZE_VERTICAL_DIRECTION {};
  struct RELAX_SYM_CONSTR {}; // Relax the symmetric constraints 
  struct DO_NOT_RELAX_SYM_CONSTR {}; // Do not relax symmetric constraints

namespace _lp_legalize_details {
  template<typename is_hor_type>
    struct cell_loc_trait {};
  template<>
    struct cell_loc_trait<LEGALIZE_HORIZONTAL_DIRECTION> {
      static LocType getLen(const Database &db, IndexType cellIdx) {
        if (cellIdx < db.numCells()) {
          return db.cell(cellIdx).cellBBox().xLen();
        }
        auto wellRectIdx = db.getWellRectIdx(cellIdx - db.numCells());
        return db.well(wellRectIdx.first).rects().at(wellRectIdx.second).xLen();
      }
      static LocType getSpacing(const Database &db, IndexType sourceCellIdx, IndexType targetCellIdx) {
        // Source cell is on the left
        return db.cellSpacing(sourceCellIdx, targetCellIdx).xLo();
      }
      static void moveWell(Well &well, LocType offset) {
        well.moveBy(offset, 0);
      }
      static LocType getCellLL(const Database &db, IndexType cellIdx) {
        return db.cell(cellIdx).xLo();
      }
      static LocType getCellLL(const Cell &cell) {
        return cell.xLo();
      }
      static LocType getPinMidLoc(const Database &db, IndexType pinIdx) {
        return db.pinOffsetToCell(pinIdx).x();
      }
      static void setCellLL(Database &db, IndexType cellIdx, LocType loc) {
        db.cell(cellIdx).setXLo(loc);
      }
      static LocType getNetVirtualPinLoc(const Database &db, IndexType netIdx) {
        return db.net(netIdx).virtualPinLoc().x();
      }
      static LocType getWellRectLL(const Well &well, IndexType rectIdx) {
        return well.rects().at(rectIdx).xLo();
      }
      static LocType getWellRectLen(const Well &well, IndexType rectIdx) {
        return well.rects().at(rectIdx).xLen();
      }
    };
  template<>
    struct cell_loc_trait<LEGALIZE_VERTICAL_DIRECTION> {
      static LocType getLen(const Database &db, IndexType cellIdx) {
        if (cellIdx < db.numCells()) {
          return db.cell(cellIdx).cellBBox().yLen();
        }
        auto wellRectIdx = db.getWellRectIdx(cellIdx - db.numCells());
        return db.well(wellRectIdx.first).rects().at(wellRectIdx.second).yLen();
      }
      static LocType getSpacing(const Database &db, IndexType sourceCellIdx, IndexType targetCellIdx) {
        // Source cell is on the bottom
        return db.cellSpacing(sourceCellIdx, targetCellIdx).yLo();
      }
      static void moveWell(Well &well, LocType offset) {
        well.moveBy(0, offset);
      }
      static LocType getCellLL(const Database &db, IndexType cellIdx) {
        return db.cell(cellIdx).yLo();
      }
      static LocType getCellLL(const Cell &cell) {
        return cell.yLo();
      }
      static LocType getPinMidLoc(const Database &db, IndexType pinIdx) {
        return db.pinOffsetToCell(pinIdx).y();
      }
      static void setCellLL(Database &db, IndexType cellIdx, LocType loc) {
        db.cell(cellIdx).setYLo(loc);
      }
      static LocType getNetVirtualPinLoc(const Database &db, IndexType netIdx) {
        return db.net(netIdx).virtualPinLoc().y();
      }
      static LocType getWellRectLL(const Well &well, IndexType rectIdx) {
        return well.rects().at(rectIdx).yLo();
      }
      static LocType getWellRectLen(const Well &well, IndexType rectIdx) {
        return well.rects().at(rectIdx).yLen();
      }
    };

  struct hpwl_trait {
    template<typename lp_type>
      static void addObj(lp_type &lp) {
        for (IndexType netIdx = 0; netIdx < lp._db.numNets(); ++netIdx) {
          if (lp._db.net(netIdx).numPinIdx() == 0) {
            continue;
          }
          if (lp._db.net(netIdx).numPinIdx() == 1 &&
              !lp._db.net(netIdx).isValidVirtualPin()) {
            continue;
          }
          auto weight = lp._db.net(netIdx).weight();
          lp._obj += weight * (lp._wlR[netIdx] - lp._wlL[netIdx]);
        }
      }
    template<typename lp_type>
      static void addVars(lp_type &lp) {
        lp._wlL.resize(lp._db.numNets());
        lp._wlR.resize(lp._db.numNets());
        for (IndexType i = 0; i < lp._db.numNets(); ++i) {
          lp._wlL.at(i) = lp_type::lp_trait::addVar(lp._solver, "WLL"+std::to_string(i));
          lp._wlR.at(i) = lp_type::lp_trait::addVar(lp._solver, "WLR"+std::to_string(i));
        }
      }
    template<typename lp_type>
      static void addConstr(lp_type &lp) {
        for (IndexType netIdx = 0; netIdx < lp._db.numNets(); ++netIdx) {
          const auto &net = lp._db.net(netIdx);
          for (IndexType pinIdxInNet = 0; pinIdxInNet < net.numPinIdx();
               ++pinIdxInNet) {
            IndexType pinIdx = net.pinIdx(pinIdxInNet);
            const auto &pin = lp._db.pin(pinIdx);
            RealType loc = static_cast<RealType>(cell_loc_trait<typename lp_type::is_hor_type>::getPinMidLoc(lp._db, pinIdx));
            // wl_l <= _loc + pin_offset for all pins in the net
            lp_type::lp_trait::addConstr(lp._solver,
                                lp._wlL.at(netIdx) - lp._locs.at(pin.cellIdx()) <= loc);
            // wl_r >= _loc + pin_offset for all pins in the net
            lp_type::lp_trait::addConstr(lp._solver,
                                lp._wlR.at(netIdx) - lp._locs.at(pin.cellIdx()) >= loc);
            lp_type::lp_trait::addConstr(lp._solver, lp._wlR.at(netIdx) - lp._wlL.at(netIdx) >= 0);
          }
          // Wirelength with virtual pin
          if (net.isValidVirtualPin()) {
            RealType loc = static_cast<RealType>(cell_loc_trait<typename lp_type::is_hor_type>::getNetVirtualPinLoc(lp._db, netIdx));
            // wl_l <= _loc + pin_offset for all pins in the net
            lp_type::lp_trait::addConstr(lp._solver, lp._wlL.at(netIdx) <= std::max(loc, 0.0));
            // wl_r >= _loc + pin_offset for all pins in the net
            lp_type::lp_trait::addConstr(lp._solver, lp._wlR.at(netIdx) >= loc);
          }
        }
      }
    };

  struct area_trait {
    template<typename lp_type>
      static void addObj(lp_type &lp) {
        lp._obj += lp._boundary;
      }
    template<typename lp_type>
      static void addVars(lp_type &lp) {
        lp._boundary = lp_type::lp_trait::addVar(lp._solver, "boundary");
      }
    template<typename lp_type>
      static void addConstr(lp_type &lp) {
        for (IndexType cellIdx = 0; cellIdx < lp._db.numCells(); ++cellIdx) {
        // 0 <= x_i <= W* - w_i
          lp_type::lp_trait::addConstr(
              lp._solver, lp._locs.at(cellIdx) - lp._boundary <=  -
                                cell_loc_trait<typename lp_type::is_hor_type>::getLen(lp._db, cellIdx));
          lp_type::lp_trait::addConstr(lp._solver,
                              lp._locs.at(cellIdx) >= 0);
        }
      }
    template<typename lp_type>
      static void fixBoundary(lp_type &lp, RealType boundary) {
        // Here shall can use either equality or <=.
        lp_type::lp_trait::addConstr(
            lp._solver, lp._boundary <=  boundary
            );
      }
  };

  /// @brief the variable type of symmetric constraint
  template<typename lp_variable_type, typename relax_sym_type>
    struct sym_variable_type {};
  template<typename lp_variable_type>
    struct sym_variable_type<lp_variable_type, RELAX_SYM_CONSTR> {
      std::vector<lp_variable_type> _symRexLeft;
      std::vector<lp_variable_type> _symRexRight;
      static constexpr RealType _largeNum = 900000.0; ///< A large number
    };
  template<typename lp_variable_type>
    struct sym_variable_type<lp_variable_type, DO_NOT_RELAX_SYM_CONSTR> {
      std::vector<lp_variable_type> _symLocs;
    };

  /// @brief Relaxation of symmetric constraints
  template<typename is_hor_type>
    struct _sym_rex_trait {};

  template<>
    struct _sym_rex_trait<LEGALIZE_HORIZONTAL_DIRECTION> {
      template<typename lp_type>
        static void addObj(lp_type &lp) {
          for (IndexType symGrpIdx = 0; symGrpIdx < lp._symVars._symRexLeft.size(); ++symGrpIdx) {
            lp._obj +=
                lp._symVars._largeNum * (lp._symVars._symRexRight.at(symGrpIdx) - lp._symVars._symRexLeft.at(symGrpIdx));
          }
        }
      template<typename lp_type>
        static void addVars(lp_type &lp) {
#ifdef MULTI_SYM_GROUP
          // Symmetric group axis variables
          lp._symVars._symRexLeft.resize(lp._db.numSymGroups());
          lp._symVars._symRexRight.resize(lp._db.numSymGroups());
          for (IndexType i = 0; i < lp._db.numSymGroups(); ++i) {
            lp._symVars._symRexLeft.at(i) = lp_type::lp_trait::addVar(lp._solver, "symRexLeft" + std::to_string(i));
            lp._symVars._symRexRight.at(i) = lp_type::lp_trait::addVar(lp._solver, "symRexRight" + std::to_string(i));
          }
#else
            lp._symVars._symRexLeft.resize(1);
            lp._symVars._symRexLeft[0] = lp_type::lp_trait::addVar(lp._solver, "symRexLeft");
            lp._symVars._symRexRight.resize(1);
            lp._symVars._symRexRight[0] = lp_type::lp_trait::addVar(lp._solver, "symRexRight");
#endif
        }
      template<typename lp_type>
        static void addConstr(lp_type &lp) {
          for (IndexType symGroupIdx = 0; symGroupIdx < lp._db.numSymGroups();
               ++symGroupIdx) {
            const auto &symGroup = lp._db.symGroup(symGroupIdx);
#ifdef MULTI_SYM_GROUP
            typename lp_type::lp_variable_type &leftSymLoc = lp._symVars._symRexLeft.at(symGroupIdx);
            typename lp_type::lp_variable_type &rightSymLoc = lp._symVars._symRexRight.at(symGroupIdx);
#else
            typename lp_type::lp_variable_type &leftSymLoc = lp._symVars._symRexLeft.at(0);
            typename lp_type::lp_variable_type &rightSymLoc = lp._symVars._symRexRight.at(0);
#endif
            // Right >= Left
            lp_type::lp_trait::addConstr(lp._solver, rightSymLoc - leftSymLoc >= 0);
            for (IndexType symPairIdx = 0; symPairIdx < symGroup.numSymPairs();
                 ++symPairIdx) {
              const auto &symPair = symGroup.symPair(symPairIdx);
              // x1 + x2 + width  <= right * 2
              // x1 + x2 + width >= left * 2
              typename lp_type::lp_variable_type x1_ = lp._locs.at(symPair.firstCell());
              typename lp_type::lp_variable_type x2_ = lp._locs.at(symPair.secondCell());
              auto width_ = cell_loc_trait<LEGALIZE_HORIZONTAL_DIRECTION>::getLen(lp._db, symPair.firstCell); // Two cells are equal in width <- assumption
              lp_type::lp_trait::addConstr(lp._solver, x1_ + x2_ - 2 * (rightSymLoc) <= -width_);
              lp_type::lp_trait::addConstr(lp._solver, x1_ + x2_ + -2 * (leftSymLoc) >= -width_);
            }
            for (IndexType ssIdx = 0; ssIdx < symGroup.numSelfSyms(); ++ssIdx) {
              IndexType cellIdx = symGroup.selfSym(ssIdx);
              auto x_ = lp._locs.at(cellIdx);
              auto width_ = cell_loc_trait<LEGALIZE_HORIZONTAL_DIRECTION>::getLen(lp._db, cellIdx);
              // x + width /2 <= right
              // x + width /2 >= left
              lp_type::lp_trait::addConstr(lp._solver, x_ - (rightSymLoc) <= -width_ / 2);
              lp_type::lp_trait::addConstr(lp._solver, x_ - (leftSymLoc) >= -width_ / 2);
            }
          }
        }
    };

  template<>
    struct _sym_rex_trait<LEGALIZE_VERTICAL_DIRECTION> {
      template<typename lp_type>
        static void addObj(lp_type &lp) {
          // Force they have the same y coordinate
          for (IndexType symGroupIdx = 0; symGroupIdx < lp._db.numSymGroups();
               ++symGroupIdx) {
            const auto &symGroup = lp._db.symGroup(symGroupIdx);
            for (IndexType symPairIdx = 0; symPairIdx < symGroup.numSymPairs();
                 ++symPairIdx) {
              const auto &symPair = symGroup.symPair(symPairIdx);
              IndexType bCellIdx = symPair.firstCell();
              IndexType tCellIdx = symPair.secondCell();
              if (lp._db.cell(bCellIdx).yLoc() > lp._db.cell(tCellIdx).yLoc()) {
                std::swap(tCellIdx, bCellIdx);
              }
              //  + M *( y_t - y_b)
              lp._obj += lp._symVars._largeNum * (lp._locs.at(tCellIdx) - lp._locs.at(bCellIdx));
            }
          }
        }
      template<typename lp_type>
        static void addVars(lp_type &lp) {
          // Same as is horizontal
          _sym_rex_trait<LEGALIZE_HORIZONTAL_DIRECTION>::addVars(lp, "symRex");
        }
      template<typename lp_type>
        static void addConstr(lp_type &lp) {
          // Force they have the same y coordinate
          for (IndexType symGroupIdx = 0; symGroupIdx < lp._db.numSymGroups();
               ++symGroupIdx) {
            const auto &symGroup = lp._db.symGroup(symGroupIdx);
            for (IndexType symPairIdx = 0; symPairIdx < symGroup.numSymPairs();
                 ++symPairIdx) {
              const auto &symPair = symGroup.symPair(symPairIdx);
              IndexType bCellIdx = symPair.firstCell();
              IndexType tCellIdx = symPair.secondCell();
              if (lp._db.cell(bCellIdx).yLoc() > lp._db.cell(tCellIdx).yLoc()) {
                std::swap(tCellIdx, bCellIdx);
              }
              // y_b - y_t <= 0
              lp_type::lp_trait::addConstr(lp._solver,
                                  lp._locs.at(bCellIdx) - lp._locs.at(tCellIdx) <= 0.0);
            }
          }
        }
    };

  template<typename is_hor_type>
    struct _do_not_relax_sym_trait {};

  template<>
    struct _do_not_relax_sym_trait<LEGALIZE_HORIZONTAL_DIRECTION> {
      template<typename lp_type>
        static void addObj(lp_type &) {}
      template<typename lp_type>
        static void addVars(lp_type &lp) {
#ifdef MULTI_SYM_GROUP
          // Symmetric group axis variables
          lp._symVars._symLocs.resize(lp._db.numSymGroups());
          for (IndexType i = 0; i < lp._db.numSymGroups(); ++i) {
            lp._symVars._symLocs.at(i) = lp_type::lp_trait::addVar(lp._solver, "sym" + std::to_string(i));
          }
#else
          lp._symVars._symLocs.resize(1);
          lp._symVars._symLocs[0] = lp_type::lp_trait::addVar(lp._solver, "sym");
#endif
        }
      
      template<typename lp_type>
        static void addConstr(lp_type &lp) {
          // Force them to be symmetric along an axis
          for (IndexType symGrpIdx = 0; symGrpIdx < lp._db.numSymGroups(); ++symGrpIdx) {
            const auto &symGrp = lp._db.symGroup(symGrpIdx);
#ifdef MULTI_SYM_GROUP
            typename lp_type::lp_variable_type symVar = lp._symVars._symLocs.at(symGroupIdx);
#else
            typename lp_type::lp_variable_type symVar = lp._symVars._symLocs.at(0);
#endif

            for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs();
                 ++symPairIdx) {
              const auto &symPair = symGrp.symPair(symPairIdx);
              // x1 + x2 + width =  2 * symAxis
              lp_type::lp_trait::addConstr(
                  lp._solver,
                  lp._locs.at(symPair.firstCell()) + lp._locs.at(symPair.secondCell()) -
                          2 * (symVar) ==
                      -cell_loc_trait<LEGALIZE_HORIZONTAL_DIRECTION>::getLen(lp._db, symPair.firstCell())); // Two cells are equal in width <- assumption
          }
          for (IndexType selfSymIdx = 0; selfSymIdx < symGrp.numSelfSyms();
               ++selfSymIdx) {
            IndexType ssCellIdx = symGrp.selfSym(selfSymIdx);
            // x1 + width + x2 = 2 * symAxis
            lp_type::lp_trait::addConstr(lp._solver,
                                2 * lp._locs.at(ssCellIdx) - 2 * (symVar) ==
                                    - cell_loc_trait<LEGALIZE_HORIZONTAL_DIRECTION>::getLen(lp._db, ssCellIdx));
          }
        }
      }
    };

  template<>
    struct _do_not_relax_sym_trait<LEGALIZE_VERTICAL_DIRECTION> {
      template<typename lp_type>
        static void addObj(lp_type &) {}
      template<typename lp_type>
        static void addVars(lp_type &lp) {
#ifdef MULTI_SYM_GROUP
          // Symmetric group axis variables
          lp._symVars._symLocs.resize(lp._db.numSymGroups());
          for (IndexType i = 0; i < lp._db.numSymGroups(); ++i) {
            lp._symVars._symLocs.at(i) = lp_type::lp_trait::addVar(lp._solver, sym + std::to_string(i));
          }
#else
          lp._symVars._symLocs.resize(1);
          lp._symVars._symLocs[0] = lp_type::lp_trait::addVar(lp._solver, "sym");
#endif
        }
      template<typename lp_type>
        static void addConstr(lp_type &lp) {
          // Force they have the same y coordinate
          for (IndexType symGroupIdx = 0; symGroupIdx < lp._db.numSymGroups();
               ++symGroupIdx) {
            const auto &symGroup = lp._db.symGroup(symGroupIdx);
            for (IndexType symPairIdx = 0; symPairIdx < symGroup.numSymPairs();
                 ++symPairIdx) {
              const auto &symPair = symGroup.symPair(symPairIdx);
              // y_i = y_j
              lp_type::lp_trait::addConstr(lp._solver, lp._locs.at(symPair.firstCell()) -
                                                   lp._locs.at(symPair.secondCell()) ==
                                               0.0);
            }
          }
        }
    };

  template<typename relax_sym_type, typename is_hor_type>
    struct _select_sym_base_trait {};
  template<typename is_hor_type>
    struct _select_sym_base_trait<DO_NOT_RELAX_SYM_CONSTR, is_hor_type> {
      typedef _do_not_relax_sym_trait<is_hor_type> base_type;
    };
  template<typename is_hor_type>
    struct _select_sym_base_trait<RELAX_SYM_CONSTR, is_hor_type> {
      typedef _sym_rex_trait<is_hor_type> base_type;
    };

  template<typename is_hor_type, typename relax_sym_type>
    struct sym_trait {
      typedef typename _select_sym_base_trait<relax_sym_type, is_hor_type>::base_type base_type;
      template<typename lp_type>
        static void addObj(lp_type &lp) { base_type::addObj(lp); }
      template<typename lp_type>
        static void addVars(lp_type &lp) { base_type::addVars(lp); }
      template<typename lp_type>
        static void addConstr(lp_type &lp) { base_type::addConstr(lp); }
    };

  struct spacing_trait {
    template<typename lp_type>
      static void addObj(lp_type &) {}
    template<typename lp_type>
      static void addVars(lp_type &) {}
    template<typename lp_type>
      static void addConstr(lp_type &lp) {
        for (const auto &edge : lp._constrs.edges()) {
          IndexType sourceIdx = edge.source();
          IndexType targetIdx = edge.target();
          if (sourceIdx == targetIdx) {
            continue;
          }
          const IndexType numCellAndRects = lp._db.numCells()  + lp._db.numWellRects();
          if (sourceIdx >= numCellAndRects or targetIdx >= numCellAndRects) {
            // the s, t constraints
            continue;
          }
          LocType spacing = cell_loc_trait<typename lp_type::is_hor_type>::getSpacing(lp._db, sourceIdx, targetIdx);
          LocType cellLen = cell_loc_trait<typename lp_type::is_hor_type>::getLen(lp._db, sourceIdx);
          if (lp._extraSpacing > 0) {
            spacing +=  lp._extraSpacing;
          }
          // Add the constraint
          // x_i + w_i + spacing <= x_j
          lp_type::lp_trait::addConstr(lp._solver, lp._locs.at(sourceIdx) - lp._locs.at(targetIdx) <=
                                           -cellLen - spacing);
        }
      }
  };

  // Ensure the splitted rectangles and cells inside are moving together
  struct well_align_trait {
    template<typename lp_type>
      static void addObj(lp_type &) {}
    template<typename lp_type>
      static void addVars(lp_type &) {}
    template<typename lp_type>
      static void addConstr(lp_type &lp) {
        using cell_loc_type = cell_loc_trait<typename lp_type::is_hor_type>;
        const Database &db = lp._db;
        IndexType locIdx = db.numCells(); // Loc starting from cells then wells
        for (IndexType wellIdx = 0; wellIdx < db.vWells().size(); ++wellIdx) {
          const auto &well = db.well(wellIdx);
          Assert(well.rects().size() > 0);
          LocType firstRectLoc = cell_loc_type::getWellRectLL(well, 0);
          IndexType firstWellIdx = locIdx;
          ++locIdx;
          // Make sure the splitted well rectangles are moving together
          for (IndexType rectIdx = 1; rectIdx < well.rects().size(); ++rectIdx) {
            LocType currentRectLoc = cell_loc_type::getWellRectLL(well, rectIdx);
            // Add constraint
            // x_first - x_cur = their original location difference
            lp_type::lp_trait::addConstr(lp._solver, lp._locs.at(firstWellIdx) - lp._locs.at(locIdx)
                == firstRectLoc - currentRectLoc);
            ++locIdx;
          }
          // Make sure the cell inside are moving together
          for (IndexType cellIdx : well.sCellIds()) {
            LocType cellLoc = cell_loc_type::getCellLL(db, cellIdx);
            // Add constraint
            // x_first - x_cell = their original location difference
            lp_type::lp_trait::addConstr(lp._solver, lp._locs.at(firstWellIdx) - lp._locs.at(cellIdx)
                == firstRectLoc - cellLoc);
          }
        }
      }
  };


  } // namespace _lp_legalize_details
} // namespace lp_legalize

namespace lp_legalize {


  /// @brief The base class for LP-based legalizer.
  /// It implememts the essential shared functions for different legalization tasks.
  /// @tparam first: LEGALIZE_HORIZONTAL_DIRECTION-> legalize horizontal x coordinates. LEGALIZE_VERTICAL_DIRECTION-> vertical y coordinates.
  template<typename is_hor=lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, typename relax_sym=lp_legalize::DO_NOT_RELAX_SYM_CONSTR>
  class _LpLegalizeBase {
      public:
        typedef is_hor is_hor_type;
        typedef relax_sym relax_sym_type;
    protected:
    typedef ::klib::lp::LpModel lp_solver_type;
    typedef ::klib::lp::LpTrait lp_trait;
    typedef lp_trait::variable_type lp_variable_type;
    typedef lp_trait::expr_type lp_expr_type;

    typedef _lp_legalize_details::hpwl_trait hpwl_trait;
    friend hpwl_trait;
    typedef _lp_legalize_details::area_trait area_trait;
    friend area_trait;
    typedef _lp_legalize_details::sym_trait<is_hor_type, relax_sym_type> sym_trait;
    friend sym_trait;
    friend typename  sym_trait::base_type;
    typedef _lp_legalize_details::sym_variable_type<lp_variable_type, relax_sym_type> sym_variable_type;
    typedef _lp_legalize_details::spacing_trait spacing_trait;
    friend spacing_trait;
    typedef _lp_legalize_details::well_align_trait well_align_trait;
    friend well_align_trait;

    public:
    explicit _LpLegalizeBase(Database &db, Constraints &constrs) : _db(db), _constrs(constrs) {
      addLocVars();
    }

    BoolType solve() {
      addVariables();
      addConstraints();
      addObjective();
      return solveLpModel();
    }

    void exportSolution() {
      using cell_loc_type = _lp_legalize_details::cell_loc_trait<is_hor_type>;
      for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx) {
        auto var = lp_trait::solution(_solver, _locs.at(cellIdx));
        // convert to cell original location
        cell_loc_type::setCellLL(_db, cellIdx, ::klib::autoRound<LocType>(var));
      }
      IndexType locIdx = _db.numCells(); // Loc starting from cells then wells
      for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
        auto &well = _db.well(wellIdx);
        Assert(well.rects().size() > 0);
        LocType firstRectLoc = cell_loc_type::getWellRectLL(well, 0);
        auto var = lp_trait::solution(_solver, _locs.at(locIdx));
        LocType offset = ::klib::autoRound<LocType>(var) - firstRectLoc;
        // Move the whole well together
        cell_loc_type::moveWell(well, offset);
        locIdx += _db.numRectInWell(wellIdx);
      }
    }

    void setExtraSpacing(LocType spacing) {
      _extraSpacing = spacing;
    }

    protected:
    void addLocVars() {
      // _locs are x or y, depending on is_hor_type
      _locs.resize(_db.numCells() + _db.numWellRects());
      for (IndexType i = 0; i < _db.numCells() + _db.numWellRects(); ++i) {
        _locs.at(i) = lp_trait::addVar(_solver, "loc"+std::to_string(i));
      }
    }
    virtual void addVariables() { AssertMsg(false, "%s: Should not call _LpLegalizeBase implementation. ", __FUNCTION__); };
    virtual void addConstraints() { AssertMsg(false, "%s: Should not call _LpLegalizeBase implementation. ", __FUNCTION__); };
    virtual void addObjective() { AssertMsg(false, "%s: Should not call _LpLegalizeBase implementation. ", __FUNCTION__); };
    BoolType solveLpModel() { 
      lp_trait::setNumThreads(_solver, _db.parameters().numThreads());
      lp_trait::setObjectiveMinimize(_solver);
      lp_trait::setObjective(_solver, _obj);
      lp_trait::solve(_solver);
      if (lp_trait::isUnbounded(_solver)) {
        ERR("LP legalization solver: LP unbounded. \n");
        return false;
      } else if (lp_trait::isOptimal(_solver)) {
        for (IndexType idx = 0; idx < _db.numCells(); ++idx) {
          auto var = lp_trait::solution(_solver, _locs.at(idx));
          if (::klib::autoRound<LocType>(var) < 0) { // Gurobi sometimes report OPTIMAL while having everything messed up
            ERR("LP legalization solver: LP infeasible. \n");
            assert(false);
            return false;
          }
        }
        INF("LP legalization solver: LP optimal \n");
        return true;
      } else if (lp_trait::isInfeasible(_solver)) {
        ERR("LP legalization solver: LP infeasible. \n");
        return false;
      } else if (lp_trait::isSuboptimal(_solver)) {
        ERR("LP legalization solver: LP suboptimal.");
        return false;
      } else {
        ERR("LP legalization solver: Unknown LP status  %s\n", lp_trait::statusStr(_solver).c_str());
        return false;
      }
    }

    protected:
    /* Configurations - Inputs */
    Database &_db;            ///< The database for the Ideaplace
    Constraints &_constrs; ///< The topology constraints
    /* Optimization supporting variables */
    lp_solver_type _solver; ///<  LP sovler
    lp_expr_type _obj;      ///< The objective function of the LP model
    std::vector<lp_variable_type> _locs; ///< The location variables of the LP model
    LocType _extraSpacing = -1; ///< Apply extra spacing
  };


  /// @brief The wire-length-driven LP-based legalizer.
  /// @tparam first: LEGALIZE_HORIZONTAL_DIRECTION-> legalize horizontal x coordinates. LEGALIZE_VERTICAL_DIRECTION-> vertical y coordinates.
  template<typename is_hor=lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, typename relax_sym=lp_legalize::DO_NOT_RELAX_SYM_CONSTR>
    class LpLegalizeWirelength : public _LpLegalizeBase<is_hor, relax_sym> {
      public:
        typedef is_hor is_hor_type;
        typedef relax_sym relax_sym_type;
      protected:
      typedef _LpLegalizeBase<is_hor_type, relax_sym_type> base_type;
      typedef typename base_type::lp_solver_type lp_solver_type;
      typedef typename base_type::lp_trait lp_trait;
      typedef typename base_type::lp_variable_type lp_variable_type;
      typedef typename base_type::lp_expr_type lp_expr_type;

      typedef typename base_type::hpwl_trait hpwl_trait;
      friend hpwl_trait;
      typedef typename base_type::area_trait area_trait;
      friend area_trait;
      typedef typename base_type::sym_trait sym_trait;
      friend sym_trait;
      friend typename  sym_trait::base_type;
      typedef typename base_type::sym_variable_type sym_variable_type;
      typedef typename base_type::spacing_trait spacing_trait;
      friend spacing_trait;
      typedef typename base_type::well_align_trait well_align_trait;
      friend well_align_trait;

      public:
      LpLegalizeWirelength(Database &db, Constraints &constrs) 
        : base_type(db, constrs) {}
      void setBoundary(LocType boundary) {
        _fixedBoundary = static_cast<RealType>(boundary);
      }

      protected:
      virtual void addVariables() override {
        hpwl_trait::addVars(*this); // Wirelength variable
        area_trait::addVars(*this); // Boundary
        sym_trait::addVars(*this); // Sym axis variable
        spacing_trait::addVars(*this); // Spacing
        well_align_trait::addVars(*this);
      }

      virtual void addConstraints() override {
        hpwl_trait::addConstr(*this);
        area_trait::addConstr(*this);
        if (_fixedBoundary > 0) {
          area_trait::fixBoundary(*this, _fixedBoundary);
        }
        sym_trait::addConstr(*this);
        spacing_trait::addConstr(*this);
        well_align_trait::addConstr(*this);
      }

      virtual void addObjective() override {
        hpwl_trait::addObj(*this);
        sym_trait::addObj(*this);
        spacing_trait::addObj(*this);
        well_align_trait::addObj(*this);
      }

      protected:
      std::vector<lp_variable_type>
          _wlL; ///< The left wirelength variables of the ILP model
      std::vector<lp_variable_type>
          _wlR;              ///< The right wirelength variables of the ILP model
      lp_variable_type _boundary; ///< The boundary of the design
      sym_variable_type _symVars;
      RealType _fixedBoundary = -1.0; ///< The fixed boundary constraint
    };

  /// @brief The area-length-driven LP-based legalizer.
  /// @tparam first: LEGALIZE_HORIZONTAL_DIRECTION-> legalize horizontal x coordinates. LEGALIZE_VERTICAL_DIRECTION-> vertical y coordinates.
  template<typename is_hor=lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, typename relax_sym=lp_legalize::DO_NOT_RELAX_SYM_CONSTR>
    class LpLegalizeArea : public _LpLegalizeBase<is_hor, relax_sym> {
      public:
        typedef is_hor is_hor_type;
        typedef relax_sym relax_sym_type;
      protected:
      typedef _LpLegalizeBase<is_hor_type, relax_sym_type> base_type;
      typedef typename base_type::lp_solver_type lp_solver_type;
      typedef typename base_type::lp_trait lp_trait;
      typedef typename base_type::lp_variable_type lp_variable_type;
      typedef typename base_type::lp_expr_type lp_expr_type;

      typedef typename base_type::hpwl_trait hpwl_trait;
      friend hpwl_trait;
      typedef typename base_type::area_trait area_trait;
      friend area_trait;
      typedef typename base_type::sym_trait sym_trait;
      friend sym_trait;
      friend typename  sym_trait::base_type;
      typedef typename base_type::sym_variable_type sym_variable_type;
      typedef typename base_type::spacing_trait spacing_trait;
      friend spacing_trait;
      typedef typename base_type::well_align_trait well_align_trait;
      friend well_align_trait;

      public:
      LpLegalizeArea(Database &db, Constraints &constrs) 
        : base_type(db, constrs) {}



      protected:
      virtual void addVariables() override {
        area_trait::addVars(*this); // Boundary
        sym_trait::addVars(*this); // Sym axis variable
        spacing_trait::addVars(*this); // Spacing
        well_align_trait::addVars(*this);
      }

      virtual void addConstraints() override {
        area_trait::addConstr(*this);
        sym_trait::addConstr(*this);
        spacing_trait::addConstr(*this);
        well_align_trait::addConstr(*this);
      }

      virtual void addObjective() override {
        area_trait::addObj(*this);
        sym_trait::addObj(*this);
        spacing_trait::addObj(*this);
        well_align_trait::addObj(*this);
      }

      protected:
      lp_variable_type _boundary; ///< The boundary of the design
      sym_variable_type _symVars;
    };

} // namespace lp_legalize

PROJECT_NAMESPACE_END
