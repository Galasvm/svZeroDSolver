// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_AUTOREGULATIONCORO_HPP_
#define SVZERODSOLVER_MODEL_AUTOREGULATIONCORO_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Open-loop coronary BC with autoregulated proximal (Ra) and
 * microvascular (Ra2) resistances.
 *
 * Circuit:
 *   Pin -> Ra1_static(30%) -> Ra1_shear(16%,As) -> Ra1_myo(54%,Am)
 *       -> [Ca] -> Ra2(100%,Amet) -> [Cim,Pim] -> Rv -> Pv
 *
 * Ra (= Ra1_static + Ra1_shear + Ra1_myo) sits before Ca.
 * Ra2 (fully metabolic) sits after Ca.
 *
 * Local unknowns:
 *   y^e = [Pin, Qin, Vim, Ashear, Amyo, Ameta, xshear, xmyo, xmeta,
 *          T, WSS, Pa, q_micro]
 *
 * Equations:
 *  (0-1)  Coronary hydraulics (Kim et al., Ra_eff = Ra1_static+Ra1_shear+Ra1_myo)
 *  (2)    T   - Pavg*(Kar_myo/Ra1_myo)^0.25 = 0
 *             Pavg = Pa + 0.5*Ra1_myo*Qin  (average across Ra1_myo)
 *  (3)    WSS - Qin*(Ra1_shear/Kar1_shear)^0.75 = 0
 *  (4)    dAshear/dt + Gshear*xshear = 0
 *  (5)    dAmyo/dt   - Gmyo*xmyo     = 0
 *  (6)    dAmeta/dt  - Gmeta*xmeta   = 0
 *  (7)    TAUshear*dxshear/dt + xshear - WSS/WSSt  + 1 = 0
 *  (8)    TAUmyo*dxmyo/dt     + xmyo   - T/Tt      + 1 = 0
 *  (9)    TAUmeta*dxmeta/dt   + xmeta  - q_micro/Qt + 1 = 0
 *  (10)   Pa - Pin + Ra1_static*Qin + Ra1_shear*Qin + Ra1_myo*Qin = 0
 *  (11)   q_micro - Qin + Ca*dPa/dt = 0
 */
class AutoregulationCoro : public Block {
 public:
  AutoregulationCoro(int id, Model *model)
      : Block(id, model,
              BlockType::autoregulation_coro,
              BlockClass::boundary_condition,
              {{"Ra1", InputParameter()},
               {"Rv1", InputParameter()},
               {"Ca", InputParameter()},
               {"Cc", InputParameter()},
               {"t", InputParameter(false, true)},
               {"Pim", InputParameter(false, true)},
               {"P_v", InputParameter()},
               {"Ra2", InputParameter()},
               {"Qt", InputParameter()},
               {"Pt", InputParameter()},
               {"Gshear", InputParameter()},
               {"taushear", InputParameter()},
               {"Gmyo", InputParameter()},
               {"taumyo", InputParameter()},
               {"Gmeta", InputParameter()},
               {"taumeta", InputParameter()},
               {"lower_frac", InputParameter(true, false, true, 0.03)},
               {"upper_frac", InputParameter(true, false, true, 1.97)}}) {}

  void setup_dofs(DOFHandler &dofhandler);

  void setup_initial_state_dependent_params(State initial_state,
                                            std::vector<double> &parameters);

  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  // F=21, E=10, D=15
  TripletsContributions num_triplets{21, 10, 15};

 private:
  bool initialized_ = false;

  double P_Cim_0 = 0.0;
  double Pim_0   = 0.0;

  double Ra1_static_ = 0.0;          ///< Fixed 30% of Ra
  double Ra1SL_ = 0.0, Ra1SU_ = 0.0; ///< Ra1_shear sigmoid bounds (16% of Ra)
  double Ra1ML_ = 0.0, Ra1MU_ = 0.0; ///< Ra1_myo  sigmoid bounds (54% of Ra)
  double Ra2L_  = 0.0, Ra2U_  = 0.0; ///< Ra2 metabolic sigmoid bounds (100% of Ra2)

  double Kar1_shear_ = 0.0; ///< Geometric constant for shear WSS (Ra1_shear)
  double Kar_myo_    = 0.0; ///< Geometric constant for myogenic tension (Ra1_myo)
  double WSSt_ = 0.0;
  double Tt_   = 0.0;
};

#endif  // SVZERODSOLVER_MODEL_AUTOREGULATIONCORONARYBC_HPP_
