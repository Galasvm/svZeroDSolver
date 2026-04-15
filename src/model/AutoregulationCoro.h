// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file AutoregulationCoro.h
 * @brief model::AutoregulationCoro header
 */
#ifndef SVZERODSOLVER_MODEL_AUTOREGULATIONCORO_HPP_
#define SVZERODSOLVER_MODEL_AUTOREGULATIONCORO_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Open-loop coronary boundary condition with autoregulated microvascular
 * resistance Ra2 (shear + myogenic + metabolic control).
 *
 * Extends the standard open-loop coronary BC (Kim et al.) by replacing the
 * fixed microvascular resistance Ra2 with a regulated Rtot composed of four
 * resistors in series (R1, R2, R3, R4), each controlled by sigmoid activation
 * functions driven by shear stress, myogenic tone, and metabolic demand.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_{a1}$, *-*] (3,0)
 * to [R, l=$R_{tot}(t)$, -] (5,0)
 * node[anchor=south]{$P_{cim}$}
 * to [R, l=$R_{v1}$, *-*] (7,0)
 * node[anchor=south]{$P_{v}$}
 * (5,0) to [C, l=$C_{im}$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a$, -*] (3,-1.5)
 * node[left]{$P_a=0$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * Coronary equations (same structure as OpenLoopCoronaryBC, with Ra2 → Rtot):
 *  (0) Cim*Rv*Qin - Vim - Cim*Rv*dVim/dt - Ca*Cim*Rv*dPin/dt + Ra*Ca*Cim*Rv*dQin/dt
 *      + Cim*(-Pim + Pv + Pim_0 - P_Cim_0) = 0
 *  (1) Cim*Rv*Pin - Cim*Rv*Ra*Qin - Rv*Vim - Cim*Rv*Pim_offset
 *      + Rtot*(-Vim + Cim*(Pv - Pim_offset) - Cim*Rv*dVim/dt) = 0
 *
 * Autoregulation equations (same as Autoregulation block):
 *  (2) T   - Pavg*(Kar2/R2)^0.25 = 0,  Pavg = Pin - Qin*(R1 + 0.5*R2)
 *  (3) WSS - Qin*(R1/Kar1)^0.75  = 0
 *  (4) dAshear/dt + Gshear*xshear = 0
 *  (5) dAmyo/dt   - Gmyo*xmyo    = 0
 *  (6) dAmeta/dt  - Gmeta*xmeta  = 0
 *  (7) TAUshear*dxshear/dt + xshear - WSS/WSSt + 1 = 0
 *  (8) TAUmyo*dxmyo/dt     + xmyo   - T/Tt     + 1 = 0
 *  (9) TAUmeta*dxmeta/dt   + xmeta  - Qin/Qt   + 1 = 0
 *
 * Local unknowns:
 *   y^e = [Pin, Qin, Vim, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS]
 *
 * ### Parameters
 *
 * * `0` Ra1:      Proximal small artery resistance (fixed)
 * * `1` Rv1:      Venous resistance (fixed)
 * * `2` Ca:       Small artery capacitance
 * * `3` Cc:       Intramyocardial capacitance (Cim)
 * * `4` Pim:      Intramyocardial pressure (time-varying, paired with t)
 * * `5` P_v:      Venous pressure
 * * `6` Ra2:      Baseline regulated microvascular resistance
 * * `7` Qt:       Target flow
 * * `8` Pt:       Target pressure at Pin
 * * `9` Gshear,  `10` taushear
 * * `11` Gmyo,   `12` taumyo
 * * `13` Gmeta,  `14` taumeta
 * * `15` lower_frac (optional, default 0.70) lower bound fraction for R1_0,R2_0,R3_0
 * * `16` upper_frac (optional, default 1.30) upper bound fraction for R1_0,R2_0,R3_0
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
               {"lower_frac", InputParameter(true, false, true, 0.70)},
               {"upper_frac", InputParameter(true, false, true, 1.30)}}) {}

  void setup_dofs(DOFHandler &dofhandler);

  void setup_initial_state_dependent_params(State initial_state,
                                            std::vector<double> &parameters);

  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  void update_time(SparseSystem &system, std::vector<double> &parameters);

  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  // F=16, E=9, D=11
  TripletsContributions num_triplets{16, 9, 11};

 private:
  bool initialized_ = false;

  double P_Cim_0 = 0.0;  ///< Pressure proximal to Cim at initial state
  double Pim_0   = 0.0;  ///< Pim at initial state

  double R1L_ = 0.0, R1U_ = 0.0;
  double R2L_ = 0.0, R2U_ = 0.0;
  double R3L_ = 0.0, R3U_ = 0.0;
  double R4_  = 0.0;

  double Kar1_ = 0.0;
  double Kar2_ = 0.0;
  double WSSt_ = 0.0;
  double Tt_   = 0.0;
};

#endif  // SVZERODSOLVER_MODEL_AUTOREGULATIONCORONARYBC_HPP_
