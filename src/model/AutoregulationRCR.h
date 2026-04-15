// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file AutoregulationRCR.h
 * @brief model::AutoregulationRCR header
 */
#ifndef SVZERODSOLVER_MODEL_AUTOREGULATIONRCR_HPP_
#define SVZERODSOLVER_MODEL_AUTOREGULATIONRCR_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief AutoregulationRCR boundary condition: RCR Windkessel with a regulated
 * distal resistance (4 resistors in series, shear + myogenic + metabolic
 * control).
 *
 * A fixed proximal resistance Rp and capacitance C sit upstream of the
 * regulated vascular segment (R1+R2+R3+R4).  Pc is the pressure at the
 * capacitor node.
 *
 * Local unknown ordering:
 *   y^e = [ Pin, Qin, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS, Pc ]
 *
 * Governing equations:
 *  (0) Pin - Pc - Rp*Qin = 0
 *  (1) Rtot*(Qin - C*dPc/dt) - Pc + Pd = 0
 *      where Rtot = R1(Ashear) + R2(Amyo) + R3(Ameta) + R4
 *  (2) T   - Pavg*(Kar2/R2(Amyo))^0.25 = 0
 *      where Pavg = Pc - Qin*(R1 + 0.5*R2)
 *  (3) WSS - Qin*(R1/Kar1)^0.75 = 0
 *  (4) dAshear/dt - Gshear*xshear = 0
 *  (5) dAmyo/dt   - Gmyo*xmyo     = 0
 *  (6) dAmeta/dt  - Gmeta*xmeta   = 0
 *  (7) TAUshear*dxshear/dt + xshear - WSS/WSSt + 1 = 0
 *  (8) TAUmyo*dxmyo/dt     + xmyo   - T/Tt     + 1 = 0
 *  (9) TAUmeta*dxmeta/dt   + xmeta  - Qin/Qt   + 1 = 0
 *
 * Parameters:
 *  - Rd       distal resistance (regulated)
 *  - Qt       target flow
 *  - Pt       target pressure at Pin
 *  - Pd       distal pressure
 *  - Gshear, taushear
 *  - Gmyo,   taumyo
 *  - Gmeta,  taumeta
 *  - Rp       proximal (fixed) resistance
 *  - C        capacitance
 */
class AutoregulationRCR : public Block {
 public:
  AutoregulationRCR(int id, Model *model)
      : Block(id, model,
              BlockType::autoregulation_rcr,
              BlockClass::boundary_condition,
              {{"Rd", InputParameter()},
               {"Qt", InputParameter()},
               {"Pt", InputParameter()},
               {"Gshear", InputParameter()},
               {"taushear", InputParameter()},
               {"Gmyo", InputParameter()},
               {"taumyo", InputParameter()},
               {"Gmeta", InputParameter()},
               {"taumeta", InputParameter()},
               {"Pd", InputParameter()},
               {"Rp", InputParameter()},
               {"C", InputParameter()},
               {"lower_frac", InputParameter(true, false, true, 0.70)},
               {"upper_frac", InputParameter(true, false, true, 1.30)}}) {}

  void setup_dofs(DOFHandler &dofhandler);

  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  // F=15, E=6, D=10
  TripletsContributions num_triplets{15, 6, 10};

 private:
  bool initialized_ = false;

  double R1L_ = 0.0, R1U_ = 0.0;
  double R2L_ = 0.0, R2U_ = 0.0;
  double R3L_ = 0.0, R3U_ = 0.0;
  double R4_  = 0.0;

  double Kar1_ = 0.0;
  double Kar2_ = 0.0;
  double WSSt_ = 0.0;
  double Tt_   = 0.0;
};

#endif  // SVZERODSOLVER_MODEL_AUTOREGULATIONRCR_HPP_
