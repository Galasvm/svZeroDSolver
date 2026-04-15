// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file autoregulation.h
 * @brief model::autoregulation header
 */
#ifndef SVZERODSOLVER_MODEL_AUTOREGULATION_HPP_
#define SVZERODSOLVER_MODEL_AUTOREGULATION_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Autoregulation boundary condition: 4 resistors in series with
 * shear + myogenic + metabolic control.
 *
 * Local unknown ordering:
 *   y^e = [ Pin, Qin, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS ]^T
 *
 * Governing equations (in this order):
 *  (0) Pin - Pd - Qin*(R1(Ashear) + R2(Amyo) + R3(Ameta) + R4) = 0
 *  (1) T   - Pavg * (Kar2 / R2(Amyo))^0.25 = 0
 *      where Pavg = (P1 + P2)/2,
 *            P1 = Pin - Qin*R1,
 *            P2 = P1  - Qin*R2  = Pin - Qin*(R1+R2),
 *            => Pavg = Pin - Qin*(R1 + 0.5*R2)
 *  (2) WSS - Qin / r1^3 = 0, with r1 = (Kar1 / R1(Ashear))^0.25
 *      => WSS = Qin * (R1/Kar1)^0.75
 *
 *  (3) dAshear/dt - Gshear*xshear = 0
 *  (4) dAmyo/dt   - Gmyo*xmyo     = 0
 *  (5) dAmeta/dt  - Gmeta*xmeta   = 0
 *
 *  (6) TAUshear*dxshear/dt + xshear - WSS/WSSt + 1 = 0
 *  (7) TAUmyo*dxmyo/dt     + xmyo   - T/Tt     + 1 = 0
 *  (8) TAUmeta*dxmeta/dt   + xmeta  - Qin/Qt   + 1 = 0
 *
 * Parameters:
 *  - R        baseline total resistance (used to define R1_0,R2_0,R3_0,R4)
 *  - Qt       target flow
 *  - Pt       target pressure (used to define Tt via baseline P1_0,P2_0)
 *  - Pd       distal pressure (for eqn 0)
 *  - Gshear, TAUshear
 *  - Gmyo,   TAUmyo
 *  - Gmeta,  TAUmeta
 *
 * One-time constants computed from parameters:
 *  - R1L,R1U,R2L,R2U,R3L,R3U,R4
 *  - Kar1 = 0.02^4 * R1_0
 *  - Kar2 = 0.01^4 * R2_0
 *  - WSSt = Qt / 0.02^3
 *  - Tt   = (P1_0+P2_0)/2 * (Kar2/R2_0)^0.25
 */
class Autoregulation : public Block {
 public:
  Autoregulation(int id, Model *model)
      : Block(id, model,
              BlockType::autoregulation,
              BlockClass::boundary_condition,
              {{"R", InputParameter()},
               {"Qt", InputParameter()},
               {"Pt", InputParameter()},
               {"Gshear", InputParameter()},
               {"taushear", InputParameter()},
               {"Gmyo", InputParameter()},
               {"taumyo", InputParameter()},
               {"Gmeta", InputParameter()},
               {"taumeta", InputParameter()},
               {"Pd", InputParameter()},
               {"lower_frac", InputParameter(true, false, true, 0.70)},
               {"upper_frac", InputParameter(true, false, true, 1.30)}}) {}

  void setup_dofs(DOFHandler &dofhandler);

  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  // F nnz ~12, E nnz ~6, dC_dy nnz ~10
  TripletsContributions num_triplets{12, 6, 10};

 private:
  bool initialized_ = false;

  // stored one-time constants
  double R1L_ = 0.0, R1U_ = 0.0;
  double R2L_ = 0.0, R2U_ = 0.0;
  double R3L_ = 0.0, R3U_ = 0.0;
  double R4_  = 0.0;

  double Kar1_ = 0.0;
  double Kar2_ = 0.0;
  double WSSt_ = 0.0;
  double Tt_   = 0.0;
};

#endif  // SVZERODSOLVER_MODEL_AUTOREGULATION_HPP_