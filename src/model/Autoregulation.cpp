// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Autoregulation.h"

void Autoregulation::setup_dofs(DOFHandler &dofhandler) {
  // Total local dofs = 9
  Block::setup_dofs_(dofhandler, 9,
                     {"Ashear", "Amyo", "Ameta", "xshear", "xmyo", "xmeta", "T", "WSS"});
}

void Autoregulation::update_constant(SparseSystem &system,
                                     std::vector<double> &parameters) {
  // Parameters (by global_param_ids)
  const double R        = parameters[global_param_ids[0]];
  const double Qt       = parameters[global_param_ids[1]];
  const double Pt       = parameters[global_param_ids[2]];
  const double Gshear   = parameters[global_param_ids[3]];
  const double TAUshear = parameters[global_param_ids[4]];
  const double Gmyo     = parameters[global_param_ids[5]];
  const double TAUmyo   = parameters[global_param_ids[6]];
  const double Gmeta    = parameters[global_param_ids[7]];
  const double TAUmeta  = parameters[global_param_ids[8]];
  (void)parameters[global_param_ids[9]];  // Pd used in update_solution

  // --------------------
  // One-time precompute
  // --------------------
  if (!initialized_) {
    const double R1_0 = 0.10 * R;
    const double R2_0 = 0.35 * R;
    const double R3_0 = 0.40 * R;

    R4_  = 0.15 * R;

    R1L_ = 0.70 * R1_0;  R1U_ = 1.30 * R1_0;
    R2L_ = 0.70 * R2_0;  R2U_ = 1.30 * R2_0;
    R3L_ = 0.70 * R3_0;  R3U_ = 1.30 * R3_0;

    Kar1_ = std::pow(0.02, 4) * R1_0;
    Kar2_ = std::pow(0.01, 4) * R2_0;

    WSSt_ = Qt / std::pow(0.02, 3);

    // Baseline pressures for target tension
    const double P1_0 = Pt - Qt * R1_0;
    const double P2_0 = P1_0 - Qt * R2_0;
    const double Pavg0 = 0.5 * (P1_0 + P2_0);

    // Tt = Pavg0 * (Kar2 / R2_0)^0.25
    Tt_ = Pavg0 * std::pow(Kar2_ / R2_0, 0.25);

    initialized_ = true;
  }

  // --------------------
  // Fill constant F entries
  // --------------------
  // Local y: [Pin, Qin, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS]

  // Eqn (0): Pin - Pd - Qin*(...) = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;  // +Pin

  // Eqn (1): T - (...) = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[8]) = 1.0;  // +T

  // Eqn (2): WSS - (...) = 0
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[9]) = 1.0;  // +WSS

  // Eqn (3): dAshear/dt + Gshear*xshear = 0
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[5]) = Gshear;

  // Eqn (4): dAmyo/dt - Gmyo*xmyo = 0
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[6]) = -Gmyo;

  // Eqn (5): dAmeta/dt - Gmeta*xmeta = 0
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[7]) = -Gmeta;

  // Eqn (6): TAUshear*dxshear/dt + xshear - WSS/WSSt + 1 = 0
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[5]) = 1.0;          // +xshear
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[9]) = -1.0 / WSSt_;  // -WSS/WSSt

  // Eqn (7): TAUmyo*dxmyo/dt + xmyo - T/Tt + 1 = 0
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[6]) = 1.0;        // +xmyo
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[8]) = -1.0 / Tt_;  // -T/Tt

  // Eqn (8): TAUmeta*dxmeta/dt + xmeta - Qin/Qt + 1 = 0
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[7]) = 1.0;       // +xmeta
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[1]) = -1.0 / Qt; // -Qin/Qt

  // --------------------
  // Fill constant E entries (time derivatives)
  // --------------------
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[2]) = 1.0;       // dAshear/dt
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[3]) = 1.0;       // dAmyo/dt
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[4]) = 1.0;       // dAmeta/dt

  system.E.coeffRef(global_eqn_ids[6], global_var_ids[5]) = TAUshear;  // TAUshear*dxshear/dt
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[6]) = TAUmyo;    // TAUmyo*dxmyo/dt
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[7]) = TAUmeta;   // TAUmeta*dxmeta/dt

  // --------------------
  // Constant C offsets for eqns (6..8): +1
  // --------------------
  system.C.coeffRef(global_eqn_ids[6]) = 1.0;
  system.C.coeffRef(global_eqn_ids[7]) = 1.0;
  system.C.coeffRef(global_eqn_ids[8]) = 1.0;
}

void Autoregulation::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {

  (void)dy;

  // States
  const double p_in  = y[global_var_ids[0]];
  const double q_in  = y[global_var_ids[1]];
  const double As    = y[global_var_ids[2]];
  const double Am    = y[global_var_ids[3]];
  const double Amet  = y[global_var_ids[4]];

  // Parameter
  const double Pd    = parameters[global_param_ids[9]];

  // exp(A) (if you see overflow, switch to a stable sigmoid implementation)
  const double eS   = std::exp(As);
  const double eM   = std::exp(Am);
  const double eMet = std::exp(Amet);

  // Effective resistors
  const double R1 = (R1L_ + R1U_ * eS)   / (1.0 + eS);
  const double R2 = (R2L_ + R2U_ * eM)   / (1.0 + eM);
  const double R3 = (R3L_ + R3U_ * eMet) / (1.0 + eMet);

  const double Rtot = R1 + R2 + R3 + R4_;

  // dRi/dAi
  const double dR1_dAs =
      (R1U_ - R1L_) * eS / ((1.0 + eS) * (1.0 + eS));
  const double dR2_dAm =
      (R2U_ - R2L_) * eM / ((1.0 + eM) * (1.0 + eM));
  const double dR3_dAmet =
      (R3U_ - R3L_) * eMet / ((1.0 + eMet) * (1.0 + eMet));

  // --------------------
  // Nonlinear C entries (eqns 0..2)
  // --------------------

  // Eqn (0): Pin - Pd - Qin*Rtot = 0
  system.C(global_eqn_ids[0]) = -Pd - q_in * Rtot;

  // Eqn (1): T - Pavg*(Kar2/R2)^0.25 = 0
  // Pavg = Pin - Qin*(R1 + 0.5*R2)
  const double Pavg = p_in - q_in * (R1 + 0.5 * R2);
  const double A = std::pow(Kar2_ / R2, 0.25);
  system.C(global_eqn_ids[1]) = -Pavg * A;

  // Eqn (2): WSS - Qin*(R1/Kar1)^0.75 = 0
  const double B = std::pow(R1 / Kar1_, 0.75);
  system.C(global_eqn_ids[2]) = -q_in * B;

  // --------------------
  // dC_dy for rows 0..2
  // --------------------

  // ----- Row 0 -----
  // C0 = -Pd - Qin*Rtot
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -Rtot;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -q_in * dR1_dAs;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[3]) = -q_in * dR2_dAm;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -q_in * dR3_dAmet;

  // ----- Row 1 -----
  // C1 = -Pavg*A,  Pavg = Pin - Qin*(R1 + 0.5*R2),  A = (Kar2/R2)^0.25

  // dC1/dPin = -A
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -A;

  // dC1/dQin = (R1 + 0.5*R2)*A
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1]) =
      (R1 + 0.5 * R2) * A;

  // dC1/dAshear = Q * dR1/dAs * A
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[2]) =
      q_in * dR1_dAs * A;

  // dC1/dAmyo:
  //   dPavg/dAm = -Q*(0.5*dR2/dAm)
  //   dA/dAm = (-1/4)*A*(dR2/dAm)/R2
  //   dC1/dAm = -(dPavg/dAm)*A - Pavg*(dA/dAm)
  //          = (Q/2)*A*dR2/dAm + Pavg*(1/4)*A*(dR2/dAm)/R2
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[3]) =
      (0.5 * q_in * A * dR2_dAm) +
      (Pavg * (0.25 * A * dR2_dAm / R2));

  // No dependence on Ameta in eqn (1) under your P1,P2 definition.

  // ----- Row 2 -----
  // C2 = -Qin*B, B = (R1/Kar1)^0.75
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[1]) = -B;

  // dC2/dAs = -Q * 0.75*B*(dR1/R1)
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[2]) =
      -q_in * (0.75 * B * (dR1_dAs / R1));
}