// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "AutoregulationRCR.h"

void AutoregulationRCR::setup_dofs(DOFHandler &dofhandler) {
  // 10 equations, 9 internal variables (+ 2 interface = 11 total DOFs)
  // Internal variable order determines global_var_ids[2..10]
  Block::setup_dofs_(dofhandler, 10,
                     {"Ashear", "Amyo", "Ameta", "xshear", "xmyo", "xmeta",
                      "T", "WSS", "Pc"});
}

void AutoregulationRCR::update_constant(SparseSystem &system,
                                        std::vector<double> &parameters) {
  // Parameter indices
  const double R        = parameters[global_param_ids[0]];
  const double Qt       = parameters[global_param_ids[1]];
  const double Pt       = parameters[global_param_ids[2]];
  const double Gshear   = parameters[global_param_ids[3]];
  const double TAUshear = parameters[global_param_ids[4]];
  const double Gmyo     = parameters[global_param_ids[5]];
  const double TAUmyo   = parameters[global_param_ids[6]];
  const double Gmeta    = parameters[global_param_ids[7]];
  const double TAUmeta  = parameters[global_param_ids[8]];
  (void)parameters[global_param_ids[9]];   // Pd  — used in update_solution
  const double Rp       = parameters[global_param_ids[10]];
  (void)parameters[global_param_ids[11]];  // C   — used in update_solution
  const double lower_frac = parameters[global_param_ids[12]];  // default 0.70
  const double upper_frac = parameters[global_param_ids[13]];  // default 1.30

  // --------------------
  // One-time precompute
  // --------------------
  if (!initialized_) {
    const double R1_0 = 0.10 * R;
    const double R2_0 = 0.35 * R;
    const double R3_0 = 0.40 * R;

    R4_  = 0.15 * R;

    R1L_ = lower_frac * R1_0;  R1U_ = upper_frac * R1_0;
    R2L_ = lower_frac * R2_0;  R2U_ = upper_frac * R2_0;
    R3L_ = lower_frac * R3_0;  R3U_ = upper_frac * R3_0;

    Kar1_ = std::pow(0.02, 4) * R1_0;
    Kar2_ = std::pow(0.01, 4) * R2_0;

    WSSt_ = Qt / std::pow(0.02, 3);

    // Baseline pressures for target tension
    // Pc_0 = Pt - Rp*Qt  (pressure after the proximal resistance at baseline)
    const double Pc_0  = Pt - Rp * Qt;
    const double P1_0  = Pc_0 - Qt * R1_0;
    const double P2_0  = P1_0 - Qt * R2_0;
    const double Pavg0 = 0.5 * (P1_0 + P2_0);

    Tt_ = Pavg0 * std::pow(Kar2_ / R2_0, 0.25);

    initialized_ = true;
  }

  // --------------------
  // Variable index mapping (for readability)
  //   global_var_ids: [0]=Pin [1]=Qin [2]=Ashear [3]=Amyo [4]=Ameta
  //                   [5]=xshear [6]=xmyo [7]=xmeta [8]=T [9]=WSS [10]=Pc
  // --------------------

  // Eqn (0): Pin - Pc - Rp*Qin = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0])  =  1.0;   // +Pin
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1])  = -Rp;    // -Rp*Qin
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[10]) = -1.0;   // -Pc

  // Eqn (1): Rtot*(Qin - C*dPc/dt) - Pc + Pd = 0
  // Only the constant linear part goes here; nonlinear terms are in update_solution
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[10]) = -1.0;   // -Pc

  // Eqn (2): T - Pavg*(Kar2/R2)^0.25 = 0
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[8])  =  1.0;   // +T

  // Eqn (3): WSS - Qin*(R1/Kar1)^0.75 = 0
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[9])  =  1.0;   // +WSS

  // Eqn (4): dAshear/dt + Gshear*xshear = 0  (same sign convention as Autoregulation)
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[5])  =  Gshear;

  // Eqn (5): dAmyo/dt - Gmyo*xmyo = 0
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[6])  = -Gmyo;

  // Eqn (6): dAmeta/dt - Gmeta*xmeta = 0
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[7])  = -Gmeta;

  // Eqn (7): TAUshear*dxshear/dt + xshear - WSS/WSSt + 1 = 0
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[5])  =  1.0;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[9])  = -1.0 / WSSt_;

  // Eqn (8): TAUmyo*dxmyo/dt + xmyo - T/Tt + 1 = 0
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[6])  =  1.0;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[8])  = -1.0 / Tt_;

  // Eqn (9): TAUmeta*dxmeta/dt + xmeta - Qin/Qt + 1 = 0
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[7])  =  1.0;
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[1])  = -1.0 / Qt;

  // --------------------
  // E matrix: time-derivative terms
  // --------------------
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[2])  =  1.0;       // dAshear/dt
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[3])  =  1.0;       // dAmyo/dt
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[4])  =  1.0;       // dAmeta/dt
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[5])  =  TAUshear;  // TAUshear*dxshear/dt
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[6])  =  TAUmyo;    // TAUmyo*dxmyo/dt
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[7])  =  TAUmeta;   // TAUmeta*dxmeta/dt

  // --------------------
  // Constant C offsets for eqns (7..9): +1
  // --------------------
  system.C.coeffRef(global_eqn_ids[7]) = 1.0;
  system.C.coeffRef(global_eqn_ids[8]) = 1.0;
  system.C.coeffRef(global_eqn_ids[9]) = 1.0;
}

void AutoregulationRCR::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {

  // States
  const double q_in  = y[global_var_ids[1]];
  const double As    = y[global_var_ids[2]];
  const double Am    = y[global_var_ids[3]];
  const double Amet  = y[global_var_ids[4]];
  const double p_c   = y[global_var_ids[10]];
  const double dp_c  = dy[global_var_ids[10]];  // dPc/dt

  // Parameters
  const double Pd    = parameters[global_param_ids[9]];
  const double C_cap = parameters[global_param_ids[11]];  // capacitance

  // exp(A) terms (sigmoid basis for bounded resistances)
  const double eS   = std::exp(As);
  const double eM   = std::exp(Am);
  const double eMet = std::exp(Amet);

  // Effective regulated resistors
  const double R1 = (R1L_ + R1U_ * eS)   / (1.0 + eS);
  const double R2 = (R2L_ + R2U_ * eM)   / (1.0 + eM);
  const double R3 = (R3L_ + R3U_ * eMet) / (1.0 + eMet);
  const double Rtot = R1 + R2 + R3 + R4_;

  // dRi/dAi
  const double dR1_dAs   = (R1U_ - R1L_) * eS   / ((1.0 + eS)   * (1.0 + eS));
  const double dR2_dAm   = (R2U_ - R2L_) * eM   / ((1.0 + eM)   * (1.0 + eM));
  const double dR3_dAmet = (R3U_ - R3L_) * eMet / ((1.0 + eMet) * (1.0 + eMet));

  // Flow through the regulated segment: Qout = Qin - C*dPc/dt
  const double q_out = q_in - C_cap * dp_c;

  // --------------------
  // Nonlinear C entries
  // --------------------

  // Eqn (1): Rtot*q_out - Pc + Pd = 0
  //   F already contributes -Pc; C gets the remaining nonlinear part
  system.C(global_eqn_ids[1]) = Rtot * q_out + Pd;

  // Eqn (2): T - Pavg*(Kar2/R2)^0.25 = 0
  //   Pavg = Pc - Qin*(R1 + 0.5*R2)
  const double Pavg = p_c - q_in * (R1 + 0.5 * R2);
  const double A    = std::pow(Kar2_ / R2, 0.25);
  system.C(global_eqn_ids[2]) = -Pavg * A;

  // Eqn (3): WSS - Qin*(R1/Kar1)^0.75 = 0
  const double B = std::pow(R1 / Kar1_, 0.75);
  system.C(global_eqn_ids[3]) = -q_in * B;

  // --------------------
  // dC_dy  (Jacobian w.r.t. y)
  // --------------------

  // --- Row 1: C[1] = Rtot*q_out + Pd ---
  // d/dQin: dq_out/dQin = 1
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1]) =  Rtot;
  // d/dAshear: dRtot/dAshear * q_out
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[2]) =  dR1_dAs   * q_out;
  // d/dAmyo:   dRtot/dAmyo   * q_out
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[3]) =  dR2_dAm   * q_out;
  // d/dAmeta:  dRtot/dAmeta  * q_out
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) =  dR3_dAmet * q_out;

  // --- Row 2: C[2] = -Pavg*A,  Pavg = Pc - Qin*(R1 + 0.5*R2),  A = (Kar2/R2)^0.25 ---
  // d/dPc = -A  (dPavg/dPc = 1)
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[10]) = -A;
  // d/dQin = (R1 + 0.5*R2)*A  (dPavg/dQin = -(R1+0.5*R2))
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[1])  = (R1 + 0.5 * R2) * A;
  // d/dAshear: dPavg/dAs = -Qin*dR1/dAs  (A doesn't depend on Ashear)
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[2])  = q_in * dR1_dAs * A;
  // d/dAmyo:   dPavg/dAm = -Qin*(0.5*dR2/dAm),  dA/dAm = (-1/4)*A*(dR2/dAm)/R2
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[3])  =
      (0.5 * q_in * A * dR2_dAm) +
      (Pavg * (0.25 * A * dR2_dAm / R2));

  // --- Row 3: C[3] = -Qin*B,  B = (R1/Kar1)^0.75 ---
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[1])  = -B;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[2])  =
      -q_in * (0.75 * B * (dR1_dAs / R1));

  // --------------------
  // dC_dydot  (Jacobian w.r.t. dy/dt)
  // --------------------

  // Eqn (1): C[1] contains -Rtot*C_cap*dPc/dt  (via q_out = Qin - C_cap*dPc/dt)
  // d C[1] / d(dPc/dt) = -Rtot * C_cap
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[10]) = -Rtot * C_cap;
}
