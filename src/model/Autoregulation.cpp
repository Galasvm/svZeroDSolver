// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Autoregulation.h"

void Autoregulation::setup_dofs(DOFHandler &dofhandler) {
  // 9 equations; 8 internal variables
  // Variable order: [0]=Pin [1]=Qin [2]=Ashear [3]=Amyo [4]=Ameta
  //                 [5]=xshear [6]=xmyo [7]=xmeta [8]=T [9]=WSS
  Block::setup_dofs_(dofhandler, 9,
                     {"Ashear", "Amyo", "Ameta", "xshear", "xmyo", "xmeta",
                      "T", "WSS"});
}

void Autoregulation::update_constant(SparseSystem &system,
                                     std::vector<double> &parameters) {
  const double R        = parameters[global_param_ids[0]];
  const double Qt       = parameters[global_param_ids[1]];
  const double Pt       = parameters[global_param_ids[2]];
  const double Gshear   = parameters[global_param_ids[3]];
  const double TAUshear = parameters[global_param_ids[4]];
  const double Gmyo     = parameters[global_param_ids[5]];
  const double TAUmyo   = parameters[global_param_ids[6]];
  const double Gmeta    = parameters[global_param_ids[7]];
  const double TAUmeta  = parameters[global_param_ids[8]];
  const double lower_frac = parameters[global_param_ids[10]];
  const double upper_frac = parameters[global_param_ids[11]];

  if (!initialized_) {
    const double R1_0 = 0.20 * R;
    const double R2_0 = 0.20 * R;
    const double R3_0 = 0.45 * R;

    R4_   = 0.05 * R;
    R1_0_ = R1_0;

    R1L_ = lower_frac * R1_0;  R1U_ = upper_frac * R1_0;
    R2L_ = lower_frac * R2_0;  R2U_ = upper_frac * R2_0;
    R3L_ = lower_frac * R3_0;  R3U_ = upper_frac * R3_0;

    Kar1_ = std::pow(0.02, 4) * R1_0;  // 0.02 = baseline radius, shear layer
    Kar2_ = std::pow(0.01, 4) * R2_0;  // 0.01 = baseline radius, myo layer
    WSSt_ = Qt / std::pow(0.02, 3);

    const double P1_0 = Pt - Qt * R1_0;
    const double P2_0 = P1_0 - Qt * R2_0;
    Tt_ = 0.5 * (P1_0 + P2_0) * std::pow(Kar2_ / R2_0, 0.25);

    // Sigmoid re-centering: k_ = exp(-C), C = -ln[(1-lower_frac)/(upper_frac-1)]
    k_ = (1.0 - lower_frac) / (upper_frac - 1.0);

    initialized_ = true;
  }

  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) =  1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[8]) =  1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[9]) =  1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[5]) =  Gshear;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[6]) = -Gmyo;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[7]) = -Gmeta;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[5]) =  1.0;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[9]) = -1.0 / WSSt_;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[6]) =  1.0;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[8]) = -1.0 / Tt_;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[7]) =  1.0;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[1]) = -1.0 / Qt;

  system.E.coeffRef(global_eqn_ids[3], global_var_ids[2]) =  1.0;
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[3]) =  1.0;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[4]) =  1.0;
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[5]) =  TAUshear;
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[6]) =  TAUmyo;
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[7]) =  TAUmeta;

  system.C.coeffRef(global_eqn_ids[6]) = 1.0;
  system.C.coeffRef(global_eqn_ids[7]) = 1.0;
  system.C.coeffRef(global_eqn_ids[8]) = 1.0;
}

void Autoregulation::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {

  const double p_in = y[global_var_ids[0]];
  const double q_in = y[global_var_ids[1]];
  const double As   = y[global_var_ids[2]];
  const double Am   = y[global_var_ids[3]];
  const double Amet = y[global_var_ids[4]];

  const double Pd   = parameters[global_param_ids[9]];

  const double eS   = k_ * std::exp(As);
  const double eM   = k_ * std::exp(Am);
  const double eMet = k_ * std::exp(Amet);

  // One-sided shear: tanh gate pins R1 at R1_0 for As >= 0 (dilation only)
  const double R1_sig = (R1L_ + R1U_ * eS) / (1.0 + eS);
  const double th     = std::tanh(50.0 * As);
  const double gate   = 0.5 * (1.0 - th);
  const double dgate  = -25.0 * (1.0 - th * th);

  const double R1   = gate * R1_sig + (1.0 - gate) * R1_0_;
  const double R2   = (R2L_ + R2U_ * eM)   / (1.0 + eM);
  const double R3   = (R3L_ + R3U_ * eMet) / (1.0 + eMet);
  const double Rtot = R1 + R2 + R3 + R4_;

  const double dR1_dAs   = dgate * (R1_sig - R1_0_) +
                           gate * (R1U_ - R1L_) * eS / ((1.0 + eS) * (1.0 + eS));
  const double dR2_dAm   = (R2U_ - R2L_) * eM   / ((1.0 + eM)   * (1.0 + eM));
  const double dR3_dAmet = (R3U_ - R3L_) * eMet / ((1.0 + eMet) * (1.0 + eMet));

  // Eqn (0): Pin - Pd - Qin*Rtot = 0
  system.C(global_eqn_ids[0]) = -Pd - q_in * Rtot;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -Rtot;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -q_in * dR1_dAs;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[3]) = -q_in * dR2_dAm;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -q_in * dR3_dAmet;

  // Eqn (1): T - Pavg*(Kar2/R2)^0.25 = 0,  Pavg = Pin - Qin*(R1 + 0.5*R2)
  const double Pavg = p_in - q_in * (R1 + 0.5 * R2);
  const double A    = std::pow(Kar2_ / R2, 0.25);
  system.C(global_eqn_ids[1]) = -Pavg * A;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -A;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1]) =  (R1 + 0.5 * R2) * A;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[2]) =  q_in * dR1_dAs * A;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[3]) =
      (0.5 * q_in * A * dR2_dAm) + (Pavg * 0.25 * A * dR2_dAm / R2);

  // Eqn (2): WSS - Qin*(R1/Kar1)^0.75 = 0
  const double B = std::pow(R1 / Kar1_, 0.75);
  system.C(global_eqn_ids[2]) = -q_in * B;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[1]) = -B;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[2]) = -q_in * 0.75 * B * dR1_dAs / R1;
}
