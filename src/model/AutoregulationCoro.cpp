// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "AutoregulationCoro.h"

#include "Model.h"

// Parameter index aliases for readability
// (matches constructor order, accounting for t+Pim pair occupying one slot)
static constexpr int P_RA1     =  0;
static constexpr int P_RV1     =  1;
static constexpr int P_CA      =  2;
static constexpr int P_CIM     =  3;
static constexpr int P_PIM     =  4;  // current time-interpolated Pim
static constexpr int P_PV      =  5;
static constexpr int P_RA2     =  6;  // baseline regulated resistance
static constexpr int P_QT      =  7;
static constexpr int P_PT      =  8;
static constexpr int P_GSHEAR  =  9;
static constexpr int P_TAUSH   = 10;
static constexpr int P_GMYO    = 11;
static constexpr int P_TAUMYO  = 12;
static constexpr int P_GMETA   = 13;
static constexpr int P_TAUMET  = 14;
static constexpr int P_LOWER   = 15;  // lower_frac, default 0.70
static constexpr int P_UPPER   = 16;  // upper_frac, default 1.30

// Variable index aliases
// y^e = [Pin, Qin, Vim, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS]
static constexpr int V_PIN    =  0;
static constexpr int V_QIN    =  1;
static constexpr int V_VIM    =  2;
static constexpr int V_AS     =  3;
static constexpr int V_AM     =  4;
static constexpr int V_AMET   =  5;
static constexpr int V_XSHEAR =  6;
static constexpr int V_XMYO   =  7;
static constexpr int V_XMET   =  8;
static constexpr int V_T      =  9;
static constexpr int V_WSS    = 10;

void AutoregulationCoro::setup_dofs(DOFHandler &dofhandler) {
  // 10 equations; 9 internal variables (Vim + 8 autoregulation vars)
  Block::setup_dofs_(dofhandler, 10,
                     {"volume_im", "Ashear", "Amyo", "Ameta",
                      "xshear", "xmyo", "xmeta", "T", "WSS"});
}

void AutoregulationCoro::setup_initial_state_dependent_params(
    State initial_state, std::vector<double> &parameters) {
  const double P_in     = initial_state.y   [global_var_ids[V_PIN]];
  const double Q_in     = initial_state.y   [global_var_ids[V_QIN]];
  const double P_in_dot = initial_state.ydot[global_var_ids[V_PIN]];
  const double Q_in_dot = initial_state.ydot[global_var_ids[V_QIN]];
  const double Ra  = parameters[global_param_ids[P_RA1]];
  const double Ra2 = parameters[global_param_ids[P_RA2]];  // baseline
  const double Ca  = parameters[global_param_ids[P_CA]];

  // Pressure proximal to Ca, distal to Ra1
  const double P_Ca     = P_in - Ra * Q_in;
  const double P_Ca_dot = P_in_dot - Ra * Q_in_dot;
  // Flow into Ra2 (inflow minus flow into Ca)
  const double Q_am     = Q_in - Ca * P_Ca_dot;
  // Pressure proximal to Cim and distal to Ra2 (at baseline, Ram = Ra2)
  P_Cim_0 = P_Ca - Ra2 * Q_am;
  // Initial intramyocardial pressure (time-interpolated value at t=0)
  Pim_0 = parameters[global_param_ids[P_PIM]];
}

void AutoregulationCoro::update_constant(
    SparseSystem &system, std::vector<double> &parameters) {

  const double Ra       = parameters[global_param_ids[P_RA1]];
  const double Rv       = parameters[global_param_ids[P_RV1]];
  const double Ca       = parameters[global_param_ids[P_CA]];
  const double Cim      = parameters[global_param_ids[P_CIM]];
  const double Ra2      = parameters[global_param_ids[P_RA2]];
  const double Qt       = parameters[global_param_ids[P_QT]];
  const double Pt       = parameters[global_param_ids[P_PT]];
  const double Gshear   = parameters[global_param_ids[P_GSHEAR]];
  const double TAUshear = parameters[global_param_ids[P_TAUSH]];
  const double Gmyo     = parameters[global_param_ids[P_GMYO]];
  const double TAUmyo   = parameters[global_param_ids[P_TAUMYO]];
  const double Gmeta    = parameters[global_param_ids[P_GMETA]];
  const double TAUmeta  = parameters[global_param_ids[P_TAUMET]];
  const double lower_frac = parameters[global_param_ids[P_LOWER]];  // default 0.70
  const double upper_frac = parameters[global_param_ids[P_UPPER]];  // default 1.30

  // ------------------------------------------------------------------
  // One-time precompute of autoregulation constants
  // ------------------------------------------------------------------
  if (!initialized_) {
    const double R1_0 = 0.10 * Ra2;
    const double R2_0 = 0.35 * Ra2;
    const double R3_0 = 0.40 * Ra2;

    R4_  = 0.15 * Ra2;

    R1L_ = lower_frac * R1_0;  R1U_ = upper_frac * R1_0;
    R2L_ = lower_frac * R2_0;  R2U_ = upper_frac * R2_0;
    R3L_ = lower_frac * R3_0;  R3U_ = upper_frac * R3_0;

    Kar1_ = std::pow(0.02, 4) * R1_0;
    Kar2_ = std::pow(0.01, 4) * R2_0;

    WSSt_ = Qt / std::pow(0.02, 3);

    // Baseline Pavg for target myogenic tension.
    // Pt is Pin from the baseline simulation (pressure entering the block).
    // Subtract Ra1*Qt to get the pressure at the node between Ra1 and the
    // regulated segment, then step through R1 and R2 to find Pavg.
    const double Pa_0  = Pt - Ra * Qt;             // pressure after Ra1
    const double P1_0  = Pa_0 - Qt * R1_0;
    const double P2_0  = P1_0 - Qt * R2_0;
    const double Pavg0 = 0.5 * (P1_0 + P2_0);
    Tt_ = Pavg0 * std::pow(Kar2_ / R2_0, 0.25);

    initialized_ = true;
  }

  // ------------------------------------------------------------------
  // Coronary equations 0 and 1
  // Different formulation for steady (initial-condition) mode
  // ------------------------------------------------------------------
  if (steady) {
    // Steady eqn 0: Vim = 0
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[V_VIM]) =  1.0;
    // Steady eqn 1: -Pin + (Ra+Rv)*Qin + C[1] = 0
    //   with C[1] = Pv + Rtot*Qin handled in update_solution
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_PIN]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN]) = Ra + Rv;
  } else {
    // Non-steady eqn 0 (Cim*Rv*Qin - Vim + ... = 0)
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[V_QIN]) =  Cim * Rv;
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[V_VIM]) = -1.0;

    system.E.coeffRef(global_eqn_ids[0], global_var_ids[V_PIN]) = -Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[V_QIN]) =  Ra * Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[V_VIM]) = -Cim * Rv;

    // Non-steady eqn 1 (constant linear parts; Rtot contributions in update_solution)
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_PIN]) =  Cim * Rv;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN]) = -Cim * Rv * Ra;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM]) = -Rv;
    // Note: no E[1, Vim] here — the Rtot*dVim/dt contribution is in dC_dydot
  }

  // ------------------------------------------------------------------
  // Autoregulation equations 2–9 (identical to Autoregulation block,
  // just shifted by 2 in equation index)
  // ------------------------------------------------------------------

  // Eqn 2: T - Pavg*(Kar2/R2)^0.25 = 0
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[V_T])     =  1.0;

  // Eqn 3: WSS - Qin*(R1/Kar1)^0.75 = 0
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[V_WSS])   =  1.0;

  // Eqn 4: dAshear/dt + Gshear*xshear = 0
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[V_XSHEAR]) =  Gshear;

  // Eqn 5: dAmyo/dt - Gmyo*xmyo = 0
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[V_XMYO])   = -Gmyo;

  // Eqn 6: dAmeta/dt - Gmeta*xmeta = 0
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[V_XMET])   = -Gmeta;

  // Eqn 7: TAUshear*dxshear/dt + xshear - WSS/WSSt + 1 = 0
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[V_XSHEAR]) =  1.0;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[V_WSS])    = -1.0 / WSSt_;

  // Eqn 8: TAUmyo*dxmyo/dt + xmyo - T/Tt + 1 = 0
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[V_XMYO])   =  1.0;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[V_T])      = -1.0 / Tt_;

  // Eqn 9: TAUmeta*dxmeta/dt + xmeta - Qin/Qt + 1 = 0
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[V_XMET])   =  1.0;
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[V_QIN])    = -1.0 / Qt;

  // E matrix: autoregulation time-derivative terms
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[V_AS])     =  1.0;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[V_AM])     =  1.0;
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[V_AMET])   =  1.0;
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[V_XSHEAR]) =  TAUshear;
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[V_XMYO])   =  TAUmyo;
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[V_XMET])   =  TAUmeta;

  // Constant +1 offsets for eqns 7–9
  system.C.coeffRef(global_eqn_ids[7]) = 1.0;
  system.C.coeffRef(global_eqn_ids[8]) = 1.0;
  system.C.coeffRef(global_eqn_ids[9]) = 1.0;
}

void AutoregulationCoro::update_time(SparseSystem &system,
                                           std::vector<double> &parameters) {
  if (steady) {
    // Nothing here: C[1] is handled entirely in update_solution
    return;
  }

  const double Pim = parameters[global_param_ids[P_PIM]];
  const double Pv  = parameters[global_param_ids[P_PV]];
  const double Cim = parameters[global_param_ids[P_CIM]];

  // C[0]: time-varying part of coronary eqn 0 (no Rtot dependence)
  system.C(global_eqn_ids[0]) = Cim * (-Pim + Pv + Pim_0 - P_Cim_0);
  // C[1]: set entirely in update_solution (needs Rtot and Vim)
}

void AutoregulationCoro::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {

  // Interface and internal states
  const double p_in = y[global_var_ids[V_PIN]];
  const double q_in = y[global_var_ids[V_QIN]];
  const double vim  = y[global_var_ids[V_VIM]];
  const double As   = y[global_var_ids[V_AS]];
  const double Am   = y[global_var_ids[V_AM]];
  const double Amet = y[global_var_ids[V_AMET]];
  const double dvim = dy[global_var_ids[V_VIM]];  // dVim/dt

  // Parameters needed for eqn 1
  const double Pim = parameters[global_param_ids[P_PIM]];
  const double Pv  = parameters[global_param_ids[P_PV]];
  const double Cim = parameters[global_param_ids[P_CIM]];
  const double Rv  = parameters[global_param_ids[P_RV1]];

  // Sigmoid activations
  const double eS   = std::exp(As);
  const double eM   = std::exp(Am);
  const double eMet = std::exp(Amet);

  // Regulated resistors
  const double R1   = (R1L_ + R1U_ * eS)   / (1.0 + eS);
  const double R2   = (R2L_ + R2U_ * eM)   / (1.0 + eM);
  const double R3   = (R3L_ + R3U_ * eMet) / (1.0 + eMet);
  const double Rtot = R1 + R2 + R3 + R4_;

  // dRi/dAi (sigmoid derivative)
  const double dR1_dAs   = (R1U_ - R1L_) * eS   / ((1.0 + eS)   * (1.0 + eS));
  const double dR2_dAm   = (R2U_ - R2L_) * eM   / ((1.0 + eM)   * (1.0 + eM));
  const double dR3_dAmet = (R3U_ - R3L_) * eMet / ((1.0 + eMet) * (1.0 + eMet));

  // ------------------------------------------------------------------
  // Eqn 1: nonlinear C and Jacobian
  //   Steady:     C[1] = Pv + Rtot*Qin
  //   Non-steady: C[1] = -Cim*Rv*Pim_offset
  //                      + Rtot*(-Vim + Cim*(Pv - Pim_offset) - Cim*Rv*dVim/dt)
  // All five dC_dy entries (Qin, Vim, Ashear, Amyo, Ameta) are always written
  // so the sparsity pattern is consistent across steady/non-steady modes.
  // ------------------------------------------------------------------
  if (steady) {
    system.C(global_eqn_ids[1]) = Pv + Rtot * q_in;

    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN])  = Rtot;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM])  = 0.0;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AS])   = dR1_dAs   * q_in;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AM])   = dR2_dAm   * q_in;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AMET]) = dR3_dAmet * q_in;

    system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM]) = 0.0;
  } else {
    const double pim_offset = Pim + P_Cim_0 - Pim_0;
    // Bracket that Rtot multiplies: captures Vim, Pv-Pim, and dVim/dt terms
    const double ram_factor = -vim + Cim * (Pv - pim_offset) - Cim * Rv * dvim;

    system.C(global_eqn_ids[1]) = -Cim * Rv * pim_offset + Rtot * ram_factor;

    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN])  = 0.0;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM])  = -Rtot;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AS])   = dR1_dAs   * ram_factor;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AM])   = dR2_dAm   * ram_factor;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AMET]) = dR3_dAmet * ram_factor;

    system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM]) = -Cim * Rv * Rtot;
  }

  // ------------------------------------------------------------------
  // Eqn 2: T - Pavg*(Kar2/R2)^0.25 = 0
  //   Pa   = Pin - Qin*Ra1          (pressure after Ra1, entering regulated seg)
  //   Pavg = Pa  - Qin*(R1 + 0.5*R2)
  //        = Pin - Qin*(Ra1 + R1 + 0.5*R2)
  // Using Pin from the baseline simulation as Pt keeps this consistent with
  // the RCR block and the user's workflow.
  // ------------------------------------------------------------------
  const double Ra   = parameters[global_param_ids[P_RA1]];
  const double Pavg = p_in - q_in * (Ra + R1 + 0.5 * R2);
  const double A    = std::pow(Kar2_ / R2, 0.25);
  system.C(global_eqn_ids[2]) = -Pavg * A;

  // dC[2]/dPin = -A  (Pa = Pin - Ra*Qin; dPa/dPin = 1)
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_PIN]) = -A;
  // dC[2]/dQin = (Ra + R1 + 0.5*R2)*A
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_QIN]) = (Ra + R1 + 0.5 * R2) * A;
  // dC[2]/dAshear: dPavg/dAs = -Qin*dR1/dAs, A doesn't depend on Ashear
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_AS])  = q_in * dR1_dAs * A;
  // dC[2]/dAmyo: dPavg/dAm = -Qin*0.5*dR2/dAm, dA/dAm = (-1/4)*A*(dR2/R2)/dAm
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_AM])  =
      (0.5 * q_in * A * dR2_dAm) + (Pavg * 0.25 * A * dR2_dAm / R2);

  // ------------------------------------------------------------------
  // Eqn 3: WSS - Qin*(R1/Kar1)^0.75 = 0
  // ------------------------------------------------------------------
  const double B = std::pow(R1 / Kar1_, 0.75);
  system.C(global_eqn_ids[3]) = -q_in * B;

  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[V_QIN]) = -B;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[V_AS])  =
      -q_in * (0.75 * B * dR1_dAs / R1);
}
