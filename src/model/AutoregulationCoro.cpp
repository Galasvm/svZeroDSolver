// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "AutoregulationCoro.h"

#include "Model.h"

void AutoregulationCoro::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 12,
                     {"volume_im", "Ashear", "Amyo", "Ameta",
                      "xshear", "xmyo", "xmeta", "T", "WSS", "Pa", "q_micro"});
}

void AutoregulationCoro::setup_initial_state_dependent_params(
    State initial_state, std::vector<double> &parameters) {
  auto P_in     = initial_state.y   [global_var_ids[0]];   // Pin
  auto Q_in     = initial_state.y   [global_var_ids[1]];   // Qin
  auto P_in_dot = initial_state.ydot[global_var_ids[0]];   // Pin
  auto Q_in_dot = initial_state.ydot[global_var_ids[1]];   // Qin
  auto Ra  = parameters[global_param_ids[0]];   // Ra1
  auto Ra2 = parameters[global_param_ids[6]];   // Ra2
  auto Ca  = parameters[global_param_ids[2]];   // Ca

  // Pa = Pin - Ra*Qin  (full Ra drop before Ca)
  auto P_Ca     = P_in - Ra * Q_in;
  auto P_Ca_dot = P_in_dot - Ra * Q_in_dot;
  auto Q_am     = Q_in - Ca * P_Ca_dot;
  P_Cim_0 = P_Ca - Ra2 * Q_am;
  Pim_0   = parameters[global_param_ids[4]];   // Pim
}

void AutoregulationCoro::update_constant(
    SparseSystem &system, std::vector<double> &parameters) {

  auto Ra       = parameters[global_param_ids[0]];    // Ra1
  auto Rv       = parameters[global_param_ids[1]];    // Rv1
  auto Ca       = parameters[global_param_ids[2]];    // Ca
  auto Cim      = parameters[global_param_ids[3]];    // Cim
  auto Ra2      = parameters[global_param_ids[6]];    // Ra2
  auto Qt       = parameters[global_param_ids[7]];    // Qt
  auto Pt       = parameters[global_param_ids[8]];    // Pt
  auto Gshear   = parameters[global_param_ids[9]];    // Gshear
  auto TAUshear = parameters[global_param_ids[10]];   // taushear
  auto Gmyo     = parameters[global_param_ids[11]];   // Gmyo
  auto TAUmyo   = parameters[global_param_ids[12]];   // taumyo
  auto Gmeta    = parameters[global_param_ids[13]];   // Gmeta
  auto TAUmeta  = parameters[global_param_ids[14]];   // taumeta
  auto lower_frac = parameters[global_param_ids[15]]; // lower_frac
  auto upper_frac = parameters[global_param_ids[16]]; // upper_frac

  if (!initialized_) {
    // Ra split: 30% static, 16% shear-regulated, 54% myo-regulated (all before Ca)
    Ra1_static_ = 0.30 * Ra;
    auto Ra1_shear_0 = 0.16 * Ra;
    auto Ra1_myo_0   = 0.54 * Ra;
    Ra1_shear_0_ = Ra1_shear_0;
    Ra1SL_ = lower_frac * Ra1_shear_0;  Ra1SU_ = upper_frac * Ra1_shear_0;
    Ra1ML_ = lower_frac * Ra1_myo_0;    Ra1MU_ = upper_frac * Ra1_myo_0;

    // Ra2: 100% metabolic (after Ca)
    Ra2L_ = lower_frac * Ra2;  Ra2U_ = upper_frac * Ra2;

    // Geometric constants
    Kar1_shear_ = std::pow(0.02, 4) * Ra1_shear_0;  // for WSS of Ra1_shear
    Kar_myo_    = std::pow(0.01, 4) * Ra1_myo_0;    // for tension of Ra1_myo
    WSSt_ = Qt / std::pow(0.02, 3);

    // Tt: target tension = Pavg0 * (Kar_myo/Ra1_myo_0)^0.25
    // Pavg0 = average pressure across Ra1_myo at baseline
    //       = Pa_0 + 0.5*Ra1_myo_0*Qt  (Pa_0 = pressure at Ca node = Pt - Ra*Qt)
    auto Pa_0   = Pt - Ra * Qt;
    auto Pavg0  = Pa_0 + 0.5 * Ra1_myo_0 * Qt;
    Tt_ = Pavg0 * std::pow(Kar_myo_ / Ra1_myo_0, 0.25);

    // Sigmoid re-centering: k_ = exp(-C), C = -ln[(1-lower_frac)/(upper_frac-1)]
    k_ = (1.0 - lower_frac) / (upper_frac - 1.0);

    initialized_ = true;
  }

  // Coronary eqns 0 & 1: Ra1_static is the constant part of Ra_eff before Ca
  if (steady) {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) =  1.0;   // Vim
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -1.0;   // Pin
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) =  Ra1_static_ + Rv;  // Qin
  } else {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) =  Cim * Rv;   // Qin
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;        // Vim
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) =  Cim * Rv;   // Pin
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = -Cim * Rv * Ra1_static_;  // Qin
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -Rv;         // Vim

    system.E.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -Ca * Cim * Rv;              // Pin
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) =  Ra1_static_ * Ca * Cim * Rv; // Qin
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -Cim * Rv;                    // Vim
  }

  system.F.coeffRef(global_eqn_ids[2], global_var_ids[9])  =  1.0;         // T
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[10]) =  1.0;         // WSS
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[6])  =  Gshear;      // xshear
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[7])  = -Gmyo;        // xmyo
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[8])  = -Gmeta;       // xmeta
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[6])  =  1.0;         // xshear
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[10]) = -1.0 / WSSt_; // WSS
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[7])  =  1.0;         // xmyo
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[9])  = -1.0 / Tt_;   // T
  // Eqn 9: metabolic error uses q_micro (flow after Ca)
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[8])  =  1.0;         // xmeta
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[12]) = -1.0 / Qt;    // q_micro

  // Eqn 10: Pa - Pin + Ra1_static*Qin + (Ra1_shear+Ra1_myo)*Qin = 0
  // Constant F part; nonlinear (Ra1_shear+Ra1_myo)*Qin handled in update_solution
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[11]) =  1.0;         // Pa
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[0])  = -1.0;         // Pin
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[1])  =  Ra1_static_; // Qin

  // Eqn 11: q_micro - Qin + Ca*dPa/dt = 0
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[12]) =  1.0;   // q_micro
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[1])  = -1.0;   // Qin
  system.E.coeffRef(global_eqn_ids[11], global_var_ids[11]) =  Ca;    // Pa

  system.E.coeffRef(global_eqn_ids[4], global_var_ids[3])  =  1.0;        // Ashear
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[4])  =  1.0;        // Amyo
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[5])  =  1.0;        // Ameta
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[6])  =  TAUshear;   // xshear
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[7])  =  TAUmyo;     // xmyo
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[8])  =  TAUmeta;    // xmeta

  system.C.coeffRef(global_eqn_ids[7]) = 1.0;
  system.C.coeffRef(global_eqn_ids[8]) = 1.0;
  system.C.coeffRef(global_eqn_ids[9]) = 1.0;
}

void AutoregulationCoro::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {

  auto q_in    = y[global_var_ids[1]];    // Qin
  auto vim     = y[global_var_ids[2]];    // Vim
  auto As      = y[global_var_ids[3]];    // Ashear
  auto Am      = y[global_var_ids[4]];    // Amyo
  auto Amet    = y[global_var_ids[5]];    // Ameta
  auto dvim    = dy[global_var_ids[2]];   // Vim
  auto dqin    = dy[global_var_ids[1]];   // Qin
  auto Pa      = y[global_var_ids[11]];   // Pa
  auto q_micro = y[global_var_ids[12]];   // q_micro

  auto Pim = parameters[global_param_ids[4]];   // Pim
  auto Pv  = parameters[global_param_ids[5]];   // Pv
  auto Ca  = parameters[global_param_ids[2]];   // Ca
  auto Cim = parameters[global_param_ids[3]];   // Cim
  auto Rv  = parameters[global_param_ids[1]];   // Rv1

  auto eS   = k_ * std::exp(As);
  auto eM   = k_ * std::exp(Am);
  auto eMet = k_ * std::exp(Amet);

  // One-sided shear: tanh gate pins Ra1_shear at Ra1_shear_0 for As >= 0
  // (dilation only)
  auto Ra1s_sig = (Ra1SL_ + Ra1SU_ * eS) / (1.0 + eS);
  auto th       = std::tanh(50.0 * As);
  auto gate     = 0.5 * (1.0 - th);
  auto dgate    = -25.0 * (1.0 - th * th);

  // Sigmoid activations
  auto Ra1_shear = gate * Ra1s_sig + (1.0 - gate) * Ra1_shear_0_;  // shear (16% baseline)
  auto Ra1_myo   = (Ra1ML_ + Ra1MU_ * eM)   / (1.0 + eM);  // myo   (54% baseline)
  auto Ra2       = (Ra2L_  + Ra2U_  * eMet) / (1.0 + eMet); // meta  (100% of Ra2)
  auto Rtot      = Ra2;  // only metabolic resistance after Ca

  // Sigmoid derivatives
  auto dRa1s_dAs  = dgate * (Ra1s_sig - Ra1_shear_0_) +
                    gate * (Ra1SU_ - Ra1SL_) * eS / ((1.0 + eS) * (1.0 + eS));
  auto dRa1m_dAm  = (Ra1MU_ - Ra1ML_) * eM   / ((1.0 + eM)   * (1.0 + eM));
  auto dRa2_dAmet = (Ra2U_  - Ra2L_)  * eMet / ((1.0 + eMet) * (1.0 + eMet));

  auto pim_offset = Pim + P_Cim_0 - Pim_0;

  // ------------------------------------------------------------------
  // Coronary eqns 0 & 1: Ra1_shear + Ra1_myo are both nonlinear parts of Ra_eff
  // ------------------------------------------------------------------
  if (steady) {
    system.C(global_eqn_ids[1]) = Pv + (Ra1_shear + Ra1_myo + Rtot) * q_in;

    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1])  =  Ra1_shear + Ra1_myo + Rtot; // Qin
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[2])  =  0.0;                        // Vim
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[3])  =  dRa1s_dAs * q_in;            // Ashear
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4])  =  dRa1m_dAm * q_in;            // Amyo
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[5])  =  dRa2_dAmet * q_in;           // Ameta

    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[3])     = 0.0;   // Ashear
    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4])     = 0.0;   // Amyo
    system.dC_dydot.coeffRef(global_eqn_ids[0], global_var_ids[1])  = 0.0;   // Qin
    system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[2])  = 0.0;   // Vim
  } else {
    auto ram_factor = -vim + Cim * (Pv - pim_offset) - Cim * Rv * dvim;

    system.C(global_eqn_ids[0]) = Cim * (-Pim + Pv + Pim_0 - P_Cim_0)
                                   + (Ra1_shear + Ra1_myo) * Ca * Cim * Rv * dqin;
    system.C(global_eqn_ids[1]) = -Cim * Rv * (Ra1_shear + Ra1_myo) * q_in
                                   - Cim * Rv * pim_offset
                                   + Rtot * ram_factor;

    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[3])     = dRa1s_dAs * Ca * Cim * Rv * dqin; // Ashear
    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4])     = dRa1m_dAm * Ca * Cim * Rv * dqin; // Amyo
    system.dC_dydot.coeffRef(global_eqn_ids[0], global_var_ids[1])  = (Ra1_shear + Ra1_myo) * Ca * Cim * Rv; // Qin

    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1])  = -Cim * Rv * (Ra1_shear + Ra1_myo); // Qin
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[2])  = -Rtot;                             // Vim
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[3])  = -Cim * Rv * dRa1s_dAs * q_in;      // Ashear
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4])  = -Cim * Rv * dRa1m_dAm * q_in;      // Amyo
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[5])  =  dRa2_dAmet * ram_factor;          // Ameta
    system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -Cim * Rv * Rtot;                // Vim
  }

  // Eqn 2 (T): Pavg = average pressure across Ra1_myo = Pa + 0.5*Ra1_myo*Qin
  // (Ra1_myo is the last resistor before Ca; Pa is the pressure at Ca)
  auto Pavg = Pa + 0.5 * Ra1_myo * q_in;
  auto A    = std::pow(Kar_myo_ / Ra1_myo, 0.25);
  system.C(global_eqn_ids[2]) = -Pavg * A;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[11]) = -A;                    // Pa
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[1])  = -0.5 * Ra1_myo * A;    // Qin
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4])  =                        // Amyo
      -0.5 * q_in * A * dRa1m_dAm + Pavg * 0.25 * A * dRa1m_dAm / Ra1_myo;

  // Eqn 3 (WSS): flow through Ra1_shear is Qin (before Ca)
  auto B = std::pow(Ra1_shear / Kar1_shear_, 0.75);
  system.C(global_eqn_ids[3]) = -q_in * B;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[1]) = -B;                                    // Qin
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[3]) = -q_in * 0.75 * B * dRa1s_dAs / Ra1_shear; // Ashear

  // Eqn 10 (Pa): nonlinear part (Ra1_shear + Ra1_myo)*Qin
  system.C(global_eqn_ids[10]) = (Ra1_shear + Ra1_myo) * q_in;
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[1]) =  Ra1_shear + Ra1_myo;  // Qin
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[3]) =  dRa1s_dAs * q_in;      // Ashear
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[4]) =  dRa1m_dAm * q_in;      // Amyo
}
