// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "AutoregulationCoro.h"

#include "Model.h"

static constexpr int P_RA1     =  0;
static constexpr int P_RV1     =  1;
static constexpr int P_CA      =  2;
static constexpr int P_CIM     =  3;
static constexpr int P_PIM     =  4;
static constexpr int P_PV      =  5;
static constexpr int P_RA2     =  6;
static constexpr int P_QT      =  7;
static constexpr int P_PT      =  8;
static constexpr int P_GSHEAR  =  9;
static constexpr int P_TAUSH   = 10;
static constexpr int P_GMYO    = 11;
static constexpr int P_TAUMYO  = 12;
static constexpr int P_GMETA   = 13;
static constexpr int P_TAUMET  = 14;
static constexpr int P_LOWER   = 15;
static constexpr int P_UPPER   = 16;

// y^e = [Pin, Qin, Vim, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS, Pa, q_micro]
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
static constexpr int V_PA     = 11;
static constexpr int V_QMICRO = 12;

void AutoregulationCoro::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 12,
                     {"volume_im", "Ashear", "Amyo", "Ameta",
                      "xshear", "xmyo", "xmeta", "T", "WSS", "Pa", "q_micro"});
}

void AutoregulationCoro::setup_initial_state_dependent_params(
    State initial_state, std::vector<double> &parameters) {
  auto P_in     = initial_state.y   [global_var_ids[V_PIN]];
  auto Q_in     = initial_state.y   [global_var_ids[V_QIN]];
  auto P_in_dot = initial_state.ydot[global_var_ids[V_PIN]];
  auto Q_in_dot = initial_state.ydot[global_var_ids[V_QIN]];
  auto Ra  = parameters[global_param_ids[P_RA1]];
  auto Ra2 = parameters[global_param_ids[P_RA2]];
  auto Ca  = parameters[global_param_ids[P_CA]];

  // Pa = Pin - Ra*Qin  (full Ra drop before Ca)
  auto P_Ca     = P_in - Ra * Q_in;
  auto P_Ca_dot = P_in_dot - Ra * Q_in_dot;
  auto Q_am     = Q_in - Ca * P_Ca_dot;
  P_Cim_0 = P_Ca - Ra2 * Q_am;
  Pim_0   = parameters[global_param_ids[P_PIM]];
}

void AutoregulationCoro::update_constant(
    SparseSystem &system, std::vector<double> &parameters) {

  auto Ra       = parameters[global_param_ids[P_RA1]];
  auto Rv       = parameters[global_param_ids[P_RV1]];
  auto Ca       = parameters[global_param_ids[P_CA]];
  auto Cim      = parameters[global_param_ids[P_CIM]];
  auto Ra2      = parameters[global_param_ids[P_RA2]];
  auto Qt       = parameters[global_param_ids[P_QT]];
  auto Pt       = parameters[global_param_ids[P_PT]];
  auto Gshear   = parameters[global_param_ids[P_GSHEAR]];
  auto TAUshear = parameters[global_param_ids[P_TAUSH]];
  auto Gmyo     = parameters[global_param_ids[P_GMYO]];
  auto TAUmyo   = parameters[global_param_ids[P_TAUMYO]];
  auto Gmeta    = parameters[global_param_ids[P_GMETA]];
  auto TAUmeta  = parameters[global_param_ids[P_TAUMET]];
  auto lower_frac = parameters[global_param_ids[P_LOWER]];
  auto upper_frac = parameters[global_param_ids[P_UPPER]];

  if (!initialized_) {
    // Ra split: 30% static, 16% shear-regulated, 54% myo-regulated (all before Ca)
    Ra1_static_ = 0.30 * Ra;
    auto Ra1_shear_0 = 0.16 * Ra;
    auto Ra1_myo_0   = 0.54 * Ra;
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

    initialized_ = true;
  }

  // Coronary eqns 0 & 1: Ra1_static is the constant part of Ra_eff before Ca
  if (steady) {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[V_VIM]) =  1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_PIN]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN]) =  Ra1_static_ + Rv;
  } else {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[V_QIN]) =  Cim * Rv;
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[V_VIM]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_PIN]) =  Cim * Rv;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN]) = -Cim * Rv * Ra1_static_;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM]) = -Rv;

    system.E.coeffRef(global_eqn_ids[0], global_var_ids[V_PIN]) = -Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[V_QIN]) =  Ra1_static_ * Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[V_VIM]) = -Cim * Rv;
  }

  system.F.coeffRef(global_eqn_ids[2], global_var_ids[V_T])      =  1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[V_WSS])    =  1.0;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[V_XSHEAR]) =  Gshear;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[V_XMYO])   = -Gmyo;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[V_XMET])   = -Gmeta;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[V_XSHEAR]) =  1.0;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[V_WSS])    = -1.0 / WSSt_;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[V_XMYO])   =  1.0;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[V_T])      = -1.0 / Tt_;
  // Eqn 9: metabolic error uses q_micro (flow after Ca)
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[V_XMET])   =  1.0;
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[V_QMICRO]) = -1.0 / Qt;

  // Eqn 10: Pa - Pin + Ra1_static*Qin + (Ra1_shear+Ra1_myo)*Qin = 0
  // Constant F part; nonlinear (Ra1_shear+Ra1_myo)*Qin handled in update_solution
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[V_PA])  =  1.0;
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[V_PIN]) = -1.0;
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[V_QIN]) =  Ra1_static_;

  // Eqn 11: q_micro - Qin + Ca*dPa/dt = 0
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[V_QMICRO]) =  1.0;
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[V_QIN])    = -1.0;
  system.E.coeffRef(global_eqn_ids[11], global_var_ids[V_PA])     =  Ca;

  system.E.coeffRef(global_eqn_ids[4], global_var_ids[V_AS])     =  1.0;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[V_AM])     =  1.0;
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[V_AMET])   =  1.0;
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[V_XSHEAR]) =  TAUshear;
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[V_XMYO])   =  TAUmyo;
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[V_XMET])   =  TAUmeta;

  system.C.coeffRef(global_eqn_ids[7]) = 1.0;
  system.C.coeffRef(global_eqn_ids[8]) = 1.0;
  system.C.coeffRef(global_eqn_ids[9]) = 1.0;
}

void AutoregulationCoro::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {

  auto q_in    = y[global_var_ids[V_QIN]];
  auto vim     = y[global_var_ids[V_VIM]];
  auto As      = y[global_var_ids[V_AS]];
  auto Am      = y[global_var_ids[V_AM]];
  auto Amet    = y[global_var_ids[V_AMET]];
  auto dvim    = dy[global_var_ids[V_VIM]];
  auto dqin    = dy[global_var_ids[V_QIN]];
  auto Pa      = y[global_var_ids[V_PA]];
  auto q_micro = y[global_var_ids[V_QMICRO]];

  auto Pim = parameters[global_param_ids[P_PIM]];
  auto Pv  = parameters[global_param_ids[P_PV]];
  auto Ca  = parameters[global_param_ids[P_CA]];
  auto Cim = parameters[global_param_ids[P_CIM]];
  auto Rv  = parameters[global_param_ids[P_RV1]];

  auto eS   = std::exp(As);
  auto eM   = std::exp(Am);
  auto eMet = std::exp(Amet);

  // Sigmoid activations
  auto Ra1_shear = (Ra1SL_ + Ra1SU_ * eS)   / (1.0 + eS);  // shear (16% baseline)
  auto Ra1_myo   = (Ra1ML_ + Ra1MU_ * eM)   / (1.0 + eM);  // myo   (54% baseline)
  auto Ra2       = (Ra2L_  + Ra2U_  * eMet) / (1.0 + eMet); // meta  (100% of Ra2)
  auto Rtot      = Ra2;  // only metabolic resistance after Ca

  // Sigmoid derivatives
  auto dRa1s_dAs  = (Ra1SU_ - Ra1SL_) * eS   / ((1.0 + eS)   * (1.0 + eS));
  auto dRa1m_dAm  = (Ra1MU_ - Ra1ML_) * eM   / ((1.0 + eM)   * (1.0 + eM));
  auto dRa2_dAmet = (Ra2U_  - Ra2L_)  * eMet / ((1.0 + eMet) * (1.0 + eMet));

  auto pim_offset = Pim + P_Cim_0 - Pim_0;

  // ------------------------------------------------------------------
  // Coronary eqns 0 & 1: Ra1_shear + Ra1_myo are both nonlinear parts of Ra_eff
  // ------------------------------------------------------------------
  if (steady) {
    system.C(global_eqn_ids[1]) = Pv + (Ra1_shear + Ra1_myo + Rtot) * q_in;

    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN])  =  Ra1_shear + Ra1_myo + Rtot;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM])  =  0.0;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AS])   =  dRa1s_dAs * q_in;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AM])   =  dRa1m_dAm * q_in;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AMET]) =  dRa2_dAmet * q_in;

    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[V_AS])     = 0.0;
    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[V_AM])     = 0.0;
    system.dC_dydot.coeffRef(global_eqn_ids[0], global_var_ids[V_QIN]) = 0.0;
    system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM]) = 0.0;
  } else {
    auto ram_factor = -vim + Cim * (Pv - pim_offset) - Cim * Rv * dvim;

    system.C(global_eqn_ids[0]) = Cim * (-Pim + Pv + Pim_0 - P_Cim_0)
                                   + (Ra1_shear + Ra1_myo) * Ca * Cim * Rv * dqin;
    system.C(global_eqn_ids[1]) = -Cim * Rv * (Ra1_shear + Ra1_myo) * q_in
                                   - Cim * Rv * pim_offset
                                   + Rtot * ram_factor;

    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[V_AS])     = dRa1s_dAs * Ca * Cim * Rv * dqin;
    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[V_AM])     = dRa1m_dAm * Ca * Cim * Rv * dqin;
    system.dC_dydot.coeffRef(global_eqn_ids[0], global_var_ids[V_QIN]) = (Ra1_shear + Ra1_myo) * Ca * Cim * Rv;

    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_QIN])  = -Cim * Rv * (Ra1_shear + Ra1_myo);
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM])  = -Rtot;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AS])   = -Cim * Rv * dRa1s_dAs * q_in;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AM])   = -Cim * Rv * dRa1m_dAm * q_in;
    system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[V_AMET]) =  dRa2_dAmet * ram_factor;
    system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[V_VIM]) = -Cim * Rv * Rtot;
  }

  // Eqn 2 (T): Pavg = average pressure across Ra1_myo = Pa + 0.5*Ra1_myo*Qin
  // (Ra1_myo is the last resistor before Ca; Pa is the pressure at Ca)
  auto Pavg = Pa + 0.5 * Ra1_myo * q_in;
  auto A    = std::pow(Kar_myo_ / Ra1_myo, 0.25);
  system.C(global_eqn_ids[2]) = -Pavg * A;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_PA])  = -A;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_QIN]) = -0.5 * Ra1_myo * A;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[V_AM])  =
      -0.5 * q_in * A * dRa1m_dAm + Pavg * 0.25 * A * dRa1m_dAm / Ra1_myo;

  // Eqn 3 (WSS): flow through Ra1_shear is Qin (before Ca)
  auto B = std::pow(Ra1_shear / Kar1_shear_, 0.75);
  system.C(global_eqn_ids[3]) = -q_in * B;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[V_QIN]) = -B;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[V_AS])  = -q_in * 0.75 * B * dRa1s_dAs / Ra1_shear;

  // Eqn 10 (Pa): nonlinear part (Ra1_shear + Ra1_myo)*Qin
  system.C(global_eqn_ids[10]) = (Ra1_shear + Ra1_myo) * q_in;
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[V_QIN]) =  Ra1_shear + Ra1_myo;
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[V_AS])  =  dRa1s_dAs * q_in;
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[V_AM])  =  dRa1m_dAm * q_in;
}
