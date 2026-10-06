// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "VarRCRBC.h"

void VarRCRBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {"pressure_c"});
}

void VarRCRBC::update_constant(SparseSystem &system,
                               std::vector<double> &parameters) {
  double Rp = parameters[global_param_ids[0]];

  // Eqn (0): Pin - Rp*Qin - Pc = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -Rp;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;

  // Eqn (1): Rd(t)*Qin - Pc + Pd - Rd(t)*C*dPc/dt = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -1.0;
}

void VarRCRBC::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  auto glob_time = model->time;

  // Get states
  double q_in = y[global_var_ids[1]];
  double dpc_dt = dy[global_var_ids[2]];

  // Get parameters
  double C = parameters[global_param_ids[1]];
  double Rd = parameters[global_param_ids[2]];
  double Pd = parameters[global_param_ids[3]];
  double A1 = parameters[global_param_ids[4]];
  double t1 = parameters[global_param_ids[5]];
  double k1 = parameters[global_param_ids[6]];
  double A2 = parameters[global_param_ids[7]];
  double t2 = parameters[global_param_ids[8]];
  double k2 = parameters[global_param_ids[9]];

  // Time-varying distal resistance: Rd is always present as a baseline;
  // A1/A2 are optional additive deviations, each independently zeroable.
  double Rd_t =
      Rd * (1.0 + A1 * (1 / (1 + exp(-(glob_time - t1) / k1))) +
           A2 * (1 / (1 + exp(-(glob_time - t2) / k2))));

  // Nonlinear term: Rd(t)*Qin - Rd(t)*C*dPc/dt + Pd
  system.C(global_eqn_ids[1]) = Rd_t * q_in - Rd_t * C * dpc_dt + Pd;

  // Derivatives of non-linear term
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1]) = Rd_t;
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[2]) =
      -Rd_t * C;
}
