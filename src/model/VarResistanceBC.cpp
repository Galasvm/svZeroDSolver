// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "VarResistanceBC.h"

void VarResistanceBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void VarResistanceBC::update_constant(SparseSystem &system,
                                       std::vector<double> &parameters) {
  double R = parameters[global_param_ids[0]];
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -R;
}


void VarResistanceBC::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  
  
  auto glob_time = model->time;
  // Get states
  double q_in = y[global_var_ids[1]];

  // Get parameters
  double R = parameters[global_param_ids[0]];
  double Pd = parameters[global_param_ids[1]];
  double A1 = parameters[global_param_ids[2]];
  double t1 = parameters[global_param_ids[3]];
  double k1 = parameters[global_param_ids[4]];
  double A2 = parameters[global_param_ids[5]];
  double t2 = parameters[global_param_ids[6]];
  double k2 = parameters[global_param_ids[7]];

  // Nonlinear term
  system.C(global_eqn_ids[0]) =
      - Pd - q_in * (  A1 * R * (1 / (1 + exp(-(glob_time - t1)/k1))) + A2 * R * (1 / (1 + exp(-(glob_time - t2)/k2)))     );

  // Derivatives of non-linear term
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      - (A1 * R * (1 / (1 + exp(-(glob_time - t1)/k1))) + A2 * R * (1 / (1 + exp(-(glob_time - t2)/k2))));
}