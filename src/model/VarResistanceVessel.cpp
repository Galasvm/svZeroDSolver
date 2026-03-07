// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "VarResistanceVessel.h"

void VarResistanceVessel::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {});
}

void VarResistanceVessel::update_constant(SparseSystem &system,
                                  std::vector<double> &parameters) {
  // Get parameters
  double R = parameters[global_param_ids[0]];

  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -R;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}


void VarResistanceVessel::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {



  auto glob_time = model->time;
  // Get states
  double q_in = y[global_var_ids[1]];

  // Get parameters
  double R = parameters[global_param_ids[0]];
  double A1 = parameters[global_param_ids[1]];
  double t1 = parameters[global_param_ids[2]];
  double k1 = parameters[global_param_ids[3]];
  double A2 = parameters[global_param_ids[4]];
  double t2 = parameters[global_param_ids[5]];
  double k2 = parameters[global_param_ids[6]];

  // Nonlinear term
  system.C(global_eqn_ids[0]) =
     - q_in * (  A1 * R * (1 / (1 + exp(-(glob_time - t1)/k1))) + A2 * R * (1 / (1 + exp(-(glob_time - t2)/k2)))     );

  // Derivatives of non-linear term
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      - (A1 * R * (1 / (1 + exp(-(glob_time - t1)/k1))) + A2 * R * (1 / (1 + exp(-(glob_time - t2)/k2))));
}
