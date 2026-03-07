// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file VarResistanceVessel.h
 * @brief model::VarResistanceVessel source file
 */
#ifndef SVZERODSOLVER_MODEL_VARRESISTANCEVESSEL_HPP_
#define SVZERODSOLVER_MODEL_VARRESISTANCEVESSEL_HPP_
#include <math.h>

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"
#include "Model.h"


class VarResistanceVessel : public Block {
 public:

  VarResistanceVessel(int id, Model *model)
      : Block(id, model, BlockType::var_resistance_vessel, BlockClass::vessel,
              {{"R", InputParameter()},
              {"A1", InputParameter()},
              {"t1", InputParameter()},
              {"k1", InputParameter()},
              {"A2", InputParameter()},
              {"t2", InputParameter()},
              {"k2", InputParameter()}}) {}


  void setup_dofs(DOFHandler &dofhandler);

  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  TripletsContributions num_triplets{5, 0, 1};
};

#endif  // SVZERODSOLVER_MODEL_VARRESISTANCEVESSEL_HPP_
