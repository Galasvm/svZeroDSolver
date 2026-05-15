// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_AUTOREGULATIONRCR2_HPP_
#define SVZERODSOLVER_MODEL_AUTOREGULATIONRCR2_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief AutoregulationRCR2: same as AutoregulationRCR but without q_out as
 * a DOF. q_out = Qin - C*dPc/dt is a derived quantity used internally for
 * Pavg, WSS, and the metabolic stimulus. Use AutoregulationRCR for debugging.
 *
 * Local unknown ordering:
 *   y^e = [ Pin, Qin, Ashear, Amyo, Ameta, xshear, xmyo, xmeta, T, WSS, Pc ]
 */
class AutoregulationRCR2 : public Block {
 public:
  AutoregulationRCR2(int id, Model *model)
      : Block(id, model,
              BlockType::autoregulation_rcr2,
              BlockClass::boundary_condition,
              {{"Rd", InputParameter()},
               {"Qt", InputParameter()},
               {"Pt", InputParameter()},
               {"Gshear", InputParameter()},
               {"taushear", InputParameter()},
               {"Gmyo", InputParameter()},
               {"taumyo", InputParameter()},
               {"Gmeta", InputParameter()},
               {"taumeta", InputParameter()},
               {"Pd", InputParameter()},
               {"Rp", InputParameter()},
               {"C", InputParameter()},
               {"lower_frac", InputParameter(true, false, true, 0.70)},
               {"upper_frac", InputParameter(true, false, true, 1.30)}}) {}

  void setup_dofs(DOFHandler &dofhandler);

  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  // F=15, E=7, D=12
  TripletsContributions num_triplets{15, 7, 12};

 private:
  bool initialized_ = false;

  double R1L_ = 0.0, R1U_ = 0.0;
  double R2L_ = 0.0, R2U_ = 0.0;
  double R3L_ = 0.0, R3U_ = 0.0;
  double R4_  = 0.0;

  double Kar1_ = 0.0;
  double Kar2_ = 0.0;
  double WSSt_ = 0.0;
  double Tt_   = 0.0;
};

#endif  // SVZERODSOLVER_MODEL_AUTOREGULATIONRCR2_HPP_
