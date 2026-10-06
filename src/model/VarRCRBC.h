// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file VarRCRBC.h
 * @brief model::VarRCRBC source file
 */
#ifndef SVZERODSOLVER_MODEL_VARRCRBC_HPP_
#define SVZERODSOLVER_MODEL_VARRCRBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"
#include "Model.h"

/**
 * @brief RCR Windkessel boundary condition with a time-varying distal
 * resistance.
 *
 * Identical circuit topology to \ref WindkesselBC (fixed proximal resistance
 * Rp and capacitance C), except the distal resistance is modulated in time
 * by the two-sigmoid envelope used in \ref VarResistanceBC. Rd is always
 * present as the baseline value; A1 and A2 are independent, optional
 * additive deviations (each defaults to a no-op when set to 0, regardless
 * of its paired t/k):
 *
 * \f[
 * R_{d}(t) = R_{d} \left( 1 + \frac{A_1}{1+e^{-(t-t_1)/k_1}} +
 * \frac{A_2}{1+e^{-(t-t_2)/k_2}} \right)
 * \f]
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_p$, *-] (3,0)
 * node[anchor=south]{$P_{c}$}
 * to [R, l=$R_d(t)$, *-*] (5,0)
 * node[anchor=south]{$P_{d}$}
 * (3,0) to [C, l=$C$, *-] (3,-1.5)
 * node[ground]{};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{in}-P_{c}-R_{p} Q_{in}=0
 * \f]
 *
 * \f[
 * R_{d}(t) Q_{in}-P_{c}+P_{d}-R_{d}(t) C \frac{d P_{c}}{d t}=0
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Proximal resistance (Rp)
 * * `1` Capacitance (C)
 * * `2` Distal resistance (Rd), baseline value modulated by the envelope
 * * `3` Distal pressure (Pd)
 * * `4` A1
 * * `5` t1
 * * `6` k1
 * * `7` A2
 * * `8` t2
 * * `9` k2
 *
 * ### Usage in json configuration file
 *
 *     "boundary_conditions": [
 *         {
 *             "bc_name": "OUT",
 *             "bc_type": "VarRCRBC",
 *             "bc_values": {
 *                 "Rp": 1000.0,
 *                 "C": 0.0001,
 *                 "Rd": 1000.0,
 *                 "Pd": 0.0,
 *                 "A1": 1.0,
 *                 "t1": 1.0,
 *                 "k1": 0.1,
 *                 "A2": 0.0,
 *                 "t2": 2.0,
 *                 "k2": 0.1
 *             }
 *         }
 *     ]
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `pressure_c`: Pressure at the capacitor
 *
 */
class VarRCRBC : public Block {
 public:
  VarRCRBC(int id, Model *model)
      : Block(id, model, BlockType::var_rcr_bc,
              BlockClass::boundary_condition,
              {{"Rp", InputParameter()},
               {"C", InputParameter()},
               {"Rd", InputParameter()},
               {"Pd", InputParameter(true)},
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

  TripletsContributions num_triplets{5, 0, 2};
};

#endif  // SVZERODSOLVER_MODEL_VARRCRBC_HPP_
