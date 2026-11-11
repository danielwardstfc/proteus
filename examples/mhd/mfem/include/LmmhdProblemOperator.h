//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "TimeDependentProblemOperator.h"

#include <fstream>
#include <memory>

// From liblmmhd (contrib/lmmhd/include).  These are .hpp, MOOSE's are .h,
// so the flat include path can't pick up a MOOSE header by mistake.
#include "parameters.hpp"
#include "tools.hpp"
#include "Discretisation.hpp"
#include "LmmhdTimeStepper.hpp"

/*
  Thin MOOSE wrapper around lmmhd::LmmhdTimeStepper.

  Init()  : build the FE spaces on MOOSE's mesh, the boundary conditions and
            the time stepper.
  Solve() : called by MFEMTransient once per MOOSE time step -> one
            SolveTimeStep().

 */
class LmmhdProblemOperator : public Moose::MFEM::TimeDependentProblemOperator
{
public:
  LmmhdProblemOperator(MFEMProblem & mfem_problem,
                       const lmmhd::Parameters & params);

  void Init(mfem::BlockVector &) override;

  void Solve() override;

private:
  lmmhd::Parameters _params;

  std::ofstream _log_file;
  std::unique_ptr<Logger> _logger;
  std::unique_ptr<Discretisation> _disc;
  std::unique_ptr<LmmhdTimeStepper> _stepper;
};

#endif