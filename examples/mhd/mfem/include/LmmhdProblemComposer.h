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

#include "MFEMProblemComposer.h"
#include "LmmhdProblemOperator.h"

/**
 * Reads the LMMHD parameters from the input file and creates the
 * LmmhdProblemOperator that MFEMTransient will drive.
 */
class LmmhdProblemComposer : public MFEMProblemComposer
{
public:
  static InputParameters validParams();

  LmmhdProblemComposer(const InputParameters & parameters);

  std::shared_ptr<Moose::MFEM::ProblemOperatorBase>
  createProblemOperator(MFEMProblem & mfem_problem) override;

private:
  lmmhd::Parameters _lmmhd_params;
};

#endif