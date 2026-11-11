//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "LmmhdProblemOperator.h"
#include "MFEMProblem.h"

#include "BoundaryConditions.hpp"

LmmhdProblemOperator::LmmhdProblemOperator(MFEMProblem & mfem_problem,
                                           const lmmhd::Parameters & params)
  : Moose::MFEM::TimeDependentProblemOperator(mfem_problem), _params(params)
{
}

void
LmmhdProblemOperator::Init(mfem::BlockVector &)
{
  mfem::ParMesh & pmesh = _problem.mesh().getMFEMParMesh();

  if (mfem::Mpi::Root())
    _log_file.open("simulation.log");
  _logger = std::make_unique<Logger>(_log_file);

  _disc = std::make_unique<Discretisation>(pmesh, _params.discretisation);

  mfem::Vector bb_min, bb_max;
  pmesh.GetBoundingBox(bb_min, bb_max);
  auto bcs = LidDrivenCavity(pmesh, bb_max(2));

  _stepper = std::make_unique<LmmhdTimeStepper>(
      _disc->Spaces(), std::move(bcs), _params, *_logger);
}

void
LmmhdProblemOperator::Solve()
{
  // The operator builds 1/(theta dt) into its matrices, so MOOSE must use the
  // same, fixed time step.
  const mfem::real_t dt = _problem.dt();
  if (std::abs(dt - _params.time.dt) > 1e-12 * _params.time.dt)
    mooseError("LmmhdProblemOperator: MOOSE time step (", dt,
               ") differs from the composer's dt (", _params.time.dt,
               "). Use a constant dt equal to the composer's.");

  *_logger << "Time step " << _stepper->StepIndex() << ", time = " << _stepper->Time()
           << ", time elapsed = " << mfem::toc() << std::endl;

  const mfem::real_t residual = _stepper->SolveTimeStep();
}

#endif