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

namespace
{
std::string describe(const mfem::ParFiniteElementSpace & f)
{
  return std::string(f.FEColl()->Name()) + ", vdim " + std::to_string(f.GetVDim()) +
         ", ordering " + (f.GetOrdering() == mfem::Ordering::byNODES ? "NODES" : "VDIM");
}

void checkSameSpace(const mfem::ParGridFunction & moose, const mfem::ParGridFunction & lib,
                    const std::string & name)
{
  const auto & mf = *moose.ParFESpace();
  const auto & lf = *lib.ParFESpace();

  const bool same_fec   = std::string(mf.FEColl()->Name()) == lf.FEColl()->Name();
  const bool same_vdim  = mf.GetVDim() == lf.GetVDim();
  // Ordering only changes the layout when there is more than one component.
  const bool same_order = lf.GetVDim() == 1 || mf.GetOrdering() == lf.GetOrdering();
  const bool same_size  = moose.Size() == lib.Size();

  if (!(same_fec && same_vdim && same_order && same_size))
    mooseError("Variable '", name, "' has FE space [", describe(mf), ", size ", moose.Size(),
               "] but the solver needs [", describe(lf), ", size ", lib.Size(), "]");
}
} // namespace

void LmmhdProblemOperator::Init(mfem::BlockVector &)
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

  checkSameSpace(*_problem.getGridFunction("current_density"), _stepper->CurrentDensity(), "current_density");
  checkSameSpace(*_problem.getGridFunction("electric_potential"), _stepper->Potential(), "electric_potential");
  checkSameSpace(*_problem.getGridFunction("velocity"), _stepper->Velocity(), "velocity");
  checkSameSpace(*_problem.getGridFunction("pressure"), _stepper->Pressure(), "pressure");

  // Copy data back to MOOSE variables for first time-step.
  *_problem.getGridFunction("current_density") = _stepper->CurrentDensity();
  *_problem.getGridFunction("electric_potential") = _stepper->Potential();
  *_problem.getGridFunction("velocity") = _stepper->Velocity();
  *_problem.getGridFunction("pressure") = _stepper->Pressure();
}

void LmmhdProblemOperator::Solve()
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

  // Copy data back to MOOSE variables.
  *_problem.getGridFunction("current_density")    = _stepper->CurrentDensity();
  *_problem.getGridFunction("electric_potential") = _stepper->Potential();
  *_problem.getGridFunction("velocity")           = _stepper->Velocity();
  *_problem.getGridFunction("pressure")           = _stepper->Pressure();
}

#endif