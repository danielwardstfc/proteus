//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "LmmhdProblemComposer.h"

registerMooseObject("LmmhdMooseApp", LmmhdProblemComposer);

InputParameters
LmmhdProblemComposer::validParams()
{
  InputParameters params = MFEMProblemComposer::validParams();
  params.addClassDescription("Liquid-metal MHD solver (liblmmhd) driven by MFEMTransient.");

  params.addRequiredParam<Real>("dt", "Time step; must equal the Executioner's constant dt.");
  params.addParam<Real>("theta", 0.6, "Theta-scheme parameter.");

  params.addRequiredParam<Real>("Re", "Reynolds number.");
  params.addRequiredParam<Real>("kappa", "Interaction parameter N = Ha^2 / Re.");
  params.addParam<Real>("alpha", 1.0, "Grad-div parameter.");
  params.addParam<std::vector<Real>>("B", {0.0, 0.0, 1.0}, "Applied magnetic field direction.");

  params.addParam<bool>("use_ads", true, "Current-density block: CG+ADS (true) or MUMPS (false).");

  params.addParam<VariableName>("current_density", "current_density", "MOOSE variable that receives j.");
  params.addParam<VariableName>("electric_potential", "velocelectric_potentiality", "MOOSE variable that receives phi.");
  params.addParam<VariableName>("velocity", "velocity", "MOOSE variable that receives u.");
  params.addParam<VariableName>("pressure", "pressure", "MOOSE variable that receives p.");
  return params;
}

LmmhdProblemComposer::LmmhdProblemComposer(const InputParameters & parameters)
  : MFEMProblemComposer(parameters)
{
  // Map MOOSE input parameters onto the library's Parameters struct.
  // Adjust member names to match your parameters.hpp.
  _lmmhd_params.time.dt = getParam<Real>("dt");
  _lmmhd_params.time.theta = getParam<Real>("theta");
  _lmmhd_params.physics.Re = getParam<Real>("Re");
  _lmmhd_params.physics.kappa = getParam<Real>("kappa");
  _lmmhd_params.physics.alpha = getParam<Real>("alpha");

  const auto & B = getParam<std::vector<Real>>("B");
  if (B.size() != 3)
    paramError("B", "must have 3 components");
  _lmmhd_params.physics.B = {B[0], B[1], B[2]};

  _lmmhd_params.solver.use_ads_for_current_density = getParam<bool>("use_ads");
}

std::shared_ptr<Moose::MFEM::ProblemOperatorBase>
LmmhdProblemComposer::createProblemOperator(MFEMProblem & mfem_problem)
{
  return std::make_shared<LmmhdProblemOperator>(
      mfem_problem, _lmmhd_params);
}

#endif