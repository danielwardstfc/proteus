//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "LmmhdMooseTestApp.h"
#include "LmmhdMooseApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
LmmhdMooseTestApp::validParams()
{
  InputParameters params = LmmhdMooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

LmmhdMooseTestApp::LmmhdMooseTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  LmmhdMooseTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

LmmhdMooseTestApp::~LmmhdMooseTestApp() {}

void
LmmhdMooseTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  LmmhdMooseApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"LmmhdMooseTestApp"});
    Registry::registerActionsTo(af, {"LmmhdMooseTestApp"});
  }
}

void
LmmhdMooseTestApp::registerApps()
{
  registerApp(LmmhdMooseApp);
  registerApp(LmmhdMooseTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
LmmhdMooseTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LmmhdMooseTestApp::registerAll(f, af, s);
}
extern "C" void
LmmhdMooseTestApp__registerApps()
{
  LmmhdMooseTestApp::registerApps();
}
