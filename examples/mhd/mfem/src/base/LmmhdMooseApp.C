#include "LmmhdMooseApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
LmmhdMooseApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

LmmhdMooseApp::LmmhdMooseApp(const InputParameters & parameters) : MooseApp(parameters)
{
  LmmhdMooseApp::registerAll(_factory, _action_factory, _syntax);
}

LmmhdMooseApp::~LmmhdMooseApp() {}

void
LmmhdMooseApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<LmmhdMooseApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"LmmhdMooseApp"});
  Registry::registerActionsTo(af, {"LmmhdMooseApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
LmmhdMooseApp::registerApps()
{
  registerApp(LmmhdMooseApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
LmmhdMooseApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LmmhdMooseApp::registerAll(f, af, s);
}
extern "C" void
LmmhdMooseApp__registerApps()
{
  LmmhdMooseApp::registerApps();
}
