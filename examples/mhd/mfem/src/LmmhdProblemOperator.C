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

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

#include "lmmhd/LmmhdOperator.hpp"
#include "lmmhd/constants.hpp"
#include "lmmhd/InputParser.hpp"
#include "lmmhd/tools.hpp"

// The custom operator constructor
LmmhdProblemOperator::LmmhdProblemOperator(MFEMProblem & mfem_problem)
  : Moose::MFEM::ProblemOperator(mfem_problem)
{
}

void LmmhdProblemOperator::Init(mfem::BlockVector &)
{
  // Get the FE-space and Variable that were just built
  /*auto fes = _problem.getProblemData().fespaces.Get("H1");
  auto gridfunction = _problem.getGridFunction("u");

  // Boundary conditions
  fes->GetBoundaryTrueDofs(_boundary_dofs);

  // Build the linear form
  _b = new mfem::ParLinearForm(fes);
  _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(_coef));
  _b->Assemble();

  // Build the bilinear form
  _a = new mfem::ParBilinearForm(fes);
  _a->AddDomainIntegrator(new mfem::DiffusionIntegrator);
  _a->Assemble();

  // Form the linear system
  _a->FormLinearSystem(_boundary_dofs, *gridfunction, *_b, _problem_operator, _X, _B);*/
}

void LmmhdProblemOperator::Solve()
{
   int debug = 0;

   std::ofstream log_file("simulation.log");
   Logger logger(log_file);

   // Initialize MPI and HYPRE.
   //Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   mfem::tic();

   // Read in parameters from file.
   InputParser input("params.in");

   // Mesh parameters
   int nx = input.GetInt("nx");
   int ny = input.GetInt("ny");
   int nz = input.GetInt("nz");

   real_t Lx = input.GetReal("Lx");
   real_t Ly = input.GetReal("Ly");
   real_t Lz = input.GetReal("Lz");

   real_t clusterX = input.GetReal("clusterX");
   real_t clusterY = input.GetReal("clusterY");
   real_t clusterZ = input.GetReal("clusterZ");

   // Timestepping
   real_t t_final = input.GetReal("t_final");
   real_t dt = input.GetReal("dt");
   int vis_steps = input.GetInt("vis_steps");

   // Dimensionless parameters
   Re = input.GetReal("Re");
   Ha = input.GetReal("Ha");

   Bx = input.GetReal("Bx");
   By = input.GetReal("By");
   Bz = input.GetReal("Bz");

   // Update constants.
   reciprocal_Re = 1.0 / Re;
   kappa_val = Ha * Ha / Re;
   neg_kappa_val = -kappa_val;
   alpha1 = alpha + reciprocal_Re;
   neg_alpha1 = -alpha1;


   if (Mpi::Root()) PrintParams(nx, ny, nz, Lx, Ly, Lz, clusterX, clusterY, clusterZ, t_final, dt, vis_steps, Re, Ha);

   // Command line options.
   /*OptionsParser args(argc, argv);
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&nx, "-nx", "--num-elements-x",
                  "Number of elements in the x direction.");
   args.AddOption(&ny, "-ny", "--num-elements-y",
                  "Number of elements in the y direction.");
   args.AddOption(&nz, "-nz", "--num-elements-z",
                  "Number of elements in the z direction.");
   args.AddOption(&Lx, "-lx", "--length-x",
                  "Length of the domain in the x direction.");
   args.AddOption(&Ly, "-ly", "--length-y",
                  "Length of the domain in the y direction.");
   args.AddOption(&Lz, "-lz", "--length-z",
                  "Length of the domain in the z direction.");
   args.AddOption(&clusterX, "-cx", "--cluster-x", "Clustering intensity in x direction.");
   args.AddOption(&clusterY, "-cy", "--cluster-y", "Clustering intensity in y direction.");
   args.AddOption(&clusterZ, "-cz", "--cluster-z", "Clustering intensity in z direction.");
   
   args.Parse();*/

   // Set fe_space orders.
   int order_currentD = 1;
   int order_electPot = 0;
   int order_pressure = 1;
   int order_velocity;
   order_velocity = order_pressure + 1;

   // Generate mesh.
   mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, ny, nz, mfem::Element::Type::HEXAHEDRON, Lx, Ly, Lz);
   //Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, mfem::Element::Type::TETRAHEDRON, Lx, Ly, Lz);
   //const char *mesh_file = "./mesh/cuboid_clustered.msh";
   //Mesh mesh = Mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   
   // Cluster vertices in y and z.
   int numVertices = mesh.GetNV();
   for (int i = 0; i < numVertices; i++)
   {
      real_t *v = mesh.GetVertex(i);
      real_t xi = v[0] / Lx; // Normalize y
      real_t yi = v[1] / Ly; // Normalize y
      real_t zi = v[2] / Lz; // Normalize z

      v[0] = cluster_symmetric(xi, clusterX) * Lx;
      v[1] = cluster_symmetric(yi, clusterY) * Ly;
      v[2] = cluster_symmetric(zi, clusterZ) * Lz;
   }

   //ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
   mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();
   

   // ----------------------------------------------------------------------------
   // Finite Element Spaces.
   // ----------------------------------------------------------------------------
   // HDiv finite elements for current density.
   RT_FECollection currentD_fec(order_currentD, dim);
   ParFiniteElementSpace currentD_fespace(&pmesh, &currentD_fec);

   // L2 finite elements for electric scalar potential.
   L2_FECollection electPot_fec(order_electPot, dim);
   ParFiniteElementSpace electPot_fespace(&pmesh, &electPot_fec);

   // H1 continuous Lagrange finite elements (order_pressure + 1)
   // for velocity.
   H1_FECollection velocity_fec(order_velocity, dim);
   ParFiniteElementSpace velocity_fespace(&pmesh, &velocity_fec, dim);

   // H1 continuous Lagrange finite elements for pressure.
   H1_FECollection pressure_fec(order_pressure, dim);
   ParFiniteElementSpace pressure_fespace(&pmesh, &pressure_fec);


   Array<ParFiniteElementSpace *> spaces(4);
   spaces[0] = &currentD_fespace;
   spaces[1] = &electPot_fespace;
   spaces[2] = &velocity_fespace;
   spaces[3] = &pressure_fespace;

   int j_space_size = currentD_fespace.GetTrueVSize();
   int phi_space_size = electPot_fespace.GetTrueVSize();
   int v_space_size = velocity_fespace.GetTrueVSize();
   int p_space_size = pressure_fespace.GetTrueVSize();

   // ----------------------------------------------------------------------------
   // Define boundaries.
   // ----------------------------------------------------------------------------
   
   // Essential (Dirichlet) boundary conditions for duct flow.
   //                  {z1, y0, x1 (outlet), y1, x0 (inlet), z0}
   //int currentD_bcs[] = {1,  1,  1,           1,  1,          1};
   //int electPot_bcs[] = {0,  0,  0,           0,  0,          0};
   //int pressure_bcs[] = {0,  0,  1,           0,  0,          0};
   //int velocity_bcs[] = {1,  1,  0,           1,  1,          1};

   // Essential (Dirichlet) boundary conditions for lid-driven cavity.
   //                  {z1, y0, x1, y1, x0, z0}
   int currentD_bcs[] = {1,  1,  1,  1,  1,  1};
   int electPot_bcs[] = {0,  0,  0,  0,  0,  0};
   int pressure_bcs[] = {0,  0,  0,  0,  0,  0};
   int velocity_bcs[] = {1,  1,  1,  1,  1,  1};

   Array<int> ess_boundary_marker_currentD(currentD_bcs, 6);
   Array<int> ess_boundary_marker_electPot(electPot_bcs, 6);
   Array<int> ess_boundary_marker_pressure(pressure_bcs, 6);
   Array<int> ess_boundary_marker_velocity(velocity_bcs, 6);

   Array<Array<int> *> ess_bdr(4);
   ess_bdr[0] = &ess_boundary_marker_currentD;
   ess_bdr[1] = &ess_boundary_marker_electPot;
   ess_bdr[2] = &ess_boundary_marker_velocity;
   ess_bdr[3] = &ess_boundary_marker_pressure;

   Array<int> pressure_ess_tdof, velocity_ess_tdof, currentD_ess_tdof, electPot_ess_tdof;
   currentD_fespace.GetEssentialTrueDofs(ess_boundary_marker_currentD, currentD_ess_tdof);
   electPot_fespace.GetEssentialTrueDofs(ess_boundary_marker_electPot, electPot_ess_tdof);
   pressure_fespace.GetEssentialTrueDofs(ess_boundary_marker_pressure, pressure_ess_tdof);
   velocity_fespace.GetEssentialTrueDofs(ess_boundary_marker_velocity, velocity_ess_tdof);

   // Print mesh statistics.
   if (Mpi::Root()) PrintFESpaces(j_space_size, phi_space_size, v_space_size, p_space_size);
   
   ParGridFunction uN_1(&velocity_fespace);
   GridFunction du(&velocity_fespace);

   // Define block structure of the solution vector (u then p).
   Array<int> block_trueOffsets(5);
   block_trueOffsets[0] = 0;
   block_trueOffsets[1] = j_space_size;
   block_trueOffsets[2] = phi_space_size;
   block_trueOffsets[3] = v_space_size;
   block_trueOffsets[4] = p_space_size;
   block_trueOffsets.PartialSum();

   // Initialise time-loop details.
   real_t t = 0.0;
   int ti_out = 0; // Time step output index.
   int ti = 0;

   // Initialise liquid-metal MHD operator.
   LmmhdOperator oper(spaces, ess_bdr, block_trueOffsets, dim, dt, debug, logger);

   ParGridFunction *j_gf = oper.GetCurrentDPointer();
   ParGridFunction *phi_gf = oper.GetPotentialPointer();
   ParGridFunction *u_gf = oper.GetVelocityPointer();
   ParGridFunction *p_gf = oper.GetPressurePointer();

   oper.SetGridFunctionsFromTrueDofs();

   // Set up visualisation in Paraview.
   ParaViewDataCollection pvdc("lmmhd", &pmesh);
   pvdc.SetPrefixPath("data");
   pvdc.SetDataFormat(VTKFormat::BINARY);
   pvdc.SetHighOrderOutput(true);
   pvdc.SetLevelsOfDetail(2);
   pvdc.SetCycle(0);
   pvdc.SetTime(t);
   pvdc.RegisterField("current density", j_gf);
   pvdc.RegisterField("electric potential", phi_gf);
   pvdc.RegisterField("velocity", u_gf);
   pvdc.RegisterField("pressure", p_gf);
   pvdc.Save();

   while (t < t_final)
   {
      if (Mpi::Root()) { std::cout << "Time step " << ti << ", time = " << t << ", time elapsed = " << mfem::toc() << std::endl; }
      logger << "Time step " << ti << ", time = " << t << ", time elapsed = " << mfem::toc() << std::endl;
            
      oper.UpdateUStar(ti); // Update value of u* in convection integrators each time-step/Picard iteration.
      oper.UpdateIntegrators(); // Propagate updated u* into F integrators.

      oper.Step(t, dt);

      oper.RemoveMeans();

      oper.SetGridFunctionsFromTrueDofs();

      oper.ReconstructPhysicalVelocityFromUBar(ti);

      oper.CalcNorms();

      oper.UpdateHistory();

      // Save data.
      pvdc.SetCycle(ti);
      pvdc.SetTime(t);
      pvdc.Save();

      ti += 1;
      t += dt;

      //if (res < tolerance) break;
   }

   if (Mpi::Root()) { std::cout << "Total simulation time: " << mfem::toc() << std::endl; }

   return;
}

LmmhdProblemOperator::~LmmhdProblemOperator()
{
  //delete _a;
  //delete _b;
}

#endif
