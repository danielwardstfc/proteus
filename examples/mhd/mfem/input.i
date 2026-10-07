[Mesh]
  type = MFEMMeshGeneratorMesh

  [mfem_mesh_generators]
      type = MFEMGeneratedMeshGenerator
      dim = 3

      nx = 12
      ny = 12
      nz = 12

      xmax = 1.0
      ymax = 1.0
      zmax = 1.0

      elem_type = HEX
  []
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
[]

[Variables]
[]

[Solvers]
[]

[ProblemComposers]
  [default_steady]
    type = LmmhdProblemComposer
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
[]

[Outputs]
[]
