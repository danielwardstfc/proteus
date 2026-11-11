[Mesh]
  type = MFEMFileMesh            # file-based MFEM mesh; name may differ in your MOOSE branch
  file = box_clustered.mesh
[]

[Problem]
  type = MFEMProblem
[]

[ProblemComposers]
  [lmmhd]
    type = LmmhdProblemComposer
    dt = 0.01
    Re = 1
    kappa = 1
    B = '1 0 0'
  []
[]

[Executioner]
  type = MFEMTransient
  dt = 0.01                  # must equal the composers dt
  end_time = 0.2
  device = cpu
[]

[Outputs]
[]