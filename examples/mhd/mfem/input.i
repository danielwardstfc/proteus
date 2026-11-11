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

[FESpaces]
  [RT0]
    type = MFEMGenericFESpace
    fec_name = RT_3D_P0
  []
  [L2_0]
    type = MFEMGenericFESpace
    fec_name = L2_3D_P0
  []
  [H1_2_vec]
    type = MFEMGenericFESpace
    fec_name = H1_3D_P2
    vdim = 3
    ordering = NODES
  []
  [H1_1]
    type = MFEMGenericFESpace
    fec_name = H1_3D_P1
  []
[]

[AuxVariables]
  [current_density]
    type = MFEMVariable
    fespace = RT0
  []
  [electric_potential]
    type = MFEMVariable
    fespace = L2_0
  []
  [velocity]
    type = MFEMVariable
    fespace = H1_2_vec
  []
  [pressure]
    type = MFEMVariable
    fespace = H1_1
  []
[]

[Executioner]
  type = MFEMTransient
  dt = 0.01                  # must equal the composers dt
  end_time = 0.2
  device = cpu
[]

[Outputs]
  [ParaView]
    type = MFEMParaViewDataCollection
    file_base = OutputData/lmmhd
  []
[]