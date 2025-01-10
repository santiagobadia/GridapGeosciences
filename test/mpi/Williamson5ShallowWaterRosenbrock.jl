module Williamson5ShallowWaterRosenbrock

  using LinearMaps
  using SparseMatricesCSR
  using PartitionedArrays
  using Test
  using FillArrays
  using Gridap
  using GridapPETSc
  using GridapGeosciences
  using GridapDistributed
  using GridapSolvers
  using GridapP4est

  include("Williamson5InitialConditions.jl")

  function petsc_gamg_options()
    """
      -ksp_type gmres -ksp_rtol 1.0e-06 -ksp_atol 0.0
      -ksp_monitor -pc_type gamg -pc_gamg_type agg
      -mg_levels_esteig_ksp_type gmres -mg_coarse_sub_pc_type lu
      -mg_coarse_sub_pc_factor_mat_ordering_type nd -pc_gamg_process_eq_limit 50
      -pc_gamg_square_graph 9 pc_gamg_agg_nsmooths 1
    """
  end
  function petsc_mumps_options()
    """
      -ksp_type preonly -ksp_error_if_not_converged true
      -pc_type lu -pc_factor_mat_solver_type mumps
    """
  end

  order  = 0
  degree = 4

  # magnitude of the descent direction of the implicit solve;
  # neutrally stable for 0.5, L-stable for 1+sqrt(2)/2
  λ = 1.0 + 0.5*sqrt(2.0)

  num_refinements = 1

  dt     = 600.0/(2^num_refinements)
  τ      = 0.5*dt
  nstep  = Int(24*60^2*20/dt) # 20 days

  function main(distribute,parts)
    ranks = distribute(LinearIndices((prod(parts),)))

    # Change directory to the location of the script, where the mesh data files are located 
    cd(@__DIR__)

    coarse_model, cell_panels, coarse_cell_wise_vertex_coordinates = parse_cubed_sphere_coarse_model("williamson-5-C12/connectivity-gridapgeo.txt",
                                                                                                     "williamson-5-C12/geometry-gridapgeo.txt")


    #GridapPETSc.with(args=split(petsc_mumps_options())) do
      model = CubedSphereDiscreteModel(ranks,
                                       coarse_model,
                                       coarse_cell_wise_vertex_coordinates,
                                       cell_panels,
                                       num_refinements;
                                       radius=Rₑ,
                                       #adaptive=true,
                                       adaptive=false,
                                       order=1)
      #model = CubedSphereDiscreteModel(ranks, 12; radius=Rₑ)

      P          = JacobiLinearSolver()
      mm_solver  = GridapSolvers.CGSolver(P;rtol=1.e-6)
      jac_solver = GMRESSolver(100;Pr=nothing,Pl=nothing,maxiter=2000,atol=1e-12,rtol=1.e-6,restart=true,m_add=20,verbose=false,name="GMRES")

      hf, uf = shallow_water_rosenbrock_time_stepper(model, order, degree,
                                                     h₀, u₀, Ωₑ, gₑ, H₀,
                                                     λ, dt, τ, nstep,
                                                     #mm_solver, jac_solver;
                                                     BackslashSolver(), BackslashSolver();
                                                     t₀=topography,
                                                     leap_frog=true,
                                                     write_solution=true,
                                                     write_solution_freq=72,
                                                     write_diagnostics=true,
                                                     write_diagnostics_freq=1,
                                                     dump_diagnostics_on_screen=true)
    #end
  end

  with_mpi() do distribute 
    main(distribute,1)
  end
end # module
