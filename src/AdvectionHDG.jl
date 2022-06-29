function advection_hdg_time_step!(
     pn, um, un, model, dΩ, ∂K, d∂K, X, Y, dt, τ,
     assem=SparseMatrixAssembler(SparseMatrixCSC{Float64,Int},Vector{Float64},X,Y))

  # Time centered second order advection
  # References:
  #   Muralikrishnan, Tran, Bui-Thanh, JCP, 2020 vol. 367
  #   Kang, Giraldo, Bui-Thanh, JCP, 2020 vol. 401

  dth = 0.5*dt

  n  = get_cell_normal_vector(∂K)
  nₒ = get_cell_owner_normal_vector(∂K)

  # First stage
  b₁(q,m) = ∫(q*pn)dΩ - 
            ∫(m*0)d∂K

  a₁((p,l),(q,m)) = ∫(q*p - dt*(∇q⋅un)*p)dΩ + ∫(dt*q*p*(un⋅n + τ*n⋅n))d∂K -   # [q,p] block
                    ∫(dt*q*l*τ*n⋅n)d∂K +                                      # [q,l] block
                    ∫(m*p*(un⋅nₒ + τ*n⋅nₒ))d∂K -                              # [m,p] block
                    ∫(m*l*τ*n⋅nₒ)d∂K                                          # [m,l] block

  op₁  = HybridAffineFEOperator((u,v)->(a₁(u,v),b₁(v)), X, Y)
  Xh   = solve(op₁)
  ph,_ = Xh

  # Second stage
  b₂(q,m) = ∫(q*pn + dth*(∇q⋅un)*ph)dΩ - ∫(dth*q*ph*(un⋅n + τ*n⋅n))d∂K - 
           ∫(0.5*m*q*(un⋅nₒ + τ*n⋅nₒ))d∂K

  a₂((p,l),(q,m)) = ∫(q*p - dth*(∇q⋅un)*p)dΩ + ∫(dth*q*p*(un⋅n + τ*n⋅n))d∂K - # [q,p] block
                    ∫(dt*q*l*τ*n⋅n)d∂K +                                      # [q,l] block
                    ∫(0.5*m*p*(un⋅nₒ + τ*n⋅nₒ))d∂K -                          # [m,p] block
                    ∫(m*l*τ*n⋅nₒ)d∂K                                          # [m,l] block

  op₂  = HybridAffineFEOperator((u,v)->(a₂(u,v),b₂(v)), X, Y)
  Xm   = solve(op₂)
  pm,_ = Xm

  get_free_dof_values(pn) .= get_free_dof_values(pm)
end

function project_initial_conditions(dΩ, P, Q, p₀, U, V, u₀, mass_matrix_solver)
  # the tracer
  b₁(q)    = ∫(q*p₀)dΩ
  a₁(p,q)  = ∫(q*p)dΩ
  rhs₁     = assemble_vector(b₁, Q)
  L2MM     = assemble_matrix(a₁, P, Q)
  L2MMchol = numerical_setup(symbolic_setup(mass_matrix_solver,L2MM),L2MM)
  pn       = FEFunction(Q, copy(rhs₁))
  pnv      = get_free_dof_values(pn)

  solve!(pnv, L2MMchol, pnv)

  # the velocity field
  b₂(v)    = ∫(v⋅u₀)dΩ
  a₂(u,v)  = ∫(v⋅u)dΩ
  rhs₂     = assemble_vector(b₂, V)
  MM       = assemble_matrix(a₂, U, V)
  MMchol   = numerical_setup(symbolic_setup(mass_matrix_solver,MM),MM)
  un       = FEFunction(V, copy(rhs₂))
  unv      = get_free_dof_values(un)

  solve!(unv, MMchol, unv)

  pn, pnv, L2MM, L2MMchol, un
end

function advection_hdg(
        model, order, degree,
        u₀, p₀, dt, N;
        mass_matrix_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        write_diagnostics=true,
        write_diagnostics_freq=1,
        dump_diagnostics_on_screen=true,
        write_solution=false,
        write_solution_freq=N/10,
        output_dir="adv_eq_ncells_$(num_cells(model))_order_$(order)")

  # Forward integration of the advection equation
  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D},model)
  Γ = Triangulation(ReferenceFE{D-1},model)
  #Ω = Triangulation(model)
  #Γ = SkeletonTriangulation(model)
  ∂K = GridapHybrid.Skeleton(model)

  reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order;space=:P)
  reffeₚ = ReferenceFE(lagrangian,Float64,order;space=:P)
  reffeₗ = ReferenceFE(lagrangian,Float64,order;space=:P)

  # Define test FESpaces
  V = TestFESpace(Ω, reffeᵤ; conformity=:L2)
  Q = TestFESpace(Ω, reffeₚ; conformity=:L2)
  M = TestFESpace(Γ, reffeₗ; conformity=:L2)
  Y = MultiFieldFESpace([Q,M])

  U = TrialFESpace(V)
  P = TrialFESpace(Q)
  L = TrialFESpace(M)
  X = MultiFieldFESpace([P,L])

  dΩ  = Measure(Ω,degree)
  d∂K = Measure(∂K,degree)

  # HDG stabilisation parameter
  τ = 1.0

  # Project the initial conditions onto the trial spaces
  pn, pnv, L2MM, L2MMchol, un = project_initial_conditions(dΩ, P, Q, p₀, U, V, u₀, mass_matrix_solver)

  # Work array
  p_tmp = copy(pnv)

  # Initial states
  mul!(p_tmp, L2MM, pnv)
  p01 = sum(p_tmp)  # total mass
  p02 = p_tmp⋅pnv   # total entropy

  function run_simulation(pvd=nothing)
    diagnostics_file = joinpath(output_dir,"adv_diagnostics.csv")
    if (write_diagnostics)
      initialize_csv(diagnostics_file,"step", "mass", "entropy")
    end

    for istep in 1:N
      advection_hdg_time_step!(pn, un, un, model, dΩ, ∂K, d∂K, X, Y, dt, τ)

      if (write_diagnostics && write_diagnostics_freq>0 && mod(istep, write_diagnostics_freq) == 0)
        # compute mass and entropy conservation
        mul!(p_tmp, L2MM, pnv)
        pn1 = sum(p_tmp)
        pn2 = p_tmp⋅pnv
	pn1 = (pn1 - p01)/p01
	pn2 = (pn2 - p02)/p02
	if dump_diagnostics_on_screen
          @printf("%5d\t%14.9e\t%14.9e\n", istep, pn1, pn2)
	end
	append_to_csv(diagnostics_file; step=istep, mass=pn1, entropy=pn2)
      end
      if (write_solution && write_solution_freq>0 && mod(istep, write_solution_freq) == 0)
        pvd[dt*Float64(istep)] = new_vtk_step(Ω,joinpath(output_dir,"n=$(istep)"),["pn"=>pn,"un"=>un])
      end
    end
    pn
  end
  if (write_diagnostics || write_solution)
    rm(output_dir,force=true,recursive=true)
    mkdir(output_dir)
  end
  if (write_solution)
    pvdfile=joinpath(output_dir,"adv_eq_ncells_$(num_cells(model))_order_$(order)")
    paraview_collection(run_simulation,pvdfile)
  else
    run_simulation()
  end
end
