function compute_potential_vorticity!(q,H1h,H1hchol,dΩ,R,S,h,u,f,n,assem,bwrk)
  a(r,s) = ∫(s*h*r)dΩ
  c(s)   = ∫(perp(n,∇(s))⋅(u) + s*f)dΩ
  Gridap.FESpaces.assemble_matrix_and_vector!(a, c, H1h, bwrk, assem, R, S)
  numerical_setup!(H1hchol,H1h)
  solve!(get_free_dof_values(q),H1hchol,bwrk)
end

function compute_velocity!(u1,dΩ,dω,V,RTMMchol,u2,qAPVM,F,ϕ,n,dt1,dt2,bwrk)
  b(v) = ∫(v⋅u2 - dt1*(qAPVM)*(v⋅⟂(F,n)))dΩ + ∫(dt2*DIV(v)*ϕ)dω
  Gridap.FESpaces.assemble_vector!(b, bwrk, V)
  solve!(get_free_dof_values(u1),RTMMchol,bwrk)
end

function compute_mass_flux!(F,dΩ,V,RTMMchol,u,bwrk)
  b(v) = ∫(v⋅u)dΩ
  Gridap.FESpaces.assemble_vector!(b, bwrk, V)
  solve!(get_free_dof_values(F),RTMMchol,bwrk)
end

function compute_depth!(h1,dΩ,dω,Q,L2MMchol,h2,F,dt,bwrk)
  b(q)  = ∫(q*h2)dΩ - ∫(dt*q*DIV(F))dω
  Gridap.FESpaces.assemble_vector!(b, bwrk, Q)
  solve!(get_free_dof_values(h1),L2MMchol,bwrk)
end

function compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMchol,uu,h,g,bwrk)
  b(q)  = ∫(q*(0.5*uu + g*h))*dΩ
  Gridap.FESpaces.assemble_vector!(b, bwrk, Q)
  solve!(get_free_dof_values(ϕ),L2MMchol,bwrk)
end

function compute_diagnostic_vorticity!(w,dΩ,S,H1MMchol,u,n,bwrk)
  b(s) = ∫(⟂(n,∇(s))⋅(u))dΩ
  Gridap.FESpaces.assemble_vector!(b, bwrk, S)
  solve!(get_free_dof_values(w),H1MMchol,bwrk)
end

function shallow_water_explicit_time_step!(
     h₂, u₂, hₚ, uₚ, ϕ, F, q₁, q₂, H1h, H1hchol,      # in/out args
     model, dΩ, dω, V, Q, R, S, f, g, h₁, u₁, hₘ, uₘ, # in args
     RTMMchol, L2MMchol, dt, τ, leap_frog, h_wrk, u_wrk, w_wrk,
     assem=SparseMatrixAssembler(SparseMatrixCSC{Float64,Int},Vector{Float64},R,S))            # more in args

  # energetically balanced explicit second order shallow water solver
  # reference: eqns (21-24) of
  # https://github.com/BOM-Monash-Collaborations/articles/blob/main/energetically_balanced_time_integration/EnergeticallyBalancedTimeIntegration_SW.tex

  n = get_normal_vector(Triangulation(model))
  # explicit step for provisional velocity, uₚ
  dt1 = dt
  if leap_frog
    dt1 = 2.0*dt
  end

  # 1.1: the mass flux
  compute_mass_flux!(F,dΩ,V,RTMMchol,u₁*h₁,u_wrk)
  # 1.2: the bernoulli function
  compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMchol,u₁⋅u₁,h₁,g,h_wrk)
  # 1.3: the potential vorticity
  compute_potential_vorticity!(q₁,H1h,H1hchol,dΩ,R,S,h₁,u₁,f,n,assem,w_wrk)
  # 1.4: solve for the provisional velocity
  compute_velocity!(uₚ,dΩ,dω,V,RTMMchol,uₘ,q₁-τ*u₁⋅∇(q₁),F,ϕ,n,dt1,dt1,u_wrk)
  # 1.5: solve for the provisional depth
  compute_depth!(hₚ,dΩ,dω,Q,L2MMchol,hₘ,F,dt1,h_wrk)

  # 2.1: the mass flux
  compute_mass_flux!(F,dΩ,V,RTMMchol,u₁*(2.0*h₁ + hₚ)/6.0+uₚ*(h₁ + 2.0*hₚ)/6.0,u_wrk)
  # 2.2: the bernoulli function
  compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMchol,(u₁⋅u₁ + u₁⋅uₚ + uₚ⋅uₚ)/3.0,0.5*(h₁ + hₚ),g,h_wrk)
  # 2.3: the potential vorticity
  compute_potential_vorticity!(q₂,H1h,H1hchol,dΩ,R,S,hₚ,uₚ,f,n,assem,w_wrk)
  # 2.4: solve for the final velocity
  compute_velocity!(u₂,dΩ,dω,V,RTMMchol,u₁,q₁-τ*u₁⋅∇(q₁)+q₂-τ*uₚ⋅∇(q₂),F,ϕ,n,0.5*dt,dt,u_wrk)
  # 2.5: solve for the final depth
  compute_depth!(h₂,dΩ,dω,Q,L2MMchol,h₁,F,dt,h_wrk)
end

function shallow_water_ssprk3_time_step!(
     h, u, ϕ, F, q, H1h, H1hchol, model, dΩ, dω, 
     V, Q, R, S, f, g, hₙ, uₙ, hₚ, uₚ,
     RTMMchol, L2MMchol, dt, τ, h_wrk, u_wrk, w_wrk,
     assem=SparseMatrixAssembler(SparseMatrixCSC{Float64,Int},Vector{Float64},R,S))

  n = get_normal_vector(Triangulation(model))
  get_free_dof_values(uₙ) .= get_free_dof_values(u)
  get_free_dof_values(hₙ) .= get_free_dof_values(h)

  # stage 1 - diagnostics
  compute_mass_flux!(F,dΩ,V,RTMMchol,uₙ*hₙ,u_wrk)
  compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMchol,uₙ⋅uₙ,hₙ,g,h_wrk)
  compute_potential_vorticity!(q,H1h,H1hchol,dΩ,R,S,hₙ,uₙ,f,n,assem,w_wrk)
  # stage 1 - prognostics
  compute_velocity!(uₚ,dΩ,dω,V,RTMMchol,uₙ,q-0.0*τ*uₙ⋅∇(q),F,ϕ,n,dt,dt,u_wrk)
  compute_depth!(hₚ,dΩ,dω,Q,L2MMchol,hₙ,F,dt,h_wrk)

  # stage 2 - diagnostics
  compute_mass_flux!(F,dΩ,V,RTMMchol,uₚ*hₚ,u_wrk)
  compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMchol,uₚ⋅uₚ,hₚ,g,h_wrk)
  compute_potential_vorticity!(q,H1h,H1hchol,dΩ,R,S,hₚ,uₚ,f,n,assem,w_wrk)
  # stage 2 - prognostics
  _dt = 0.25*dt
  get_free_dof_values(uₚ) .= 0.25*get_free_dof_values(uₚ) .+ 0.75*get_free_dof_values(uₙ)
  get_free_dof_values(hₚ) .= 0.25*get_free_dof_values(hₚ) .+ 0.75*get_free_dof_values(hₙ)
  compute_velocity!(uₚ,dΩ,dω,V,RTMMchol,uₚ,q-0.0*τ*uₙ⋅∇(q),F,ϕ,n,_dt,_dt,u_wrk)
  compute_depth!(hₚ,dΩ,dω,Q,L2MMchol,hₚ,F,_dt,h_wrk)

  # stage 3 - diagnostics
  compute_mass_flux!(F,dΩ,V,RTMMchol,uₚ*hₚ,u_wrk)
  compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMchol,uₚ⋅uₚ,hₚ,g,h_wrk)
  compute_potential_vorticity!(q,H1h,H1hchol,dΩ,R,S,hₚ,uₚ,f,n,assem,w_wrk)
  # stage 3 - prognostics
  _dt = (2.0/3.0)*dt
  get_free_dof_values(uₚ) .= (2.0/3.0)*get_free_dof_values(uₚ) .+ (1.0/3.0)*get_free_dof_values(uₙ)
  get_free_dof_values(hₚ) .= (2.0/3.0)*get_free_dof_values(hₚ) .+ (1.0/3.0)*get_free_dof_values(hₙ)
  compute_velocity!(u,dΩ,dω,V,RTMMchol,uₚ,q-0.0*τ*uₙ⋅∇(q),F,ϕ,n,_dt,_dt,u_wrk)
  compute_depth!(h,dΩ,dω,Q,L2MMchol,hₚ,F,_dt,h_wrk)
end

function project_shallow_water_initial_conditions(dΩ, Q, V, S, L2MMchol, RTMMchol, H1MMchol, h₀, u₀, f₀)
  b₁(q)   = ∫(q*h₀)dΩ
  rhs1    = assemble_vector(b₁, Q)
  hn      = FEFunction(Q, copy(rhs1))

  b₂(v)   = ∫(v⋅u₀)dΩ
  rhs2    = assemble_vector(b₂, V)
  un      = FEFunction(V, copy(rhs2))

  b₃(s)   = ∫(s*f₀)*dΩ
  rhs3    = assemble_vector(b₃, S)
  f       = FEFunction(S, copy(rhs3))

  hnv,unv,fv=get_free_dof_values(hn,un,f)

  solve!(hnv, L2MMchol, hnv)
  solve!(unv, RTMMchol, unv)
  solve!(fv , H1MMchol, fv)

  hn, un, f, hnv, unv, fv
end

function shallow_water_explicit_time_stepper(
        model, order, degree,
        h₀, u₀, f₀, g,
        dt, τ, N;
        mass_matrix_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        write_diagnostics=true,
        write_diagnostics_freq=1,
        dump_diagnostics_on_screen=true,
        write_solution=false,
        write_solution_freq=N/10,
        output_dir="nswe_eq_ncells_$(num_cells(model))_order_$(order)_explicit")

  # Forward integration of the shallow water equations
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  dω = Measure(Ω, degree, ReferenceDomain())

  # Setup the trial and test spaces
  R, S, U, V, P, Q = setup_mixed_spaces(model, order)

  # assemble the mass matrices
  H1MM, _, L2MM, H1MMchol, RTMMchol, L2MMchol =
    setup_and_factorize_mass_matrices(dΩ, R, S, U, V, P, Q;
                                      mass_matrix_solver=mass_matrix_solver)

  # Project the initial conditions onto the trial spaces
  hn, un, f, hnv, unv, fv =  project_shallow_water_initial_conditions(dΩ, Q, V, S,
                               L2MMchol, RTMMchol, H1MMchol, h₀, u₀, f₀)

  # work arrays
  h_tmp = copy(hnv)
  w_tmp = copy(fv)
  u_tmp = copy(unv)
  # build the potential vorticity lhs operator once just to initialise
  bmm(a,b) = ∫(a*hn*b)dΩ
  H1h      = assemble_matrix(bmm, R, S)
  H1hchol  = numerical_setup(symbolic_setup(mass_matrix_solver,H1h),H1h)

  function run_simulation(pvd=nothing)
    diagnostics_file = joinpath(output_dir,"nswe_diagnostics.csv")

    hm1    = clone_fe_function(Q,hn)
    hm2    = clone_fe_function(Q,hn)
    hp     = clone_fe_function(Q,hn)
    ϕ      = clone_fe_function(Q,hn)

    um1    = clone_fe_function(V,un)
    um2    = clone_fe_function(V,un)
    up     = clone_fe_function(V,un)
    F      = clone_fe_function(V,un)

    wn     = clone_fe_function(S,f)
    q1     = clone_fe_function(S,f)
    q2     = clone_fe_function(S,f)

    # first step, no leap frog integration
    #shallow_water_explicit_time_step!(hn, un, hp, up, ϕ, F, q1, q2, H1h, H1hchol,
    #                                  model, dΩ, dω, V, Q, R, S, f, g, hm1, um1, hm2, um2,
    #                                  RTMMchol, L2MMchol, dt, τ, false, h_tmp, u_tmp, w_tmp)
    #if (write_diagnostics)
    #  initialize_csv(diagnostics_file,"time", "mass", "vorticity", "kinetic", "potential", "power")
    #end
    #if (write_diagnostics && write_diagnostics_freq==1)
    #  compute_diagnostic_vorticity!(wn, dΩ, S, H1MMchol, un, get_normal_vector(Ω),w_tmp)
    #  dump_diagnostics_shallow_water!(h_tmp, w_tmp,
    #                                  model, dΩ, dω, S, L2MM, H1MM,
    #                                  hn, un, wn, ϕ, F, g, 1, dt,
    #                                  diagnostics_file,
    #                                  dump_diagnostics_on_screen)
    #  writevtk(Ω,"output_00001",cellfields=["hn"=>hn,"un"=>un,"wn"=>wn])
    #end

    # subsequent steps, do leap frog integration
    # (now that we have the state at two previous time levels)
    #for istep in 2:N
    #  h_aux = hm2
    #  hm2   = hm1
    #  hm1   = hn
    #  hn    = h_aux
    #  u_aux = um2
    #  um2   = um1
    #  um1   = un
    #  un    = u_aux

    #  shallow_water_explicit_time_step!(hn, un, hp, up, ϕ, F, q1, q2, H1h, H1hchol,
    #                                    model, dΩ, dω, V, Q, R, S, f, g, hm1, um1, hm2, um2,
    #                                    RTMMchol, L2MMchol, dt, τ, true, h_tmp, u_tmp, w_tmp)

    #  if (write_diagnostics && write_diagnostics_freq>0 && mod(istep, write_diagnostics_freq) == 0)
    #    compute_diagnostic_vorticity!(wn, dΩ, S, H1MMchol, un, get_normal_vector(Ω),w_tmp)
    #    dump_diagnostics_shallow_water!(h_tmp, w_tmp,
    #                                    model, dΩ, dω, S, L2MM, H1MM,
    #                                    hn, un, wn, ϕ, F, g, istep, dt,
    #                                    diagnostics_file,
    #                                    dump_diagnostics_on_screen)
    #  end
    #  if (write_solution && write_solution_freq>0 && mod(istep, write_solution_freq) == 0)
    #    compute_diagnostic_vorticity!(wn, dΩ, S, H1MMchol, un, get_normal_vector(Ω),w_tmp)
    #    writevtk(Ω,"output_$(lpad(istep,5,"0"))",cellfields=["hn"=>hn,"un"=>un,"wn"=>wn])
    #  end
    #end
    for istep in 1:N
      shallow_water_ssprk3_time_step!(hn, un, ϕ, F, q1, H1h, H1hchol,
                                      model, dΩ, dω, V, Q, R, S, f, g, hm1, um1, hm2, um2,
                                      RTMMchol, L2MMchol, dt, τ, h_tmp, u_tmp, w_tmp)
    end
    hn, un
  end
  if (write_diagnostics || write_solution)
    rm(output_dir,force=true,recursive=true)
    mkdir(output_dir)
  end
  #if (write_solution)
  #  pvdfile=joinpath(output_dir,"nswe_eq_ncells_$(num_cells(model))_order_$(order)_explicit")
  #  paraview_collection(run_simulation,pvdfile)
  #else
  run_simulation()
  #end
end
