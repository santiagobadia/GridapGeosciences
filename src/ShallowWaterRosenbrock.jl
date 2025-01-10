function assemble_residuals!(duh, dΩ, dω, Y, qAPVM, ϕ, F, n)
  bᵤ(v) = ∫(DIV(v)*ϕ)dω - ∫(qAPVM*(v⋅⟂(F,n)))dΩ 
  bₕ(q) = ∫(-q*DIV(F))dω
  bₕᵤ((v,q)) = bᵤ(v) + bₕ(q)
  Gridap.FESpaces.assemble_vector!(bₕᵤ, duh, Y)
end

function _compute_potential_vorticity!(q,dΩ,R,S,h,u,f,n,solver)
  a(r,s) = ∫(s*h*r)dΩ
  c(s)   = ∫(perp(n,∇(s))⋅(u) + s*f)dΩ
  A      = assemble_matrix(a,R,S)
  Ass    = symbolic_setup(solver,A)
  Ans    = numerical_setup(Ass,A)
  rhs    = assemble_vector(c,S)
  solve!(get_free_dof_values(q),Ans,rhs)
  consistent!(get_free_dof_values(q)) |> wait
end

function shallow_water_rosenbrock_time_step!(
  y₂, ϕ, F, q₁, q₂, duh₁, duh₂, H1h, H1hns, y_wrk,  # in/out args
  model, dΩ, dω, Y, V, Q, R, S, f, g, y₁, y₀,       # in args
  RTMMns, L2MMns, Amat, Bns, Blfns,                 # more in args
  dt, τ, leap_frog, mm_solver, h_wrk, u_wrk, w_wrk, y_wrk_2, topog=nothing)
  # energetically balanced second order rosenbrock shallow water solver
  # reference: eqns (24) and (39) of
  # https://github.com/BOM-Monash-Collaborations/articles/blob/main/energetically_balanced_time_integration/EnergeticallyBalancedTimeIntegration_SW.tex

  n = get_normal_vector(Triangulation(model))
  dt₁ = dt
  if leap_frog
    dt₁ = 2.0*dt
  end

  y₀v, y₁v, y₂v = get_free_dof_values(y₀,y₁,y₂)
  # multifield terms
  u₁, h₁ = y₁

  # 1.1: the mass flux
  compute_mass_flux!(F,dΩ,V,RTMMns,u₁*h₁,u_wrk)
  # 1.2: the bernoulli function
  if topog==nothing
    compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMns,u₁⋅u₁,h₁,g,h_wrk)
  else
    compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMns,u₁⋅u₁,h₁+topog,g,h_wrk)
  end
  # 1.3: the potential vorticity
  _compute_potential_vorticity!(q₁,dΩ,R,S,h₁,u₁,f,n,mm_solver)
  # 1.4: assemble the momentum and continuity equation residuals
  assemble_residuals!(y_wrk_2, dΩ, dω, Y, q₁ - τ*u₁⋅∇(q₁), ϕ, F, n)

  # Solve for du₁, dh₁ over a MultiFieldFESpace
  if leap_frog
    solve!(duh₁, Blfns, y_wrk_2)
  else
    solve!(duh₁, Bns, y_wrk_2)
  end
  consistent!(duh₁) |> wait

  # update
  y₂v   .= y₀v .+ dt₁ .* duh₁
  u₂, h₂ = y₂

  # 2.1: the mass flux
  #compute_mass_flux!(F,dΩ,V,RTMMns,u₁*(2.0*h₁ + h₂)/6.0+u₂*(h₁ + 2.0*h₂)/6.0,u_wrk)
  compute_mass_flux!(F,dΩ,V,RTMMns,u₂*h₂,u_wrk)
  # 2.2: the bernoulli function
  if topog==nothing
    #compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMns,(u₁⋅u₁ + u₁⋅u₂ + u₂⋅u₂)/3.0,0.5*(h₁ + h₂),g,h_wrk)
    compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMns,u₂⋅u₂,h₂,g,h_wrk)
  else
    #compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMns,(u₁⋅u₁ + u₁⋅u₂ + u₂⋅u₂)/3.0,0.5*(h₁ + h₂)+topog,g,h_wrk)
    compute_bernoulli_potential!(ϕ,dΩ,Q,L2MMns,u₂⋅u₂,h₂+topog,g,h_wrk)
  end
  # 2.3: the potential vorticity
  _compute_potential_vorticity!(q₂,dΩ,R,S,h₂,u₂,f,n,mm_solver)
  # 2.4: assemble the momentum and continuity equation residuals
  #assemble_residuals!(y_wrk_2, dΩ, dω, Y, 0.5*(q₁ - τ*u₁⋅∇(q₁) + q₂ - τ*u₂⋅∇(q₂)), ϕ, F, n)
  assemble_residuals!(y_wrk_2, dΩ, dω, Y, q₂ - τ*u₂⋅∇(q₂), ϕ, F, n)

  # subtract A*[du₁,dh₁] from [du₂,dh₂] vector
  mul!(y_wrk, Amat, duh₁)
  y_wrk_2 .= y_wrk_2 .- y_wrk

  # solve for du₂, dh₂
  solve!(duh₂, Bns, y_wrk_2)
  consistent!(duh₂) |> wait

  # update yⁿ⁺¹
  y₂v .= y₁v .+ dt .* duh₂
end

function compute_mean_depth!(wrk, L2MM, h)
  # compute the mean depth over the sphere, for use in the approximate Jacobian
  mul!(wrk, L2MM, get_free_dof_values(h))
  h_int = sum(wrk)
  a_int = sum(L2MM,dims=[1,2])
  h_avg = h_int/a_int
  h_avg
end

function shallow_water_rosenbrock_time_stepper(
  model, order, degree,
  h₀, u₀, f₀, g, H₀,
  λ, dt, τ, N,
  mass_matrix_solver,
  jacobian_matrix_solver;
  t₀=nothing,
  leap_frog=false,
  write_diagnostics=true,
  write_diagnostics_freq=1,
  dump_diagnostics_on_screen=true,
  write_solution=false,
  write_solution_freq=N/10,
  output_dir="nswe_eq_ncells_$(num_cells(model))_order_$(order)_rosenbrock")

  # Forward integration of the shallow water equations
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  dω = Measure(Ω, degree, ReferenceDomain())

  # Setup the trial and test spaces
  R, S, U, V, P, Q = setup_mixed_spaces(model, order)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  # assemble the mass matrices (RTMM mass matrix not needed)
  H1MM, _, L2MM, H1MMns, RTMMns, L2MMns = setup_and_factorize_mass_matrices(dΩ, R, S, U, V, P, Q;
                                                                            mass_matrix_solver=mass_matrix_solver)

  # Project the initial conditions onto the trial spaces
  b₁(q)   = ∫(q*h₀)dΩ
  b₂(v)   = ∫(v⋅u₀)dΩ
  b₃(s)   = ∫(s*f₀)*dΩ
  rhs1    = assemble_vector(b₃, S)
  f       = FEFunction(S, copy(rhs1))
  w_wrk   = copy(get_free_dof_values(f))
  w_wrk  .= get_free_dof_values(f)
  solve!(get_free_dof_values(f),H1MMns,w_wrk)
  consistent!(get_free_dof_values(f)) |> wait

  # assemble the approximate MultiFieldFESpace Jacobian
  n = get_normal_vector(Ω)

  Amat((u,p),(v,q)) = ∫(-dt*λ*f*(v⋅⟂(u,n)))dΩ + ∫(dt*λ*g*(DIV(v)*p))dω - ∫(dt*λ*H₀*(q*DIV(u)))dω
  Mmat((u,p),(v,q)) = ∫(u⋅v)dΩ + ∫(p*q)dΩ # block mass matrix
  A = assemble_matrix(Amat, X,Y)
  M = assemble_matrix(Mmat, X,Y)
  Bmat((u,p),(v,q)) = ∫(u⋅v)dΩ + ∫(p*q)dΩ + ∫(dt*λ*f*(v⋅⟂(u,n)))dΩ - ∫(dt*λ*g*(DIV(v)*p))dω + ∫(dt*λ*H₀*(q*DIV(u)))dω
  B = assemble_matrix(Bmat, X,Y)

  Bns = numerical_setup(symbolic_setup(jacobian_matrix_solver,B),B)
  Mns = numerical_setup(symbolic_setup(mass_matrix_solver,M),M)
  # leap frog matrices, using 2×dt
  lf = 1.0
  if leap_frog
    lf = 2.0
  end
  Blfmat((u,p),(v,q)) = ∫(u⋅v)dΩ + ∫(p*q)dΩ + ∫(lf*dt*λ*f*(v⋅⟂(u,n)))dΩ - ∫(lf*dt*λ*g*(DIV(v)*p))dω + ∫(lf*dt*λ*H₀*(q*DIV(u)))dω
  Blf   = assemble_matrix(Blfmat, X,Y)
  Blfns = numerical_setup(symbolic_setup(jacobian_matrix_solver,Blf),Blf)

  # multifield initial condtions
  b₄((v,q)) = b₁(q) + b₂(v)
  rhs2      = assemble_vector(b₄, Y)
  y_wrk_2   = copy(rhs2)
  y_wrk_2  .= rhs2
  solve!(rhs2, Mns, y_wrk_2)
  consistent!(rhs2) |> wait
  yn  = FEFunction(Y, rhs2)

  # project the bottom topography onto the L2 space
  topog = nothing
  if t₀ != nothing
    b₅(q) = ∫(q*t₀)dΩ
    rhs3  = assemble_vector(b₅,Q)
    topog = FEFunction(Q, copy(rhs3))
    t_tmp = copy(get_free_dof_values(topog))
    t_tmp .= get_free_dof_values(topog)
    solve!(get_free_dof_values(topog),L2MMns,t_tmp)
    consistent!(get_free_dof_values(topog)) |> wait
  end

  un, hn       = yn
  hnv, fv, ynv = get_free_dof_values(hn,f,yn)

  # work arrays
  h_tmp = copy(hnv)
  w_tmp = copy(fv)
  h_wrk   = copy(hnv)
  u_wrk   = copy(get_free_dof_values(un))

  # build the potential vorticity lhs operator once just to initialise
  bmm(a,b) = ∫(a*hn*b)dΩ
  H1h      = assemble_matrix(bmm, R, S)
  H1hns    = numerical_setup(symbolic_setup(mass_matrix_solver,H1h),H1h)

  function run_simulation()
    diagnostics_file = joinpath(output_dir,"nswe__rosenbrock_diagnostics.csv")

    ϕ      = clone_fe_function(Q,hn)
    F      = clone_fe_function(V,un)
    wn     = clone_fe_function(S,f)
    q1     = clone_fe_function(S,f)
    q2     = clone_fe_function(S,f)

    # mulifield fe functions
    ym1     = clone_fe_function(Y,yn)
    ym2     = clone_fe_function(Y,yn)
    duh1    = copy(ynv)
    duh2    = copy(ynv)
    y_wrk   = copy(ynv)

    if (write_diagnostics)
      initialize_csv(diagnostics_file,"time", "mass", "vorticity", "kinetic", "potential", "power")
    end
    if (write_solution)
      nv = get_normal_vector(Ω)
      compute_diagnostic_vorticity!(wn, dΩ, S, H1MMns, un, nv, w_wrk)
      writevtk(Ω,"output_0000",cellfields=["hn"=>hn,"un"=>un,"wn"=>wn,"nv"=>nv])
    end

    # dummy step
    shallow_water_rosenbrock_time_step!(yn, ϕ, F, q1, q2, duh1, duh2, H1h, H1hns, y_wrk,
                                        model, dΩ, dω, Y, V, Q, R, S, f, g, ym1, ym2,
                                        RTMMns, L2MMns, A, Bns, Blfns, dt, τ, false, 
                                        mass_matrix_solver, h_wrk, u_wrk, w_wrk, y_wrk_2, topog)

    # first step, no leap frog
    istep = 1
    shallow_water_rosenbrock_time_step!(yn, ϕ, F, q1, q2, duh1, duh2, H1h, H1hns, y_wrk,
                                        model, dΩ, dω, Y, V, Q, R, S, f, g, ym1, ym2,
                                        RTMMns, L2MMns, A, Bns, Blfns, dt, τ, false, 
                                        mass_matrix_solver, h_wrk, u_wrk, w_wrk, y_wrk_2, topog)
    if (write_diagnostics && write_diagnostics_freq>0 && mod(istep, write_diagnostics_freq) == 0)
      compute_diagnostic_vorticity!(wn, dΩ, S, H1MMns, un, get_normal_vector(Ω), w_wrk)
      dump_diagnostics_shallow_water!(h_tmp, w_tmp,
                                      model, dΩ, dω, S, L2MM, H1MM,
                                      hn, un, wn, ϕ, F, g, istep, dt,
                                      diagnostics_file,
                                      dump_diagnostics_on_screen)
    end
    if (write_solution && write_solution_freq>0 && mod(istep, write_solution_freq) == 0)
      compute_diagnostic_vorticity!(wn, dΩ, S, H1MMns, un, get_normal_vector(Ω), w_wrk)
      writevtk(Ω,"output_$(lpad(istep,4,"0"))",cellfields=["hn"=>hn,"un"=>un,"wn"=>wn])
    end
    # time step iteration loop
    for istep in 2:N
      aux=ym2
      ym2=ym1
      ym1=yn
      yn=aux
      shallow_water_rosenbrock_time_step!(yn, ϕ, F, q1, q2, duh1, duh2, H1h, H1hns, y_wrk,
                                          model, dΩ, dω, Y, V, Q, R, S, f, g, ym1, ym2,
                                          RTMMns, L2MMns, A, Bns, Blfns, dt, τ, leap_frog, 
                                          mass_matrix_solver, h_wrk, u_wrk, w_wrk, y_wrk_2, topog)

      # IMPORTANT NOTE: We need to extract un, hn out of yn at each iteration because
      #                 the association of yn with its object instance changes at the beginning of
      #                 each iteration
      un, hn = yn
      if (write_diagnostics && write_diagnostics_freq>0 && mod(istep, write_diagnostics_freq) == 0)
        compute_diagnostic_vorticity!(wn, dΩ, S, H1MMns, un, get_normal_vector(Ω), w_wrk)
        dump_diagnostics_shallow_water!(h_tmp, w_tmp,
                                        model, dΩ, dω, S, L2MM, H1MM,
                                        hn, un, wn, ϕ, F, g, istep, dt,
                                        diagnostics_file,
                                        dump_diagnostics_on_screen)
      end
      if (write_solution && write_solution_freq>0 && mod(istep, write_solution_freq) == 0)
        compute_diagnostic_vorticity!(wn, dΩ, S, H1MMns, un, get_normal_vector(Ω), w_wrk)
        writevtk(Ω,"output_$(lpad(istep,4,"0"))",cellfields=["hn"=>hn,"un"=>un,"wn"=>wn])
      end
    end
    hn, un
  end
  if (write_diagnostics || write_solution)
    rm(output_dir,force=true,recursive=true)
    mkdir(output_dir)
  end
  run_simulation()
end
