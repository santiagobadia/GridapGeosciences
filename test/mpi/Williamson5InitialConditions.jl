# Initial conditions for the zonal flow over an isolated mountain ("Williamson 5")
# shallow water test case. Involves a balanced geostrophic flow with a Gaussian
# hill in the northern hemisphere.
# reference:
#   D.L. Williamson, J.B. Drake, J.J. Hack, R. Jakob, P.N. Swarztrauber, 
#   J. Comput. Phys. 102 (1992) 211–224.

const Rₑ = 6371220.0
const gₑ = 9.80616
const U₀ = 20.0
const Ωₑ = 7.292e-5
const H₀ = 5960.0
const α₀ = 0.0

function uθ(θϕr)
  θ,ϕ,r = θϕr
  U₀*(cos(ϕ)*cos(α₀) + cos(θ)*sin(ϕ)*sin(α₀))
end

function uϕ(θϕr)
  θ,ϕ,r = θϕr
  -U₀*sin(θ)*sin(α₀);
end

# Initial velocity
function u₀(xyz)
  θϕr = xyz2θϕr(xyz)
  u   = uθ(θϕr)
  v   = uϕ(θϕr)
  spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)
end

# Topography
function topography(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  #θc    = -π/2.0
  θc    =  0.0
  ϕc    =  π/6.0
  bo    = 2000.0
  rad   = π/9.0
  rsq   = (ϕ - ϕc)*(ϕ - ϕc) + (θ - θc)*(θ - θc)
  r     = sqrt(rsq)
  b     = 0.0
  if(r < rad) 
    b = bo*(1.0 - r/rad)
  end
  b
end

# Initial fluid depth
function h₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  b     = -cos(θ)*cos(ϕ)*sin(α₀) + sin(ϕ)*cos(α₀)
  bt    = topography(xyz)
  H₀ - (Rₑ*Ωₑ*U₀ + 0.5*U₀*U₀)*b*b/gₑ - bt
end
